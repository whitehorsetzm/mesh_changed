#include <cassert>
#include <algorithm>
#include <set>
#include <map>
#include <vector>
#include "boundary.h"
#include "lookup_table.h"
#include <iostream>
using std::map;
using std::multimap;
using std::set;
using std::vector;

template <int NoNodes>
struct OrderedNodes
{
	OrderedNodes() {}
	OrderedNodes(const idx_t indices[]) { init(indices); }
	OrderedNodes(const OrderedNodes& other) { init(other.indices); }

	OrderedNodes& operator=(const OrderedNodes& other)
	{
		init(other.indices);
		return *this;
	}

	bool operator<(const OrderedNodes& other) const
	{
		for (int i = 0; i < NoNodes; ++i)
			if (indices[i] != other.indices[i])
				return indices[i] < other.indices[i];
		return false;
	}

	bool operator==(const OrderedNodes& other) const
	{
		for (int i = 0; i < NoNodes; ++i)
			if (indices[i] != other.indices[i])
				return false;
		return true;
	}

	idx_t indices[NoNodes];

private:
	void init(const idx_t indices[])
	{
		for (int i = 0; i < NoNodes; ++i)
			this->indices[i] = indices[i];
		std::sort(this->indices, this->indices + NoNodes);
#ifdef DEBUG
		if (this->indices[0] < 0)
		{
			char buf[255];
			sprintf(buf, "OrderedNodes<%d>:: -1 index", NoNodes);
			throw std::out_of_range(buf);
		}
#endif
	}
};

template <int t_nIndexNodes>
MPI_Datatype LookupTable<t_nIndexNodes>::s_mpiType;

typedef OrderedNodes<2> OrderedEdge;
typedef OrderedNodes<3> OrderedTri;
typedef OrderedNodes<4> OrderedQuad;

static bool unifyCells(HYBRID_MESH& refinedMesh, MPI_Comm comm)
{
	int rank;
	int comm_size;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &comm_size);

	if (comm_size == 1)
	{
		idx_t globalID = 0;
		for (int i = 0; i < refinedMesh.NumTetras; ++i)
			refinedMesh.pTetras[i].index = globalID++;
		return true;
	}

	if (rank == 0)
		printf("=========================\n");
	printf("#%d: unifyCells started.\n", rank);
	fflush(stdout);
	idx_t numCells = refinedMesh.numOfCells();
	printf("#%d: %lld Cells, %lld tetras.\n", rank, numCells, refinedMesh.NumTetras);

    // MPI_Barrier(comm);

	idx_t offset = numCells;
    MPI_Scan(&numCells, &offset, 1, MPI_LONG_LONG_INT, MPI_SUM, comm);
	offset -= numCells;

#ifdef DEBUG
	printf("#%d: MPI_Scan called, offset = %lld.\n", rank, offset);
	fflush(stdout);
#endif

	int localID = 0;
	for (int i = 0; i < refinedMesh.NumTetras; ++i, ++localID)
		refinedMesh.pTetras[i].index = localID + offset;
		// refinedMesh.pTetras[i].index = i + offset;

	/*
	for (int i = 0; i < refinedMesh.NumHexes; ++i, ++localID)
		refinedMesh.pHexes[i].index = localID + offset;
	*/

	// TODO Unify the other kinds of cell.
	printf("#%d: unifyCells finished.\n", rank);
	return true;
}

static void buildEdgeToNode(
		const HYBRID_MESH& originalMesh,
		HYBRID_MESH& refinedMesh,
		map<OrderedEdge, int>& edgeToNode)
{
	for (int i = 0; i < originalMesh.NumTris; ++i)
	{
		const TRI& tri = originalMesh.pTris[i];
		if (tri.iOppoProc != -1)
		{
			// The added nodes of the j-th edge.
			for (int j = 0; j < 3; ++j)
			{
#ifdef DEBUG
				assert(tri.addedNodes[j] >= originalMesh.NumNodes &&
					   tri.addedNodes[j] < refinedMesh.NumNodes);
#endif
				const int iAddedNode = tri.addedNodes[j];
				const int edgeNodes[] = {
					tri.vertices[(j + 1) % 3],
					tri.vertices[(j + 2) % 3]
				};
				const idx_t edgeGlobalIDs[] = {
					originalMesh.nodes[edgeNodes[0]].index,
					originalMesh.nodes[edgeNodes[1]].index
				};
				OrderedEdge edge(edgeGlobalIDs);
				edgeToNode[edge] = iAddedNode;

				const set<int>& node0Procs = refinedMesh.nodes[edgeNodes[0]].procs;
				const set<int>& node1Procs = refinedMesh.nodes[edgeNodes[1]].procs;
				set<int>& addedNodeProcs = refinedMesh.nodes[iAddedNode].procs;
				std::insert_iterator<set<int> > addedNodeProcsIter =
						std::inserter(addedNodeProcs, addedNodeProcs.begin());
				std::set_intersection(node0Procs.begin(), node0Procs.end(),
									  node1Procs.begin(), node1Procs.end(),
									  addedNodeProcsIter);
				// edgeToProcs[edge] = &addedNodeProcs;
			}
		}
	}
}

static void buildEdgeLookups(
		HYBRID_MESH& refinedMesh,
		const map<OrderedEdge, int>& edgeToNode,
		vector<EdgeLookup> edgeLookups[])
{
	for (map<OrderedEdge, int>::const_iterator mapIter = edgeToNode.begin();
		 mapIter != edgeToNode.end(); ++mapIter)
	{
		const OrderedEdge& edge = mapIter->first;
		const int iNode = mapIter->second;
		const set<int>& procs = refinedMesh.nodes[iNode].procs;
		for (set<int>::const_iterator setIter = procs.begin();
			 setIter != procs.end(); ++setIter)
		{
			EdgeLookup lookup(edge.indices, iNode);
			edgeLookups[*setIter].push_back(lookup);
		}
	}
}

#ifdef DEBUG
static void checkLocalEdgeLookups(
		const EdgeLookup edgeLookups[],
		const int nEdgeLookups,
		const int nOriginalNodes,
		const int nRefinedNodes,
		const int rank)
{
	for (int i = 0; i < nEdgeLookups; ++i)
	{
		int iNode = edgeLookups[i].lookupIndex;
		// assert(iNode >= nOriginalNodes && iNode < nRefinedNodes);
		if (iNode < nOriginalNodes || iNode >= nRefinedNodes)
			printf("#%d: invalid edgeLookup: {%lld, %lld, %lld}\n", rank,
				   edgeLookups[i].nodes[0], edgeLookups[i].nodes[1], iNode);
	}
}

static void checkSendRecvEdgeLookups(
		const EdgeLookup sendEdgeLookups[],
		const int nSendEdgeLookups,
		const EdgeLookup recvEdgeLookups[],
		const int nRecvEdgeLookups,
		const int rank)
{
	assert(nSendEdgeLookups == nRecvEdgeLookups);
	for (int i = 0; i < nSendEdgeLookups; ++i)
	{
		if (sendEdgeLookups[i] != recvEdgeLookups[i])
		{
			printf("#%d conflicts found: sendEdgeLookups = {%lld, %lld, %lld}\n",
				   rank,
				   sendEdgeLookups[i].nodes[0],
				   sendEdgeLookups[i].nodes[1],
				   sendEdgeLookups[i].lookupIndex);
			printf("#%d                  recvEdgeLookups = {%lld, %lld, %lld}\n",
				   rank,
				   recvEdgeLookups[i].nodes[0],
				   recvEdgeLookups[i].nodes[1],
				   recvEdgeLookups[i].lookupIndex);
		}
	}
}
#endif

static void buildEdgeToProcs(
		HYBRID_MESH& refinedMesh,
		const EdgeLookup edgeLookups[],
		const int nEdgeLookups[],
		const map<OrderedEdge, int>& edgeToNode,
		map<OrderedEdge, map<int, int> >& edgeToProcs,
		const idx_t nOriginalNodes,
		const int comm_size, const int rank)
{
	for (int i = nOriginalNodes; i < refinedMesh.NumNodes; ++i)
		refinedMesh.nodes[i].procs.clear();
	int currentProcOffset = 0;
	for (int proc = 0, index = 0; proc < comm_size; ++proc)
	{
		if (proc == rank)
		{
			currentProcOffset = index;
			index += nEdgeLookups[proc];
		}
		else for (int i = 0; i < nEdgeLookups[proc]; ++i, ++index)
		{
			OrderedEdge edge(edgeLookups[index].nodes);
			int localID = edgeLookups[index].lookupIndex;
			edgeToProcs[edge][proc] = localID;

			map<OrderedEdge, int>::const_iterator iter = edgeToNode.find(edge);
			if (iter != edgeToNode.end())
				refinedMesh.nodes[iter->second].procs.insert(proc);
		}
	}
}

static void countDuplicatedNodes(
		const HYBRID_MESH& refinedMesh,
		idx_t &nDupOriginalNodes, idx_t &nDupAddedNodes,
		const idx_t nOriginalNodes, const int rank)
{
	nDupOriginalNodes = 0;
	for (int i = 0; i < nOriginalNodes; ++i)
	{
		const Node& node = refinedMesh.nodes[i];
		if (!node.procs.empty() && *node.procs.begin() < rank)
			++nDupOriginalNodes;
	}
	nDupAddedNodes = 0;
	for (int i = nOriginalNodes; i < refinedMesh.NumNodes; ++i)
	{
		const Node& node = refinedMesh.nodes[i];
		if (!node.procs.empty() && *node.procs.begin() < rank)
			++nDupAddedNodes;
	}
}

static idx_t calcNodeIndexOffset(
		const idx_t nUniqueOriginalNodes, const idx_t nUniqueAddedNodes, MPI_Comm comm)
{
	int rank;
	MPI_Comm_rank(comm, &rank);
	// Sum up to get the global offset due to original nodes.
	// The global IDs of the original nodes will not be changed.
	// const int nUniqueOriginalNodes = nOriginalNodes - nDupOriginalNodes;
	idx_t globalOffset = nUniqueOriginalNodes;
	MPI_Allreduce(&nUniqueOriginalNodes, &globalOffset, 1, MPI_LONG_LONG_INT, MPI_SUM, comm);

	// Sum up to get the local offset due to the added nodes.
	// const int nUniqueAddedNodes = nAddedNodes - nDupAddedNodes;
	idx_t offset = nUniqueAddedNodes;
	MPI_Scan(&nUniqueAddedNodes, &offset, 1, MPI_LONG_LONG_INT, MPI_SUM, comm);
	offset -= nUniqueAddedNodes;
	offset += globalOffset;

	printf("#%d: nUniqueOriginalNodes = %lld, nUniqueAddedNodes = %lld\n",
		   rank, nUniqueOriginalNodes, nUniqueAddedNodes);

	printf("#%d: globalOffset = %lld, offset = %lld\n", rank, globalOffset, offset);

	return offset;
}

static void buildNodeLookups(
		HYBRID_MESH& refinedMesh,
		const map<OrderedEdge, map<int, int> >& edgeToProcs,
		const map<OrderedEdge, int>& edgeToNode,
		vector<NodeLookup> nodeLookups[],
		const idx_t nOriginalNodes, const idx_t offset,
		const int rank)
{
#ifdef DEBUG
	int *theoricalRecvCounts = new int[rank];
	for (int i = 0; i < rank; ++i)
		theoricalRecvCounts[i] = 0;
#endif
	int localID = 0;
	for (int i = nOriginalNodes; i < refinedMesh.NumNodes; ++i)
		if (refinedMesh.nodes[i].procs.empty())
			refinedMesh.nodes[i].index = localID++ + offset;

	// Remap the non-duplicated nodes, and get prepared to
	// send the global IDs of those duplicated added nodes.

	// for (int i = 0; i < nEdgeLookups; ++i)
	for (map<OrderedEdge, int>::const_iterator edgeNodeIter = edgeToNode.begin();
		 edgeNodeIter != edgeToNode.end(); ++edgeNodeIter)
	{
		const OrderedEdge& edge = edgeNodeIter->first;
		const int iNode = edgeNodeIter->second;
		// OrderedEdge edge(edgeLookups[i].nodes);
		// int iNode = edgeLookups[i].lookupIndex;
#ifdef DEBUG
		if (iNode < nOriginalNodes || iNode >= refinedMesh.NumNodes)
			printf("#%d: iNode = %d, nOriginalNodes = %d, nRefinedNodes = %d\n",
				   rank, iNode, nOriginalNodes, refinedMesh.NumNodes);
		assert(iNode >= nOriginalNodes && iNode < refinedMesh.NumNodes);
#endif
		const map<int, int>& procs = edgeToProcs.at(edge);
		map<int, int>::const_iterator iter = procs.begin();
#ifdef DEBUG
		if (rank == 5 && iNode == 1297 ||
				rank == 0 && iNode == 1491 ||
				rank == 6 && iNode == 1404)
		{
			printf("#%d: node %d = {", rank, iNode);
			for (int j = 0; j < 8; ++j)
				if (procs.find(j) != procs.end())
					printf("%d, ", j);
			printf("}\n");
		}
#endif
		if (iter->first >= rank)
		{
			const idx_t globalID =
					refinedMesh.nodes[iNode].index = localID++ + offset;
			// Prepare the global IDs of the duplicated added nodes to send.
			// OrderedEdge edge = addedNodeToEdge.at(i);
			// NodeLookup lookup(edge.indices, addedNodeGlobalID);
			while (iter != procs.end())
			{
				NodeLookup lookup(&(static_cast<const idx_t&>(iter->second)), globalID);
				nodeLookups[iter++->first].push_back(lookup);
			}
			// for (int k = rank + 1; k < comm_size; ++k)
			// 	nodeLookups[k].push_back(lookup);
		}
#ifdef DEBUG
		else
		{
			for (; iter != procs.end() && iter->first < rank; ++iter)
			{
				assert(iter->first >= 0 && iter->first < rank);
				++theoricalRecvCounts[iter->first];
			}
		}
#endif
	}
#ifdef DEBUG
	printf("#%d: theorical recvCounts = {", rank);
	for (int i = 0; i < rank; ++i)
		printf("%d, ", theoricalRecvCounts[i]);
	printf("}\n");
	delete [] theoricalRecvCounts;
#endif
}

static void remapDuplicatedNodes(
		HYBRID_MESH& refinedMesh,
		const NodeLookup nodeLookups[],
		const int nNodeLookups)
{
	for (int i = 0; i < nNodeLookups; ++i)
	{
		int localID = nodeLookups[i].nodes[0];
		idx_t globalID = nodeLookups[i].lookupIndex;
		refinedMesh.nodes[localID].index = globalID;
	}
}

#ifdef DEBUG
static void checkNodeIndices(
		const HYBRID_MESH& refinedMesh, const int nUniqueNodes, MPI_Comm comm)
{
	int rank;
	MPI_Comm_rank(comm, &rank);

	idx_t nGlobalNodes = 0;
	MPI_Allreduce(&nUniqueNodes, &nGlobalNodes, 1, MPI_LONG_LONG_INT, MPI_SUM, comm);

	printf("#%d: nUniqueNodes = %d, nGlobalNodes = %lld\n", rank, nUniqueNodes, nGlobalNodes);

	for (int i = 0; i < refinedMesh.NumNodes; ++i)
		if (refinedMesh.nodes[i].index < 0 ||
				refinedMesh.nodes[i].index >= nGlobalNodes)
		{
			printf("#%d: node index %d (%lf, %lf, %lf) out of range: (%lld)!\n",
				   rank, i,
				   refinedMesh.nodes[i].coord.x,
				   refinedMesh.nodes[i].coord.y,
				   refinedMesh.nodes[i].coord.z,
				   refinedMesh.nodes[i].index);
			char buf[255];
			sprintf(buf, "#%d: node index %d (%lf, %lf, %lf) out of range: (%lld)!",
					rank, i,
					refinedMesh.nodes[i].coord.x,
					refinedMesh.nodes[i].coord.y,
					refinedMesh.nodes[i].coord.z,
					refinedMesh.nodes[i].index);
			// throw std::out_of_range(buf);
		}
}
#endif

static bool unifyNodes(const HYBRID_MESH& originalMesh,
					   HYBRID_MESH &refinedMesh, MPI_Comm comm)
{
	int rank;
	int comm_size;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &comm_size);

	if (comm_size == 1)
	{
		for (int i = 0; i < refinedMesh.NumNodes; ++i)
			refinedMesh.nodes[i].index = i;
		return true;
	}

	if (rank == 0)
		printf("=========================\n");
	printf("#%d: unifyNodes started.\n", rank);
	fflush(stdout);

	const idx_t nOriginalNodes = originalMesh.NumNodes;
	const idx_t nAddedNodes = refinedMesh.NumNodes - nOriginalNodes;

	printf("#%d: nOriginalNodes = %lld, nRefinedNodes = %lld\n",
		   rank, nOriginalNodes, refinedMesh.NumNodes);

	map<OrderedEdge, int> edgeToNode;
	buildEdgeToNode(originalMesh, refinedMesh, edgeToNode);

	vector<EdgeLookup> *edgeLookups = new vector<EdgeLookup>[comm_size];
	buildEdgeLookups(refinedMesh, edgeToNode, edgeLookups);

	int *sendCounts = nullptr;
	int *sendOffsets = nullptr;
	idx_t *sendBuf = nullptr;
	EdgeLookup::packLookupTable(
				edgeLookups, sendCounts, sendOffsets, sendBuf, comm_size);

#ifdef DEBUG
	printf("#%d: size of edgeLookups = {", rank);
	for (int i = 0; i < comm_size; ++i)
		printf("%lu, ", edgeLookups[i].size());
	printf("}\n");
#endif

	idx_t *recvBuf = nullptr;
	int *recvCounts = nullptr;
	int *recvOffsets = nullptr;
	int recvCountSum = EdgeLookup::sendDataAllToAll(
				sendBuf, sendCounts, sendOffsets,
				recvBuf, recvCounts, recvOffsets, comm);

	map<OrderedEdge, map<int, int> > edgeToProcs;
	idx_t nDupOriginalNodes = 0;
	idx_t nDupAddedNodes = 0;
	buildEdgeToProcs(refinedMesh, EdgeLookup::constPtr(recvBuf), recvCounts,
					 edgeToNode, edgeToProcs, nOriginalNodes, comm_size, rank);

	countDuplicatedNodes(refinedMesh, nDupOriginalNodes, nDupAddedNodes,
						 nOriginalNodes, rank);

	printf("#%d: nDupOriginalNodes = %lld, nDupAddedNodes = %lld\n",
		   rank, nDupOriginalNodes, nDupAddedNodes);

	idx_t offset = calcNodeIndexOffset(nOriginalNodes - nDupOriginalNodes,
									   nAddedNodes - nDupAddedNodes, comm);

	vector<NodeLookup> *nodeLookups = new vector<NodeLookup>[comm_size];
	buildNodeLookups(refinedMesh, edgeToProcs, edgeToNode,
					 nodeLookups, nOriginalNodes, offset, rank);

	delete [] sendCounts;
	sendCounts = nullptr;
	delete [] sendOffsets;
	sendOffsets = nullptr;
	delete [] sendBuf;
	sendBuf = nullptr;
	NodeLookup::packLookupTable(nodeLookups, sendCounts, sendOffsets, sendBuf, comm_size);

#ifdef DEBUG
	printf("#%d: size of nodeLookups = {", rank);
	for (int i = 0; i < comm_size; ++i)
		printf("%lu, ", nodeLookups[i].size());
	printf("}\n");
#endif

	delete [] recvBuf;
	recvBuf = nullptr;
	delete [] recvCounts;
	recvCounts = nullptr;
	delete [] recvOffsets;
	recvOffsets = nullptr;
	recvCountSum = NodeLookup::sendDataAllToAll(
				sendBuf, sendCounts, sendOffsets,
				recvBuf, recvCounts, recvOffsets, comm);

	remapDuplicatedNodes(refinedMesh, NodeLookup::constPtr(recvBuf), recvCountSum);

#ifdef DEBUG
	const int nUniqueNodes = refinedMesh.NumNodes - nDupOriginalNodes - nDupAddedNodes;
	checkNodeIndices(refinedMesh, nUniqueNodes, comm);
#endif

	delete [] edgeLookups;
	delete [] nodeLookups;
	delete [] sendCounts;
	delete [] sendOffsets;
	delete [] sendBuf;
	delete [] recvCounts;
	delete [] recvOffsets;
	delete [] recvBuf;
	printf("#%d: unifyNodes finished.\n", rank);
	return true;
}

static idx_t countSurfFacets(HYBRID_MESH& refinedMesh, MPI_Comm comm)
{
	int rank;
	MPI_Comm_rank(comm, &rank);

	idx_t nSurfFacets = 0;
	for (int i = 0; i < refinedMesh.NumTris; ++i)
		if (refinedMesh.pTris[i].iSurf >= 0)
			++nSurfFacets;

	printf("#%d: nSurfFacets = %lld\n", rank, nSurfFacets);

	idx_t offset = nSurfFacets;
	// idx_t nGlobalSurfFacets = nSurfFacets;
	MPI_Scan(&nSurfFacets, &offset, 1, MPI_LONG_LONG_INT, MPI_SUM, comm);
	MPI_Allreduce(&nSurfFacets, &refinedMesh.NumUniqueSurfFacets, 1,
				  MPI_LONG_LONG_INT, MPI_SUM, comm);
	offset -= nSurfFacets;

	int localID = 0;
	for (int i = 0; i < refinedMesh.NumTris; ++i)
		if (refinedMesh.pTris[i].iSurf >= 0)
			refinedMesh.pTris[i].index = localID++ + offset;
	return nSurfFacets;
}

static idx_t countDuplicatedFacets(const HYBRID_MESH& refinedMesh, int rank)
{
	idx_t nDuplicatedFacets = 0;
	for (int i = 0; i < refinedMesh.NumTris; ++i)
	{
		const TRI& tri = refinedMesh.pTris[i];
		if (tri.iSurf < 0 && tri.iOppoProc != -1 && tri.iOppoProc < rank)
			++nDuplicatedFacets;
	}

	for (int i = 0; i < refinedMesh.NumQuads; ++i)
	{
		// TODO: cope with the tetrahedrons
	}

	return nDuplicatedFacets;
}

static void buildFacetLookups(
		const HYBRID_MESH &refinedMesh,
		vector<TriLookup> triLookups[],
		const idx_t offset, const int rank)
{
#ifdef DEBUG
	printf("#%d: buildFacetLookups started.\n", rank);
#endif
	for (int i = 0, localID = 0; i < refinedMesh.NumTris; ++i)
	{
		const TRI& tri = refinedMesh.pTris[i];
		if (tri.iSurf >= 0)
			continue;
		if (tri.iOppoProc == -1)
			refinedMesh.pTris[i].index = localID++ + offset;
		else if (tri.iOppoProc >= rank)
		{
			refinedMesh.pTris[i].index = localID++ + offset;
			idx_t triGlobalID = refinedMesh.pTris[i].index;
			idx_t nodeGlobalIDs[3];
#ifdef DEBUG
			for (int j = 0; j < 3; ++j)
			{
				if (tri.vertices[j] < 0 || tri.vertices[j] >= refinedMesh.NumNodes)
					throw std::out_of_range("buildFacetLookups: triangle indices error!");
				if (refinedMesh.nodes[tri.vertices[j]].index < 0)
					throw std::out_of_range("buildFacetLookups: triangle globalID -1 error!");
			}
#endif
			for (int j = 0; j < 3; ++j)
				nodeGlobalIDs[j] = refinedMesh.nodes[tri.vertices[j]].index;
			OrderedTri sorteTri(nodeGlobalIDs);
			triLookups[tri.iOppoProc].push_back(TriLookup(sorteTri.indices, triGlobalID));
		}
	}

	// for (int i = 0; i < refinedMesh.NumQuads; ++i)
		// refinedMesh.pQuads[i].index = i + offset;
#ifdef DEBUG
	printf("#%d: buildFacetLookups finished.\n", rank);
#endif
}

static void remapDuplicatedFacets(
		HYBRID_MESH& refinedMesh,
		// const map<OrderedEdge, int>& edgeToNode,
		// const vector<NodeLookup> nodeLookups[],
		const TriLookup triLookups[],
		const int nTriLookups)
{
	map<OrderedTri, int> triHash;
	// map<OrderedQuad, int> quadHash;

	for (int i = 0; i < refinedMesh.NumTris; ++i)
	{
		const TRI& tri = refinedMesh.pTris[i];
		idx_t nodeGlobalIDs[3];
		for (int j = 0; j < 3; ++j)
			nodeGlobalIDs[j] = refinedMesh.nodes[tri.vertices[j]].index;
		triHash[OrderedTri(nodeGlobalIDs)] = i;
	}

	for (int i = 0; i < nTriLookups; ++i)
	{
		OrderedTri tri(triLookups[i].nodes);
#ifdef DEBUG
		try {
#endif
		int iTri = triHash.at(tri);
		idx_t globalID = triLookups[i].lookupIndex;
		refinedMesh.pTris[iTri].index = globalID;
#ifdef DEBUG
		} catch (const std::exception& e) {
			printf("#?: Triangle {%lld, %lld, %lld} not found!\n",
				   tri.indices[0], tri.indices[1], tri.indices[2]);
			throw;
		}
#endif
	}
}

static bool unifyFacets(HYBRID_MESH &refinedMesh, MPI_Comm comm)
{
	int rank;
	int comm_size;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &comm_size);

	if (comm_size == 0)
	{
		for (int i = 0, localID = 0; i < refinedMesh.NumTris; ++i)
			refinedMesh.pTris[i].index = localID++;
		refinedMesh.NumUniqueSurfFacets = refinedMesh.numOfFacets();
		refinedMesh.NumUniqueInterfFacets = 0;

		printf("##: nGlobalUniqueInterfFacets = %lld, nGlobalUniqueSurfFacets = %lld\n",
			   refinedMesh.NumUniqueInterfFacets, refinedMesh.NumUniqueSurfFacets);

		return true;
	}

	if (rank == 0)
		printf("=========================\n");
	printf("#%d: unifyFacets started.\n", rank);
	fflush(stdout);

	const idx_t nSurfFacets = countSurfFacets(refinedMesh, comm);
	const idx_t nInterfFacets = refinedMesh.numOfFacets() - nSurfFacets;

	const idx_t nDuplicatedFacets = countDuplicatedFacets(refinedMesh, rank);

	printf("#%d: nDuplicatedFacets = %lld\n", rank, nDuplicatedFacets);

	// Calculate the global indices of the refined mesh
	const idx_t nUniqueInterfFacets = nInterfFacets - nDuplicatedFacets;
	idx_t offset = nUniqueInterfFacets;
	MPI_Scan(&nUniqueInterfFacets, &offset, 1, MPI_LONG_LONG_INT, MPI_SUM, comm);
	MPI_Allreduce(&nUniqueInterfFacets, &refinedMesh.NumUniqueInterfFacets,
				  1, MPI_LONG_LONG_INT, MPI_SUM, comm);

	offset -= nUniqueInterfFacets;
	offset += refinedMesh.NumUniqueSurfFacets;

	printf("#%d: nUniqueInterfFacets = %lld, nSurfFacets = %lld, offset = %lld\n",
		   rank, nUniqueInterfFacets, nSurfFacets, offset);
	if (rank == 0)
		printf("##: nGlobalUniqueInterfFacets = %lld, nGlobalUniqueSurfFacets = %lld\n",
			   refinedMesh.NumUniqueInterfFacets, refinedMesh.NumUniqueSurfFacets);

	vector<TriLookup> *triLookups = new vector<TriLookup>[comm_size];
	buildFacetLookups(refinedMesh, triLookups, offset, rank);

	// for (int i = 0; i < refinedMesh.NumTris; ++i)
		// refinedMesh.pTris[i].index = i + offset;

	// for (int i = 0; i < refinedMesh.NumQuads; ++i)
		// refinedMesh.pQuads[i].index = i + offset;

	int *sendCounts = nullptr;
	int *sendOffsets = nullptr;
	idx_t *sendBuf = nullptr;
	TriLookup::packLookupTable(triLookups, sendCounts, sendOffsets,
							   sendBuf, comm_size);

#ifdef DEBUG
	for (int r = 0; r < comm_size; ++r)
	{
		MPI_Barrier(comm);
		if (rank == r)
		{
			printf("#%d: size of triLookups = {", rank);
			for (int i = 0; i < comm_size; ++i)
				printf("%lu, ", triLookups[i].size());
			printf("}\n");
		}
	}
#endif

	idx_t *recvBuf = nullptr;
	int *recvCounts = nullptr;
	int *recvOffsets = nullptr;
	int recvCountSum = TriLookup::sendDataAllToAll(
				sendBuf, sendCounts, sendOffsets,
				recvBuf, recvCounts, recvOffsets, comm);

	remapDuplicatedFacets(refinedMesh, TriLookup::constPtr(recvBuf), recvCountSum);

	delete [] triLookups;
	delete [] sendCounts;
	delete [] sendOffsets;
	delete [] sendBuf;
	delete [] recvCounts;
	delete [] recvOffsets;
	delete [] recvBuf;
	printf("#%d: unifyFacets finished.\n", rank);
	return true;
}

bool unifyBoundaries(
		const HYBRID_MESH &originalMesh, HYBRID_MESH &refinedMesh, MPI_Comm comm)
{
#ifdef DEBUG
	try {
#endif
	return unifyCells(refinedMesh, comm)
			&& unifyNodes(originalMesh, refinedMesh, comm)
			&& unifyFacets(refinedMesh, comm);
#ifdef DEBUG
	} catch (const std::exception& e) {
		int rank;
		MPI_Comm_rank(comm, &rank);
		printf("#%d: error!%s\n", rank, e.what());
		throw;
	}
#endif
}

static void buildFacetGlobalToLocal(
		const HYBRID_MESH meshParts[], const int nParts,
		map<idx_t, std::pair<int, int> >& triGlobalToLocal)
{
	for (int i = 0; i < nParts; ++i)
		for (int j = 0; j < meshParts[i].NumTris; ++j)
		{
			TRI& tri = meshParts[i].pTris[j];
			const idx_t globalID = tri.index;
			triGlobalToLocal[globalID] = std::make_pair(i, j);
		}
}

static void buildFacetIndexLookups(
		const HYBRID_MESH meshParts[], const int nParts,
		vector<IndexLookup> indexLookups[],
		const int rank, const int comm_size)
{
	for (int i = 0; i < nParts; ++i)
		for (int j = 0; j < meshParts[i].NumTris; ++j)
		{
			TRI& tri = meshParts[i].pTris[j];
			if (tri.iOppoProc >= 0 && tri.iOppoProc < comm_size &&
					tri.iOppoProc != rank)
			{
				const idx_t triGlobalID = tri.index;
				// const int newPartID = meshParts[i].pTetras[tri.iCell].partMarker;
				const int newPartID = i * comm_size + rank;
				IndexLookup lookup(&triGlobalID, newPartID);
				indexLookups[tri.iOppoProc].push_back(lookup);
			}
		}
}


static void remapFacetOppoProcs(
		HYBRID_MESH meshParts[], const int nParts,
		const map<idx_t, std::pair<int, int> >& triGlobalToLocal,
		const IndexLookup indexLookups[],
		const int nIndexLookups)
{
	for (int i = 0; i < nIndexLookups; ++i)
	{
		idx_t triGlobalID = indexLookups[i].nodes[0];
		int iOppoProc = indexLookups[i].lookupIndex;
		std::pair<int, int> pair = triGlobalToLocal.at(triGlobalID);
		int iPart = pair.first;
		int iTri = pair.second;
		meshParts[iPart].pTris[iTri].iOppoProc = iOppoProc;
	}
}

bool unifyRepartitionedFacets(
		HYBRID_MESH meshParts[], const int nParts, MPI_Comm comm)
{
	int rank;
	int comm_size;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &comm_size);

	if (rank == 0)
		printf("=========================\n");
	printf("#%d: unifyRepartitionedFacets started.\n", rank);
	fflush(stdout);

	map<idx_t, std::pair<int, int> > triGlobalToLocal;
	buildFacetGlobalToLocal(meshParts, nParts, triGlobalToLocal);

	vector<IndexLookup> *indexLookups = new vector<IndexLookup>[comm_size];
	buildFacetIndexLookups(meshParts, nParts, indexLookups, rank, comm_size);

	int *sendCounts = nullptr;
	int *sendOffsets = nullptr;
	idx_t *sendBuf = nullptr;
	IndexLookup::packLookupTable(indexLookups, sendCounts, sendOffsets,
								 sendBuf, comm_size);

#ifdef DEBUG
	for (int r = 0; r < comm_size; ++r)
	{
		MPI_Barrier(comm);
		if (rank == r)
		{
			printf("#%d: size of indexLookups = {", rank);
			for (int i = 0; i < comm_size; ++i)
				printf("%lu, ", indexLookups[i].size());
			printf("}\n");
		}
	}
#endif


	idx_t *recvBuf = nullptr;
	int *recvCounts = nullptr;
	int *recvOffsets = nullptr;
	int recvCountSum = IndexLookup::sendDataAllToAll(
				sendBuf, sendCounts, sendOffsets,
				recvBuf, recvCounts, recvOffsets, comm);

	remapFacetOppoProcs(meshParts, nParts, triGlobalToLocal,
						IndexLookup::constPtr(recvBuf), recvCountSum);

	delete [] indexLookups;
	delete [] sendCounts;
	delete [] sendOffsets;
	delete [] sendBuf;
	delete [] recvCounts;
	delete [] recvOffsets;
	delete [] recvBuf;

	printf("#%d: unifyRepartitionedFacets finished.\n", rank);
	fflush(stdout);
	return true;
}
