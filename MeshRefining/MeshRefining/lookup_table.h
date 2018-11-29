#ifndef LOOKUP_TABLE_H
#define LOOKUP_TABLE_H

#include <vector>
#include "mpi.h"
#include "dataclass.h"

using std::vector;

struct MPITypeAssistant;

template <int t_nIndexNodes>
struct LookupTable
{
	static const int N_NODES = t_nIndexNodes + 1;
	static MPI_Datatype s_mpiType;
	friend class MPITypeAssistant;
	static LookupTable<t_nIndexNodes> *&ptrRef(idx_t *&buf) {
		return reinterpret_cast<LookupTable<t_nIndexNodes>*&>(buf);
	}
	static const LookupTable<t_nIndexNodes> *constPtr(const idx_t *buf) {
		return reinterpret_cast<const LookupTable<t_nIndexNodes>*>(buf);
	}
	static void packLookupTable(
			const vector<LookupTable> lookupTable[],
			int *&sendCounts, int *&sendOffsets, idx_t *&sendBuf,
			const int comm_size) {
		packLookupTable(lookupTable, sendCounts, sendOffsets, ptrRef(sendBuf), comm_size);
	}
	static int sendDataAllToAll(
			const idx_t sendBuf[], const int sendCounts[], const int sendOffsets[],
			idx_t *&recvBuf, int *&recvCounts, int *&recvOffsets,
			MPI_Comm comm) {
		return sendDataAllToAll(
					constPtr(sendBuf), sendCounts, sendOffsets,
					ptrRef(recvBuf), recvCounts, recvOffsets, comm);
	}
#ifdef DEBUG
	static void printArray(const int buf[], const int count, const int rank) {
		printArray(constPtr(buf), count, rank);
	}
#endif
	// static int getSize(int count) { return count * N_NODES; }
	// static int getCount(int size) { return size / N_NODES; }

	LookupTable() {
		for (int i = 0; i < N_NODES; ++i)
			allNodes[i] = -1;
	}
	LookupTable(const idx_t nodes[], idx_t lookupIndex) {
		for (int i = 0; i < t_nIndexNodes; ++i)
			this->nodes[i] = nodes[i];
		this->lookupIndex = lookupIndex;
#ifdef DEBUG
		for (int i = 0; i < N_NODES; ++i)
			if (allNodes[i] < 0)
			{
				char buf[255];
				sprintf(buf, "LookupTable<%d>:: -1 index", t_nIndexNodes);
				throw std::out_of_range(buf);
			}
#endif
	}
	LookupTable& operator=(const LookupTable& other) {
		for (int i = 0; i < N_NODES; ++i)
			allNodes[i] = other.allNodes[i];
		return *this;
	}

#ifdef DEBUG
	bool operator==(const LookupTable& other) const {
		for (int i = 0; i < N_NODES; ++i)
			if (allNodes[i] != other.allNodes[i])
				return false;
		return true;
	}
	bool operator!=(const LookupTable& other) const {
		return !(*this == other);
	}
#endif

	union {
		struct {
			// The global IDs of the original nodes for looking up.
			idx_t nodes[t_nIndexNodes];
			// The global ID of the lookup item.
			idx_t lookupIndex;
		};
		idx_t allNodes[N_NODES];
	};

private:
	static void mpiTypeConstructor() {
		MPI_Type_contiguous(N_NODES, MPI_LONG_LONG_INT, &s_mpiType);
		MPI_Type_commit(&s_mpiType);
	}
	static void mpiTypeDestructor() {
		MPI_Type_free(&s_mpiType);
	}

	// template <int t_nIndexNodes>
	static void packLookupTable(
			const vector<LookupTable> lookupTable[],
			int *&sendCounts,
			int *&sendOffsets,
			// int *&sendBuf,
			LookupTable *&sendBuf,
			const int comm_size) {
		// Prepare the sending data.
		sendCounts = new int[comm_size];
		sendOffsets = new int[comm_size];
		sendOffsets[0] = 0;
		int sendCountSum = 0;
		for (int i = 0; i < comm_size; ++i) {
			sendCountSum += sendCounts[i] =
					// LookupTable<t_nIndexNodes>::getSize(lookupTable[i].size());
					lookupTable[i].size();
			if (i > 0)
				sendOffsets[i] = sendOffsets[i - 1] + sendCounts[i - 1];
		}

		// sendBuf = new int[sendCountSum == 0 ? 1 : sendCountSum];
		// sendBuf = reinterpret_cast<LookupTable*>(new int[sendCountSum * N_NODES]);
		sendBuf = new LookupTable[sendCountSum];
		int index = 0;
		for (int i = 0; i < comm_size; ++i)
			// for (int j = 0; j < lookupTable[i].size(); ++j)
			for (int j = 0; j < sendCounts[i]; ++j)
				// for (int k = 0; k < LookupTable<t_nIndexNodes>::N_NODES; ++k)
				// sendBuf[index++] = lookupTable[i][j].allNodes[k];
				sendBuf[index++] = lookupTable[i][j];

#ifdef DEBUG
		assert(index == sendCountSum);
#endif

		// printf("#?: packed node lookups, sum of sendCounts = %d\n", sendCountSum);
	}
	// template <int t_nIndexNodes>
	static int sendDataAllToAll(
			// const int *sendBuf,
			const LookupTable sendBuf[],
			const int sendCounts[],
			const int sendOffsets[],
			// int *&recvBuf,
			LookupTable *&recvBuf,
			int *&recvCounts,
			int *&recvOffsets,
			MPI_Comm comm) {
		int rank;
		int comm_size;
		MPI_Comm_rank(comm, &rank);
		MPI_Comm_size(comm, &comm_size);

#ifdef DEBUG
		for (int r = 0; r < comm_size; ++r) {
			MPI_Barrier(comm);
			if (rank == r)
			{
				printf("#%d: allToAll: sendCounts = {", rank);
				for (int i = 0; i < comm_size; ++i)
					printf("%d, ", sendCounts[i]);
				printf("}\n");
				// fflush(stdout);
			}
		}
#endif

		recvCounts = new int[comm_size];
		MPI_Alltoall(sendCounts, 1, MPI_INT, recvCounts, 1, MPI_INT, comm);

#ifdef DEBUG
		for (int r = 0; r < comm_size; ++r) {
			MPI_Barrier(comm);
			if (rank == r) {
				printf("#%d: allToAll: recvCounts = {", rank);
				for (int i = 0; i < comm_size; ++i)
					printf("%d, ", recvCounts[i]);
				printf("}\n");
				// fflush(stdout);
			}
		}
#endif

		int recvCountSum = calcRecvOffsets(recvCounts, recvOffsets, comm_size);

		// recvBuf = new int[recvCountSum];
		// recvBuf = reinterpret_cast<LookupTable<t_nIndexNodes>*>(
		// 			new int[recvCountSum * LookupTable<t_nIndexNodes>::N_NODES]);
		recvBuf = new LookupTable[recvCountSum];
		// MPI_Alltoallv(sendBuf, sendCounts, sendOffsets, MPI_INT,
		// 			  recvBuf, recvCounts, recvOffsets, MPI_INT, comm);
		MPI_Alltoallv(sendBuf, sendCounts, sendOffsets, s_mpiType,
					  recvBuf, recvCounts, recvOffsets, s_mpiType, comm);

		return recvCountSum;
	}
	static int calcRecvOffsets(
			const int recvCounts[], int *&recvOffsets, int comm_size) {
		recvOffsets = new int[comm_size];
		recvOffsets[0] = 0;
		int recvCountSum = 0;
		for (int i = 0; i < comm_size; ++i) {
			recvCountSum += recvCounts[i];
			if (i > 0)
				recvOffsets[i] = recvOffsets[i - 1] + recvCounts[i - 1];
		}
		return recvCountSum;
	}
#ifdef DEBUG
	static void printArray(const LookupTable buf[], const int count, const int rank) {
		printf("{\n");
		for (int i = 0; i < count; ++i) {
			printf("              {");
			for (int j = 0; j < N_NODES; ++j)
				printf("%d, ", buf[i].allNodes[j]);
			printf("},\n");
		}
		printf("}\n");
	}
#endif
};

typedef LookupTable<1> NodeLookup;
typedef LookupTable<2> EdgeLookup;
typedef LookupTable<3> TriLookup;
typedef LookupTable<1> IndexLookup;

struct MPITypeAssistant
{
	MPITypeAssistant() {
		NodeLookup::mpiTypeConstructor();
		EdgeLookup::mpiTypeConstructor();
		TriLookup::mpiTypeConstructor();
        // IndexLookup::mpiTypeConstructor();
	}
	~MPITypeAssistant() {
		NodeLookup::mpiTypeDestructor();
		EdgeLookup::mpiTypeDestructor();
		TriLookup::mpiTypeDestructor();
        // IndexLookup::mpiTypeDestructor();
	}
};

#endif // LOOKUP_TABLE_H
