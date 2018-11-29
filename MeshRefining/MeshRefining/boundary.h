#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <vector>
#include "mpi.h"
#include "dataclass.h"
#include "lookup_table.h"

using std::vector;

// bool buildBoundaryFacets(HYBRID_MESH& mesh, MPI_Comm comm);
bool unifyBoundaries(const HYBRID_MESH& originalMesh, HYBRID_MESH& refinedMesh, MPI_Comm comm);
bool unifyRepartitionedFacets(HYBRID_MESH meshParts[], const int nParts, MPI_Comm comm);
#endif // BOUNDARY_H
