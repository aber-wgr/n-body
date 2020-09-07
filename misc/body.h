
#ifndef _BODY_H_N_BODY_
#define _BODY_H_N_BODY_
#include <vector>
#include <array>
#include "../tree/tree.h"
#include "mpi_types.h"

#include <mpi.h>

/**
typedef struct {
    double pos[3];
    double vel[3];
    double m;
    double work;
} Body;
*/

// specific exception to SoA. Bounds are nearly always accessed in very close proximity (often in the same comparator)
// and therefore should be kept in the same cache lines
struct Bounds
{
    double min_bounds[4]; // alignment again
    double max_bounds[4];
};

// So the above is replaced with SOA and some of the functions are filtered into this new class. We'll also look at some cache line work too.
// (Example - pos and vel are both 24-byte values, which means they will routinely cross a 128 byte boundary. It's not critical but it's a cost nonetheless.)

struct BodySet
{
	std::vector<double> position;
    std::vector<double> velocity;
    std::vector<double> mass;
    std::vector<double> work;
    int nextIndex = 0;
};

class BodyManager
{
	
public:
	BodySet localBodies; // essentially this is where we keep the transform results each time
	
    int AddBody(double pos[3], double vel[3], double m, double work);   // returns index

    void RemoveBodies(std::vector<size_t> indicesToDelete); // returns number actually removed
    void AddBodies(std::vector<double> positions, std::vector<double> velocities, std::vector<double> masses, std::vector<double> works); // returns number added

    void GenerateBodies(int n, std::array<double, 3> min, std::array<double, 3> max, int rank);

	// Bounds utils
    void GetLocalMinMax(double* min, double* max);
    void GetGlobalMinMax(double* min, double* max);
	
	// ORB functions
	bool IsAboveSplit(int rank, int n_processes);
	int GetPartner(int rank, int n_processes);
    double WeightFrac(double split, int coord, const MPI_Comm& comm);
	void Orb(std::vector<Bounds>& bounds, std::vector<Bounds>& other_bounds, std::vector<std::pair<int, bool> > & partners, const double * global_min, const double * global_max, int rank, int n_processors);

    void WriteBodies(const char* filename, MPI_Comm comm, bool overwrite);
    void ReadBodies(const char* filename, MPI_Comm comm);
};

#endif // _BODY_H_N_BODY_
