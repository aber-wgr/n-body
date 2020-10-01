#include "body.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <functional>
#include <iostream>

#include "../misc/mpi_types.h"
#include "../misc/utils.h"

#include <string>
#include <fstream>
#include <sstream>

#include <sys/stat.h>
#include <iso646.h>
#include "readwrite.h"
#include <memory.h>

int BodyManager::AddBody(double pos[3], double vel[3], double m, double w)
{
    int n = localBodies.nextIndex++;
	
    for(int i=0;i<3;i++)
    {
        localBodies.position.push_back(pos[i]);
        localBodies.velocity.push_back(vel[i]);
    }

    localBodies.position.push_back(0.0);
    localBodies.velocity.push_back(0.0);
	
	
	localBodies.mass.push_back(m);
	localBodies.work.push_back(w);

	return n;
}

void BodyManager::RemoveBodies(std::vector<size_t> indicesToDelete)
{
    // a table transform version of the above function, expanded out and implemented per column
	// Should be roughly O(4M)

    if (indicesToDelete.empty())
        return;

	BodySet out;

    size_t s = localBodies.position.size() - indicesToDelete.size();
    out.position.reserve(s*4);
    out.velocity.reserve(s*4);
    out.mass.reserve(s);
    out.work.reserve(s);
	
	std::sort(indicesToDelete.begin(), indicesToDelete.end());

    // new we can assume there is at least 1 element to delete. copy blocks at a time.
    std::vector<double>::iterator itPosBlockBegin = localBodies.position.begin();
    std::vector<double>::iterator itVelBlockBegin = localBodies.velocity.begin();
    std::vector<double>::iterator itMassBlockBegin = localBodies.mass.begin();
	std::vector<double>::iterator itWorkBlockBegin = localBodies.work.begin();
	
    for (std::vector<size_t>::const_iterator it = indicesToDelete.begin(); it != indicesToDelete.end(); ++it)
    {
        std::vector<double>::iterator itPosBlockEnd = localBodies.position.begin() + (*it * 4);
        std::vector<double>::iterator itVelBlockEnd = localBodies.velocity.begin() + (*it * 4);
        std::vector<double>::iterator itMassBlockEnd = localBodies.mass.begin() + *it;
        std::vector<double>::iterator itWorkBlockEnd = localBodies.work.begin() + *it;
    	
        if (itPosBlockBegin != itPosBlockEnd)
        {
            std::copy(itPosBlockBegin, itPosBlockEnd, std::back_inserter(out.position));
            std::copy(itVelBlockBegin, itVelBlockEnd, std::back_inserter(out.velocity));
            std::copy(itMassBlockBegin, itMassBlockEnd, std::back_inserter(out.mass));
            std::copy(itWorkBlockBegin, itWorkBlockEnd, std::back_inserter(out.work));
        }
        itPosBlockBegin = itPosBlockEnd + 1;
        itVelBlockBegin = itVelBlockEnd + 1;
        itMassBlockBegin = itMassBlockEnd + 1;
        itWorkBlockBegin = itWorkBlockEnd + 1;
    }

    // copy last block.
    if (itPosBlockBegin != localBodies.position.end())
    {
        std::copy(itPosBlockBegin, localBodies.position.end(), std::back_inserter(out.position));
        std::copy(itVelBlockBegin, localBodies.velocity.end(), std::back_inserter(out.velocity));
        std::copy(itMassBlockBegin, localBodies.mass.end(), std::back_inserter(out.mass));
        std::copy(itWorkBlockBegin, localBodies.work.end(), std::back_inserter(out.work));
    }

    localBodies = out;
}

void BodyManager::AddBodies(std::vector<double> positions, std::vector<double> velocities, std::vector<double> masses, std::vector<double> works)
{
    localBodies.position.insert(localBodies.position.end(), positions.begin(), positions.end());
    localBodies.velocity.insert(localBodies.velocity.end(), velocities.begin(), velocities.end());
    localBodies.mass.insert(localBodies.mass.end(), masses.begin(), masses.end());
	localBodies.work.insert(localBodies.work.end(), works.begin(), works.end());
}

double f_rand(double f_min, double f_max)
{
    double f = (double)rand() / RAND_MAX;
    return f_min + f * (f_max - f_min);
}

void BodyManager::GetLocalMinMax(double* min, double* max)
{
    // This is one example of why I haven't joined the position and velocity values together.
	// Since this only deals with positions, I can maximise cache coherence for this comparison.
	
    // Initialize bounds to infinity
    min[0] = min[1] = min[2] = std::numeric_limits<double>::infinity();
    max[0] = max[1] = max[2] = -std::numeric_limits<double>::infinity();
	
    // Find local min and max bounds
    for (auto it = localBodies.position.begin(); it != localBodies.position.end(); it +=4)
    {
        auto b = &it[0];
        if (b[0] < min[0]) min[0] = b[0];
        if (b[1] < min[1]) min[1] = b[1];
        if (b[2] < min[2]) min[2] = b[2];
        if (b[0] > max[0]) max[0] = b[0];
        if (b[1] > max[1]) max[1] = b[1];
        if (b[2] > max[2]) max[2] = b[2];
    }
}

void BodyManager::GetGlobalMinMax(double* min, double* max)
{
    GetLocalMinMax(min, max);
	
	// Find global min and max bounds
    for(int c = 0; c < 3; c++){
        double global_min;
        double global_max;
        MPI_Allreduce(&min[c], &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&max[c], &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        min[c] = global_min;
        max[c] = global_max;
    }
}

double BodyManager::WeightFrac(double split, int coord, const MPI_Comm& comm) {

    int i;
    double work_above, work_below, tot_work_above, tot_work_below;

    // Compute local contribution to work

	
    work_above = work_below = 0;
	
    for (i = 0; i < localBodies.work.size(); i++) 
    {
        if (localBodies.position.at((i*4)+coord) > split) 
        {
            work_above += localBodies.work[i];
        }
        else 
        {
            work_below += localBodies.work[i];
        }
    }
    
    // Sum up over all processes
    MPI_Allreduce(&work_above, &tot_work_above, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&work_below, &tot_work_below, 1, MPI_DOUBLE, MPI_SUM, comm);

    // Return fraction of work below
    return tot_work_below / (tot_work_above + tot_work_below);
}

bool BodyManager::IsAboveSplit(int rank, int n_processes)
{
    int rel_bit = (int) log2(n_processes);
    bool above = (rank >> (rel_bit - 1)) & 1;
    // gray code
    //above = above ^ ((rank >> (rel_bit)) & 1); 
    return above;
}

int BodyManager::GetPartner(int rank, int n_processes)
{
    int rel_bit = (int) log2(n_processes);
    int partner = rank ^ (1 << (rel_bit - 1));
    return partner; 
}


int BodyManager::SplitCoord(double* min, double* max)
{
    double max_diff = -1;
    int coord;
    for (int c = 0; c < 3; c++) {
        // split along maximum dimension
        if (max[c] - min[c] > max_diff) {
            coord = c;
            max_diff = max[c] - min[c];
        }
    }
    return coord;
}

void BodyManager::Orb(std::vector<Bounds>& bounds, std::vector<Bounds>& other_bounds, std::vector<std::pair<int, bool> >& partners, const double* global_min, const double* global_max, int rank, int n_processors)
{
    DebugOutput("ENTERING ORB FUNCTION", rank);
    // Declare variables
    int coord, partner, n_recv_bodies, color, this_side_size;
    int n_splits, n_proc_left;
    bool above;
    double split;
    BodySet other_side, this_side;
    std::function<double(double)> f;
    MPI_Status status;
    
    array<double, 3> min_bounds, max_bounds, min, max, other_min, other_max;

    // Copy over global min and max
    std::copy(global_min, &global_min[3], std::begin(min));
    std::copy(global_max, &global_max[3], std::begin(max));

    // ORB
    color = 0;
    above = false;

    // split until number of processes in subset is 1
    n_splits = log2(n_processors);
    for (int i = 0; i < n_splits; i++) 
    {
        MPI_Comm subset_comm;
        DebugOutput("Running split " + std::to_string(i+1) + "/" + std::to_string(n_splits), rank);
        n_proc_left = n_processors / pow(2, i);

        // Find split in current subset of processors

        // color processors such that each subset gets unique color
        // color of processor depends on how many times it has been 'above'
        color = color << 1;
        if (above) {
            color |= 1;
        }
        DebugOutput("Comm color set to " + std::to_string(color), rank);
        // split the world comm to processor subset
        DebugOutput("Just before MPI_Comm_split", rank);
        MPI_Comm_split(MPI_COMM_WORLD, color, rank, &subset_comm);
        DebugOutput("Just after MPI_Comm_split", rank);

        // choose cartesian coordinate to split
        coord = SplitCoord(min.data(), max.data());
        std::string crd = "x";
    	if (coord > 0) crd = "y";
        if (coord > 1) crd = "z";
        DebugOutput("Chosen split on " + crd + "/" + std::to_string(coord), rank);

        // find split through the bisection method
        // which divides workload equally
        f = [&](double split) {return WeightFrac(split, coord, subset_comm) - 0.5; };
        split = bisection(min[coord], max[coord], f);

        DebugOutput("Splitting at " + std::to_string(split), rank);
    	
        // free comm
        DebugOutput("Just before MPI_Comm_free", rank);
        MPI_Comm_free(&subset_comm);
        DebugOutput("Just after MPI_Comm_free", rank);

        // Decide if process is above or below split
        DebugOutput("Just before IsAboveSplit", rank);
        above = IsAboveSplit(rank, n_proc_left);
        DebugOutput("Just after IsAboveSplit", rank);
        std::string ab = above ? "above" : "below";
        DebugOutput("This process is " + ab + " the split", rank);

        // Save bounds
        other_min = min;
        other_max = max;
        if (above) {
            min[coord] = split;
            other_max[coord] = split;
        }
        else {
            max[coord] = split;
            other_min[coord] = split;
        }
        // bounds on this side of split
        Bounds here,there;
        for (int i = 0; i < 3; i++)
        {
            here.min_bounds[i] = min[i];
            here.max_bounds[i] = max[i];
            there.min_bounds[i] = other_min[i];
            there.max_bounds[i] = other_max[i];
        }
        here.min_bounds[3] = 0.0;
        here.max_bounds[3] = 0.0;
        there.min_bounds[3] = 0.0;
        there.max_bounds[3] = 0.0;

        std::string here_b = "(" + std::to_string(here.min_bounds[0]) + "," + std::to_string(here.min_bounds[1]) + "," + std::to_string(here.min_bounds[2]) + "," + std::to_string(here.min_bounds[3]) + ")";
        std::string there_b = "(" + std::to_string(there.min_bounds[0]) + "," + std::to_string(there.min_bounds[1]) + "," + std::to_string(there.min_bounds[2]) + "," + std::to_string(there.min_bounds[3]) + ")";
    	
        DebugOutput("Bounds here " + here_b, rank);
        DebugOutput("Bounds there " + there_b, rank);
    	
        bounds.push_back(here);
        other_bounds.push_back(there);

        // Partition bodies between split
        other_side.position.clear();
        other_side.velocity.clear();
        other_side.mass.clear();
        other_side.work.clear();

        this_side.position.clear();
        this_side.velocity.clear();
        this_side.mass.clear();
        this_side.work.clear();

        for (int i = 0; i < localBodies.mass.size(); i++)
        {
            std::string ppp = "(" + std::to_string(localBodies.position[i*4+0]) + "," + std::to_string(localBodies.position[i * 4 + 1]) + "," + std::to_string(localBodies.position[i * 4 + 2]) + "," + std::to_string(localBodies.position[i * 4 + 3]) + ")";
            DebugOutput("Locating Body #" + std::to_string(i) + " position " + ppp, rank);
            if ((localBodies.position[i*4+coord] - split > 0) == above)
            {
                DebugOutput("is this side",rank);
                for (int j = 0; j < 4; j++)
                {
                    this_side.position.push_back(localBodies.position[i * 4 + j]);
                    this_side.velocity.push_back(localBodies.velocity[i * 4 + j]);
                }
                this_side.mass.push_back(localBodies.mass.at(i));
                this_side.work.push_back(localBodies.work.at(i));
            }
            else
            {
                DebugOutput("is other side", rank);
                for (int j = 0; j < 4; j++)
                {
                    other_side.position.push_back(localBodies.position[i * 4 + j]);
                    other_side.velocity.push_back(localBodies.velocity[i * 4 + j]);
                }
                other_side.mass.push_back(localBodies.mass.at(i));
                other_side.work.push_back(localBodies.work.at(i));
            }
        }

        this_side_size = this_side.mass.size();

        DebugOutput("This side size:" + std::to_string(this_side_size), rank);
        DebugOutput("Other side size:" + std::to_string(other_side.mass.size()), rank);
    	
        // Exchange bodies with other side
        // get communication partner
        partner = GetPartner(rank, n_proc_left);
        partners.push_back(std::make_pair(partner, above));

        // pack the sending array
        //size_t send_size = other_side.position.size() + other_side.velocity.size() + other_side.mass.size() + other_side.work.size();
        //std::vector<double> send = other_side.position;
        //send.resize(send_size);
        //send.insert(send.end(), other_side.velocity.begin(), other_side.velocity.end());
        //send.insert(send.end(), other_side.mass.begin(), other_side.mass.end());
        //send.insert(send.end(), other_side.work.begin(), other_side.work.end());

        // DebugOutput("just before send array build", rank);

        size_t pos_size = other_side.position.size();
    	size_t vel_size = other_side.velocity.size();
        size_t mass_size = other_side.mass.size();
        size_t work_size = other_side.work.size();

        size_t send_size = pos_size + vel_size + mass_size + work_size;

        double* mem_send = (double*)malloc(send_size * sizeof(double));
        std::fill(mem_send, mem_send+(send_size), 1.0);
        
        std::copy(other_side.position.begin(), other_side.position.end(), mem_send);
        std::copy(other_side.velocity.begin(), other_side.velocity.end(), mem_send + pos_size);
        std::copy(other_side.mass.begin(), other_side.mass.end(), mem_send + pos_size + vel_size);
        std::copy(other_side.work.begin(), other_side.work.end(), mem_send + pos_size + vel_size + mass_size);
    	
        DebugOutput("send array built", rank);
    	
        if (above) 
        {
        	// receive then send
            // since we're sending as a flat double array, we know we'll be sending n * 10 doubles (pos * 4, vel * 4, mass, work).

        	DebugOutput("above: receive then send", rank);
            DebugOutput("to partner " + std::to_string(partner), rank);
            MPI_Probe(partner, 0, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &n_recv_bodies); // length of "send"
            DebugOutput("Receive Data Count:" + std::to_string(n_recv_bodies), rank);
            int n_actual_bodies = n_recv_bodies / 10; //actual body count

            std::vector<double> recv(n_recv_bodies);
            MPI_Recv(recv.data(), n_recv_bodies, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD, &status);
            DebugOutput("Receive Done", rank);
            MPI_Send(mem_send, send_size, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD);
            DebugOutput("Send Done", rank);
            auto posEnd = recv.begin() + n_actual_bodies * 4;
            auto velEnd = recv.begin() + n_actual_bodies * 8;
            auto massEnd = recv.begin() + n_actual_bodies * 9;
        	
            this_side.position.insert(this_side.position.end(), recv.begin(), posEnd);
            this_side.velocity.insert(this_side.velocity.end(), posEnd, velEnd);
            this_side.mass.insert(this_side.mass.end(), velEnd, massEnd);
            this_side.work.insert(this_side.work.end(), massEnd, recv.end());  
        }
        else 
        {
            DebugOutput("below: send then receive", rank);
            DebugOutput("to partner " + std::to_string(partner), rank);
        	// send then receive
            MPI_Send(mem_send, send_size, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD);
            DebugOutput("Send Done", rank);
            // probe recv size and resize vector (same caveat as above)
            MPI_Probe(partner, 0, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &n_recv_bodies);
            DebugOutput("Receive Data Count:" + std::to_string(n_recv_bodies), rank);
            int n_actual_bodies = n_recv_bodies / 10;
            std::vector<double> recv(n_recv_bodies);
            MPI_Recv(recv.data(), n_recv_bodies, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD, &status);
            DebugOutput("Receive Done", rank);

            auto posEnd = recv.begin() + n_actual_bodies * 4;
            auto velEnd = recv.begin() + n_actual_bodies * 8;
            auto massEnd = recv.begin() + n_actual_bodies * 9;

            this_side.position.insert(this_side.position.end(), recv.begin(), posEnd);
            this_side.velocity.insert(this_side.velocity.end(), posEnd, velEnd);
            this_side.mass.insert(this_side.mass.end(), velEnd, massEnd);
            this_side.work.insert(this_side.work.end(), massEnd, recv.end());
        }
    	DebugOutput("Data resorted", rank);
        DebugOutput("this side size after appending:" + std::to_string(this_side.mass.size()), rank);
        DebugOutput("this side positions after appending:" + std::to_string(this_side.position.size()), rank);
    	for(int n=0;n<this_side.mass.size();n++)
    	{
            std::string ppos = "Position #" + std::to_string(n) + "(" + std::to_string(this_side.position[n * 4 + 0]) + "," + std::to_string(this_side.position[n * 4 + 1]) + "," + std::to_string(this_side.position[n * 4 + 2]) + "," + std::to_string(this_side.position[n * 4 + 3]) + ")";
            DebugOutput(ppos, rank);
    	}
        DebugOutput("Body sending done", rank);
        
    	free(mem_send);
        //delete(mem_send);
        
        // Update bodies
        localBodies = this_side;
    }

    DebugOutput("LEAVING ORB FUNCTION", rank);
}

void BodyManager::GenerateBodies(int n, std::array<double, 3> min, std::array<double, 3> max, int rank)
{
    srand(100 * rank);

    for (int i = 0; i < n; i++) {
        // pre-unrolled
        double newPos[4];
        newPos[0] = f_rand(min[0], max[0]);
        newPos[1] = f_rand(min[1], max[1]);
        newPos[2] = f_rand(min[2], max[2]);
        newPos[3] = 0.0;

        double newVel[4];
        newVel[0] = 0.0;
        newVel[1] = 0.0;
        newVel[2] = 0.0;
        newVel[3] = 0.0;
        for (int j = 0; j < 4; j++)
        {
            localBodies.position.push_back(newPos[j]);
            localBodies.velocity.push_back(newVel[j]);
        }
        localBodies.mass.push_back(f_rand(0.5, 1));
        localBodies.work.push_back(1.0);
        localBodies.nextIndex++;
    }
}

int BodyManager::ReadBodies(const char* filename, MPI_Comm comm)
{
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Read bodies from file */
    if (rank != 0) {
        MPI_Status status;
        int x;
        MPI_Recv(&x, 1, MPI_INT, rank - 1, 0,
            MPI_COMM_WORLD, &status);
    }

    std::ifstream infile;
    infile.open(filename);

    std::string line;
    int i = 0;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        double x, y, z, vx, vy, vz, m;
        if (!(iss >> x >> y >> z >> vx >> vy >> vz >> m)) { break; } // error
        if ((i % size) == rank) {
            localBodies.position.push_back(x);
            localBodies.position.push_back(y);
            localBodies.position.push_back(z);
            localBodies.position.push_back(0.0);

            localBodies.velocity.push_back(vx);
            localBodies.velocity.push_back(vy);
            localBodies.velocity.push_back(vz);
            localBodies.velocity.push_back(0.0);

            localBodies.mass.push_back(m);

            localBodies.work.push_back(1.0);

            localBodies.nextIndex++;
        }
        i++;
    }
    infile.close();

    if (rank != size - 1) {
        int x = 1;
        MPI_Send(&x, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }

    return i;
}

void BodyManager::WriteBodies(const char* filename, MPI_Comm comm, bool overwrite)
{
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::ofstream myfile;
    if (rank != 0) {
        MPI_Status status;
        int x;
        MPI_Recv(&x, 1, MPI_INT, rank - 1, 0,
            MPI_COMM_WORLD, &status);
        myfile.open(filename, std::ios::app);
    }
    else {
        if (!file_exists(filename) or overwrite) {
            myfile.open(filename);
        }
        else {
            myfile.open(filename, std::ios::app);
        }
    }

    int index;

    for (int i = 0; i < localBodies.mass.size(); i++)
    {
        for (int j = 0; j < 4; j++)
        {
            myfile << localBodies.position[i * 4 + j] << " ";
        }
        myfile << std::endl;
    }

    if (rank != size - 1) {
        myfile.close();
        int x = 1;
        MPI_Send(&x, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }
    else {
        myfile << std::endl;
        myfile.close();
    }
}
