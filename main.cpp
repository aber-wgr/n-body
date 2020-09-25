/* input/output */
#include <iostream>

/* math */
#include <cmath>

/* srand */
#include <stdlib.h>

/* STL types */
#include <vector>
#include <array>

/* MPI library */
#include <mpi.h>

/* Custom MPI structs */
#include "misc/mpi_types.h"

/* Body class */
#include "misc/body.h"

/* Reading and writing */
#include "misc/readwrite.h"

/* Tree building */
#include "tree/tree.h"
#include "tree/build_tree.h"

/* Input parsing */
#include "misc/inputparser.h"
#include <iso646.h>
#include "misc/utils.h"


int main(int argc, char * argv[]){
    
    int size, rank, tmax, N, nbodies;
    BodyManager* bodyManager = new BodyManager();
    std::vector<pair<double, int> > splits;
    //vector<pair<array<double, 3>, array<double, 3> > > bounds, other_bounds;
    vector<Bounds> bounds, other_bounds;
    vector<pair<int, bool> > partners;
    vector<double> comp_time;
    double dt, min[3], max[3];
    double start_time, stop_time;
    InputParser ip;
    bool overwrite;
    
    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Number of processes must be a power of 2 */
    if(fmod(log2(size), 1) != 0){
        if(rank == 0){
            std::cerr << "Error: Number of processes must be a power of 2!" << std::endl;
        }
        MPI_Finalize();
        return -1;
    }


    if(!ip.parse(argc, argv)){
        if(rank == 0){
            std::cout << "Error: Invalid command line input" << std::endl;
            print_usage(argc, argv);
        }
        return -1; 
    }

    /* Initialize custom MPI structures */
    init_mpi_types();

    if(ip.read_bodies()){
        /* Read bodies from file */
        bodyManager->ReadBodies(ip.in_file().c_str(), MPI_COMM_WORLD);
        //auto p = read_bodies(ip.in_file().c_str(), MPI_COMM_WORLD);
        // bodies = p.first;
        N = bodyManager->localBodies.mass.size();
    }
    else{
        N = ip.n_bodies();
        nbodies = ip.n_bodies() / size;
        if (rank <= ip.n_bodies() % size - 1){
            nbodies++;
        }
        double lim = (double) 10 * N;
        bodyManager->GenerateBodies(nbodies, { {-lim, -lim, -lim} }, { {lim, lim, lim} }, rank);
        // bodies = generate_bodies(nbodies, {{-lim, -lim, -lim}}, {{lim, lim, lim}}, rank);
    }

    /* Write initial positions to file */
    overwrite = true;
    if(ip.write_positions())
    {
        bodyManager->WriteBodies(ip.out_file().c_str(), MPI_COMM_WORLD, overwrite);
        // write_bodies(ip.out_file().c_str(), bodies, MPI_COMM_WORLD, overwrite);
    }

    
    tmax = ip.n_steps(); // number of time steps
    dt = ip.time_step(); // time step
    
    if(ip.clock_run())
    {
        MPI_Barrier(MPI_COMM_WORLD);
        start_time = MPI_Wtime();
    }

    std::vector<double> timings;
    timings.resize(tmax);
	
    for(int t = 0; t < tmax; t++){

        /* Reset variables */
        bounds.clear();
        other_bounds.clear();
        partners.clear();
        
        /* Domain composition and transfer of bodies */
        bodyManager->GetGlobalMinMax(min, max);
        //global_minmax(bodies, min, max);
        DebugOutput("BEFORE ORB", rank);
        bodyManager->Orb(bounds, other_bounds, partners, min, max, rank, size);
    	//orb(bodies, bounds, other_bounds, partners, min, max, rank, size);
        DebugOutput("AFTER ORB", rank);
        
        /* Build the local tree */
        Tree tree(bodyManager, min, max, ip.bh_approx_constant());
        //DebugOutput("TREE CONSTRUCTED,BUILDING", rank);
        build_tree(bodyManager, bounds, other_bounds, partners, tree, rank);
        //DebugOutput("TREE BUILT", rank);
        /* Compute forces */
        std::vector<array<double, 3> > forces;
    	for(int b=0;b<bodyManager->localBodies.mass.size();b++)
    	{
            /* time the computation */
            double compute_time = MPI_Wtime();
            array<double, 3> f = tree.compute_force(b);
            /* update the workload for the body */
            bodyManager->localBodies.work[b] = MPI_Wtime() - compute_time;
            forces.push_back(f);
    	}
        
        
        /* Update positions */
        for (int i = 0; i < bodyManager->localBodies.mass.size(); i++) 
        {
            // Body & b = bodies[i];
        	
            for(int c = 0; c < 3; c++)
            {
                bodyManager->localBodies.velocity[i * 4 + c] = bodyManager->localBodies.velocity[i * 4 + c] + forces[i][c] * dt / bodyManager->localBodies.mass[i];
            	// b.vel[c] = b.vel[c] + forces[i][c] * dt / b.m;
                bodyManager->localBodies.position[i * 4 + c] = bodyManager->localBodies.position[i * 4 + c] + bodyManager->localBodies.velocity[i * 4 + c] * dt;
            	//b.pos[c] = b.pos[c] + b.vel[c] * dt;
            }
        }
        
        /* Output */
        /* Print time step to stdout */
        if(rank == 0 and ip.verbose()){
            std::cout << "\rTime step: " << t + 1 << "/" << tmax;
            if(t == tmax - 1){
                std::cout << std::endl;
            }
        }


    	
        if(ip.sampling_interval() == 1 or (t % ip.sampling_interval() == 0 and t != 0)){

            /* Stop the time */
            if(ip.clock_run()){
                MPI_Barrier(MPI_COMM_WORLD);
                stop_time = MPI_Wtime(); 
            }

            /* Write tree */
            if(rank == 0 and ip.write_tree()){
                write_tree(ip.out_tree_file().c_str(), tree, true, overwrite);
            }
            
            /* Write tree size*/
            if(rank == 0 and ip.write_tree_size()){
                write_to_file(ip.out_tree_size_file().c_str(), tree.size(), overwrite);
            }

            /* Write positions */
            if(ip.write_positions())
            {
                bodyManager->WriteBodies(ip.out_file().c_str(), MPI_COMM_WORLD, false);
                // write_bodies(ip.out_file().c_str(), bodies, MPI_COMM_WORLD, false);
            }

        	if (rank == 0)
				timings.push_back(stop_time - start_time);

            if(overwrite){
                overwrite = false;
            }

            if(ip.clock_run()){
                // start the time
                MPI_Barrier(MPI_COMM_WORLD);
                start_time = MPI_Wtime();        
            }
        }
    }

    /* Finalize */
    free_mpi_types();
    MPI_Finalize();

    /* Write running time */
    if (ip.clock_run()) {
        if (rank == 0)
        {
            write_vector(ip.out_time_file().c_str(), timings, overwrite);
        }
    }
   
    if(ip.write_summary()){
        if(rank == 0){
            write_summary(ip, N, size);
        }
    }
    
    delete bodyManager;
}

