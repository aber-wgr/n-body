#include "readwrite.h"

#include <string>
#include <fstream>
#include <sstream>

#include <sys/stat.h>
#include <iso646.h>

bool file_exists(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}


void write_to_file(const char* filename, double x, bool overwrite) {
    std::ofstream myfile;
    if (!file_exists(filename) or overwrite) {
        myfile.open(filename);
    }
    else {
        myfile.open(filename, std::ios::app);
    }
    myfile << std::to_string(x) << std::endl;
    myfile.close();
}

void write_vector(const char * filename, std::vector<double> x, bool overwrite){
    std::ofstream myfile;
    if(!file_exists(filename) or overwrite){
        myfile.open(filename);
    }else{
        myfile.open(filename, std::ios::app);
    }
	for(double d : x)
		myfile << std::to_string(d) << std::endl;
    myfile.close();
}

void write_tree(const char * filename, const Tree & tree, bool fulltree, bool overwrite){
    std::ofstream myfile;
    if(!file_exists(filename) or overwrite){
        myfile.open(filename);
    }else{
        myfile.open(filename, std::ios::app);
    }
    myfile << tree.to_string(fulltree) << std::endl;
    myfile.close();
}

void write_summary(InputParser ip, int N, int P){

    std::ofstream file;
    file.open(ip.out_sum_file());

    file << "number of processes: " << P << std::endl;
    file << "number of bodies: " << N << std::endl;
    if(ip.read_bodies()){
        file << "read bodies from file: " << ip.in_file() << std::endl;
    }
    else{
        file << "read bodies from file: none" << std::endl;
    }
    file << "barnes hut approximaiton constant (theta): " << ip.bh_approx_constant() << std::endl;
    file << "gravitational constant (G): " << ip.grav_constant() << std::endl;
    file << "number of time steps: " << ip.n_steps() << std::endl;
    file << "time step: " << ip.time_step() << std::endl;
    file << "sampling interval: " << ip.sampling_interval() << std::endl;
    if(ip.write_positions()){
        file << "output file with positions: " << ip.out_file() << std::endl;
    }
    else{
        file << "output file with positions: none" << std::endl;
    }
    if(ip.write_tree()){
        file << "output file with tree: " << ip.out_tree_file() << std::endl;
    }
    else{
        file << "output file with tree: none" << std::endl;
    }
    if(ip.clock_run()){
        file << "output file with running times: " << ip.out_time_file() << std::endl;
    }
    else{
        file << "output file with running_times: none" << std::endl;
    }

}
