#pragma once
#ifndef _SUBDOMAIN_H_
#define _SUBDOMAIN_H_

#include <vector>
#include "Node.h"
#include "HaloLayer.h"
#include "mpi.h"
#include <omp.h>
#include "typedef.h"
#include "Struct.h"
#include "Boundary.h"
#include "Utility.h"
#include "ComplexGeometry.h"


class SubDomain {

private:
	int myrank_, num_procs_;
	std::shared_ptr<VelocitySet> velSet_;
public:

	Utility Utility_;
	Geometry Geometry_;
	int start_x_;
	int	start_y_;
	int	start_z_;
	int	end_x_;
	int	end_y_;
	int	end_z_;
	int my_Nx_;
	int my_Ny_;
	int my_Nz_;
	

	int halo_size;
	std::vector<int> gridSize_;
	std::vector<int> NewgridSize_;
	int SubdomainSize_;
	std::vector<std::vector<std::vector<std::shared_ptr<Node>>>> lattice_;
	std::vector<std::shared_ptr<Node>> InletBoundaryNode_;
	std::vector<std::shared_ptr<Node>> OutletBoundaryNode_;
	std::vector<std::shared_ptr<Node>> WallBoundaryNode_;


	MPI_Comm cart_comm_;
	int my_up_, my_down_, my_left_, my_right_, my_front_, my_back_;
	std::vector<int>coords_;
	std::vector<int>dims_;
	std::vector<std::pair<std::tuple<int, int, int>, std::tuple<int, int, int>>> local2global_;
	SubDomain()=default;
	SubDomain(std::shared_ptr<VelocitySet> velSet, Geometry &Geometry);
	~SubDomain();
	void Create_cart_comm();
	void CreateNodes();
	void Connectivity();
	void Initialize();
	bool IsInSubDomain(int coord_x, int coord_y, int coord_z);
	bool IsFluidInSubDomain(int coord_x, int coord_y, int coord_z);
	bool IsNotFluidInSubDomain(int coord_x, int coord_y, int coord_z);
	void GetBoundary();
	void WriteLatticeToVTK(const std::string& filename, const std::vector<std::shared_ptr<Node>>& lattice);
	void WriteNotNullLatToVTKFile(const std::string& filename);
	const std::shared_ptr<VelocitySet> GetVelSet();

	friend class Node;
	friend class Simulation;

};
#endif