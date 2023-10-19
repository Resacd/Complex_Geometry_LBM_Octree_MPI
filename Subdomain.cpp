#include "Subdomain.h"

SubDomain::SubDomain(std::shared_ptr<VelocitySet> velSet, Geometry& Geometry) {

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank_);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs_);
    Geometry_ = Geometry;
    velSet_ = velSet;
    gridSize_ = Geometry_.GeomSize;
    Create_cart_comm();
    CreateNodes();
    Connectivity();
    Initialize();
    GetBoundary();
    MPI_Comm_free(&cart_comm_);

    //WriteNotNullLatToVTKFile("Lattice.vtk");
    //WriteLatticeToVTK("inletboundary.vtk", InletBoundaryNode_);
    //WriteLatticeToVTK("outletboundary.vtk", OutletBoundaryNode_);
    //WriteLatticeToVTK("wallboundary.vtk", WallBoundaryNode_);




   
}

SubDomain::~SubDomain() {

}

void SubDomain::Create_cart_comm() {
    
    int const ndims = 3;
    dims_ = { 0, 0, 0 };
    int periods[ndims] = { 0, 0, 0 };
    int reorder = 1;
    coords_ = { 0, 0, 0 };
    int remain_dims[3] = { 1, 1, 1 };
    int Nx = gridSize_[0];
    int Ny = gridSize_[1];
    int Nz = gridSize_[2];
    MPI_Dims_create(num_procs_, ndims, dims_.data());
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims_.data(), periods, reorder, &cart_comm_);
    MPI_Cart_shift(cart_comm_, 0, 1, &my_left_, &my_right_);
    MPI_Cart_shift(cart_comm_, 1, 1, &my_down_, &my_up_);
    MPI_Cart_shift(cart_comm_, 2, 1, &my_front_, &my_back_);
    MPI_Comm_rank(cart_comm_, &myrank_);
    MPI_Cart_coords(cart_comm_, myrank_, ndims, coords_.data());
    int cart_coords[ndims];
    MPI_Cart_get(cart_comm_, ndims, dims_.data(), periods, cart_coords);

        halo_size = 1;
        my_Nx_ = Nx / dims_[0];
        my_Ny_ = Ny / dims_[1];
        my_Nz_ = Nz / dims_[2];
        start_x_ = coords_[0] * my_Nx_;
        start_y_ = coords_[1] * my_Ny_;
        start_z_ = coords_[2] * my_Nz_;
        end_x_ = start_x_ + my_Nx_ - 1;
        end_y_ = start_y_ + my_Ny_ - 1;
        end_z_ = start_z_ + my_Nz_ - 1;

        if (Nx % dims_[0] > 0 && coords_[0] == dims_[0] - 1) {
            my_Nx_ += Nx % dims_[0];
        }
        if (Ny % dims_[1] > 0 && coords_[1] == dims_[1] - 1) {
            my_Ny_ += Ny % dims_[1];
        }
        if (Nz % dims_[2] > 0 && coords_[2] == dims_[2] - 1) {
            my_Nz_ += Nz % dims_[2];
        }

        if (dims_[0] > 1) {
            if (coords_[0] == 0) {
                my_Nx_ += halo_size;
                end_x_ += halo_size;
            }
            else if (coords_[0] == dims_[0] - 1) {
                my_Nx_ += halo_size;
                start_x_ += 2 * halo_size * coords_[0] - halo_size;
                end_x_ = start_x_ + my_Nx_ - 1;
            }
            else {
                my_Nx_ += 2 * halo_size;
                start_x_ += coords_[0] * 2 * halo_size - halo_size;
                end_x_ = start_x_ + my_Nx_ - 1;
            }
        }

        if (dims_[1] > 1) {
            if (coords_[1] == 0) {
                my_Ny_ += halo_size;
                end_y_ += halo_size;
            }
            else if (coords_[1] == dims_[1] - 1) {
                my_Ny_ += halo_size;
                start_y_ += 2 * halo_size * coords_[1] - halo_size;
                end_y_ = start_y_ + my_Ny_ - 1;
            }
            else {
                my_Ny_ += 2 * halo_size;
                start_y_ += coords_[1] * 2 * halo_size - halo_size;
                end_y_ = start_y_ + my_Ny_ - 1;
            }
        }

        if (dims_[2] > 1) {
            if (coords_[2] == 0) {
                my_Nz_ += halo_size;
                end_z_ += halo_size;
            }
            else if (coords_[2] == dims_[2] - 1) {
                my_Nz_ += halo_size;
                start_z_ += 2 * halo_size * coords_[2] - halo_size;
                end_z_ = start_z_ + my_Nz_ - 1;
            }
            else {
                my_Nz_ += 2 * halo_size;
                start_z_ += coords_[2] * 2 * halo_size - halo_size;
                end_z_ = start_z_ + my_Nz_ - 1;
            }
        }

        NewgridSize_.resize(3);

        NewgridSize_[0] = (dims_[0] - 1) * 2 * halo_size + gridSize_[0];
        NewgridSize_[1] = (dims_[1] - 1) * 2 * halo_size + gridSize_[1];
        NewgridSize_[2] = (dims_[2] - 1) * 2 * halo_size + gridSize_[2];

    for (int p = 0; p < num_procs_; p++) {
        if (myrank_ == p) {

            printf("Processor %d:\n", myrank_);
            printf("my_Nx_ = %d, my_Ny_ = %d, my_Nz_ = %d\n", my_Nx_, my_Ny_, my_Nz_);
            printf("start_x_ = %d, end_x_ = %d\n", start_x_, end_x_);
            printf("start_y_ = %d, end_y_ = %d\n", start_y_, end_y_);
            printf("start_z_ = %d, end_z_ = %d\n", start_z_, end_z_);

            for (int k = 0; k < my_Nz_; k++) {
                for (int j = 0; j < my_Ny_; j++) {
                    for (int i = 0; i < my_Nx_; i++) {
                        int global_i = i + start_x_;
                        int global_j = j + start_y_;
                        int global_k = k + start_z_;
                        std::tuple<int, int, int> t1 = std::make_tuple(i, j, k);
                        std::tuple<int, int, int> t2 = std::make_tuple(global_i, global_j, global_k);
                        std::pair<std::tuple<int, int, int>, std::tuple<int, int, int>> p = std::make_pair(t1, t2);
                        local2global_.push_back(p);
                        printf("Process %d [%d,%d,%d]: local index = (%d, %d, %d), global index = (%d, %d, %d)\n",
                            myrank_, coords_[0], coords_[1], coords_[2], i, j, k,
                            global_i, global_j, global_k);
                    }
                }
            }
            fflush(stdout);
        }
        MPI_Barrier(cart_comm_); 
    }
}


bool SubDomain::IsInSubDomain(int coord_x, int coord_y, int coord_z) {

    
    if (coord_x >= 0 && coord_x < my_Nx_
        && coord_y >= 0 && coord_y < my_Ny_
        && coord_z >= 0 && coord_z < my_Nz_)
        return true;
    else
        return false;

}

bool SubDomain::IsNotFluidInSubDomain(int coord_x, int coord_y, int coord_z)
{
    
    int global_coordX = coord_x + start_x_;
    int global_coordY = coord_y + start_y_;
    int global_coordZ = coord_z + start_z_;

    if (coords_[0] > 0) {
        global_coordX -= 2 * coords_[0] * halo_size;
    }
    if (coords_[1] > 0) {
        global_coordY -= 2 * coords_[1] * halo_size;
    }
    if (coords_[2] > 0) {
        global_coordZ -= 2 * coords_[2] * halo_size;
    }

    if (IsInSubDomain(coord_x, coord_y, coord_z) && 
        Geometry_.pointsBox[global_coordX][global_coordY][global_coordZ].Type==NodeType::NotFluid)
        return true;
    else
        return false;
}


bool SubDomain::IsFluidInSubDomain(int coord_x, int coord_y, int coord_z)
{

    int global_coordX = coord_x + start_x_;
    int global_coordY = coord_y + start_y_;
    int global_coordZ = coord_z + start_z_;

    if (coords_[0] > 0) {
        global_coordX -= 2 * coords_[0] * halo_size;
    }
    if (coords_[1] > 0) {
        global_coordY -= 2 * coords_[1] * halo_size;
    }
    if (coords_[2] > 0) {
        global_coordZ -= 2 * coords_[2] * halo_size;
    }

    if (IsInSubDomain(coord_x, coord_y, coord_z) && 
        Geometry_.pointsBox[global_coordX][global_coordY][global_coordZ].Type == NodeType::Fluid)
        return true;

    else
        return false;
}


void SubDomain::CreateNodes() {

    lattice_.resize(my_Nx_);

    for (int i = 0; i < my_Nx_; i++) {
        lattice_[i].resize(my_Ny_);
        for (int j = 0; j < my_Ny_; j++) {
            lattice_[i][j].resize(my_Nz_);

        }
    }


    for (int i = 1; i < my_Nx_-1 ; i++) {

        for (int j = 1; j < my_Ny_-1 ; j++) {

            for (int k = 1; k < my_Nz_-1 ; k++) {

                if (IsFluidInSubDomain(i, j, k)) {
                    lattice_[i][j][k] = std::make_shared<Node>(i, j, k,NodeType::Fluid);
                }
                else {
                    lattice_[i][j][k] = nullptr;
                }

            }
        }
    }


    if (my_left_ != MPI_PROC_NULL) {

        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {
                lattice_[0][j][k] = std::make_shared<Node>(0, j, k,NodeType::Halo);
            }
        }
    }

    if (my_left_ == MPI_PROC_NULL) {

        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {
                if (IsFluidInSubDomain(0, j, k)) {
                    lattice_[0][j][k] = std::make_shared<Node>(0, j, k, NodeType::Fluid);
                }
            }
        }
    }

    if (my_right_ != MPI_PROC_NULL) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {
                lattice_[my_Nx_ - 1][j][k] = std::make_shared<Node>(my_Nx_ - 1, j, k, NodeType::Halo);
            }
        }
    }

    if (my_right_ == MPI_PROC_NULL) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {
                if (IsFluidInSubDomain(my_Nx_ - 1, j, k)) {
                    lattice_[my_Nx_ - 1][j][k] = std::make_shared<Node>(my_Nx_ - 1, j, k, NodeType::Fluid);
                }
            }
        }
    }


    if (my_front_ != MPI_PROC_NULL) {
        for (int i = 0; i < my_Nx_; i++) {
            for (int k = 0; k < my_Nz_; k++) {
                lattice_[i][0][k] = std::make_shared<Node>(i, 0, k, NodeType::Halo);
            }
        }
    }

    if (my_front_ == MPI_PROC_NULL) {
        for (int i = 0; i < my_Nx_; i++) {
            for (int k = 0; k < my_Nz_; k++) {
                if (IsFluidInSubDomain(i, 0, k)) {
                    lattice_[i][0][k] = std::make_shared<Node>(i, 0, k, NodeType::Fluid);
                }
            }
        }
    }

    if (my_back_ != MPI_PROC_NULL) {
        for (int i = 0; i < my_Nx_; i++) {
            for (int k = 0; k < my_Nz_; k++) {
                lattice_[i][my_Ny_ - 1][k] = std::make_shared<Node>(i, my_Ny_ - 1, k, NodeType::Halo);
            }
        }
    }

    if (my_back_ == MPI_PROC_NULL) {
        for (int i = 0; i < my_Nx_; i++) {
            for (int k = 0; k < my_Nz_; k++) {
                if (IsFluidInSubDomain(i, my_Ny_ - 1, k)) {
                    lattice_[i][my_Ny_ - 1][k] = std::make_shared<Node>(i, my_Ny_ - 1, k, NodeType::Fluid);
                }
            }
        }
    }



    if (my_down_ != MPI_PROC_NULL) {

        for (int i = 0; i < my_Nx_; i++) {
            for (int j = 0; j < my_Ny_; j++) {

                lattice_[i][j][0] = std::make_shared<Node>(i, j, 0, NodeType::Halo);
            }
        }
    }

    if (my_down_ == MPI_PROC_NULL) {

        for (int i = 0; i < my_Nx_; i++) {
            for (int j = 0; j < my_Ny_; j++) {
                if (IsFluidInSubDomain(i, j, 0)) {
                    lattice_[i][j][0] = std::make_shared<Node>(i, j, 0, NodeType::Fluid);
                }
            }
        }
    }


    if (my_up_ != MPI_PROC_NULL) {

        for (int i = 0; i < my_Nx_; i++) {
            for (int j = 0; j < my_Ny_; j++) {

                lattice_[i][j][my_Nz_ - 1] = std::make_shared<Node>(i, j, my_Nz_ - 1, NodeType::Halo);
            }
        }
    }

    if (my_up_ == MPI_PROC_NULL) {

        for (int i = 0; i < my_Nx_; i++) {
            for (int j = 0; j < my_Ny_; j++) {
                if (IsFluidInSubDomain(i, j, my_Nz_ - 1)) {
                    lattice_[i][j][my_Nz_ - 1] = std::make_shared<Node>(i, j, my_Nz_ - 1, NodeType::Fluid);
                }
            }
        }
    }


    for (int rank = 0; rank < num_procs_; ++rank) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == myrank_) {

            for (int i = 0; i < my_Nx_; i++) {
                for (int j = 0; j < my_Ny_; j++) {
                    for (int k = 0; k < my_Nz_; k++) {
                        if (lattice_[i][j][k]) {
                          //  std::cout << "Processor " << myrank_ << ": Node at (" << i << ", " << j << ", " << k << ")\n";
                        }
                        else {
                           // std::cout << "Processor " << myrank_ << ": nullptr at (" << i << ", " << j << ", " << k << ")\n";
                        }
                    }
                }
            }
        }
    }

}


void SubDomain::Connectivity()
{
    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {
                if (lattice_[i][j][k] != nullptr) {
                    lattice_[i][j][k]->m_neighbours[0] = nullptr;

                    for (int dir = 1; dir < _nLatNodes; dir++) {
                        int x_neighbor = lattice_[i][j][k]->x_position + velSet_->GetDirection(dir)[0];
                        int y_neighbor = lattice_[i][j][k]->y_position + velSet_->GetDirection(dir)[1];
                        int z_neighbor = lattice_[i][j][k]->z_position + velSet_->GetDirection(dir)[2];
                        if (IsFluidInSubDomain(x_neighbor, y_neighbor, z_neighbor)) {
                            lattice_[i][j][k]->m_neighbours[dir] = lattice_[x_neighbor][y_neighbor][z_neighbor];
                        }
                        else
                            lattice_[i][j][k]->m_neighbours[dir] = nullptr;
                    }
                }
            }
        }
    }
}


void SubDomain::Initialize()
{
    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {

                if (lattice_[i][j][k] != nullptr) {

                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        lattice_[i][j][k]->Initialize(velSet_, dir);
                    }
                }
            }
        }
    }
}



void SubDomain::GetBoundary() {
    std::vector<Point> InletBoundary = Geometry_.InletboundaryPoints;
    std::vector<Point> OutletBoundary = Geometry_.OutletboundaryPoints;
    std::vector<Point> WallBoundary = Geometry_.WallboundaryPoints;

    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {
                int global_coordX = i + start_x_;
                int global_coordY = j + start_y_;
                int global_coordZ = k + start_z_;


                if (coords_[0] > 0) {
                    global_coordX -= 2 * coords_[0] * halo_size;
                }
                if (coords_[1] > 0) {
                    global_coordY -= 2 * coords_[1] * halo_size;
                }
                if (coords_[2] > 0) {
                    global_coordZ -= 2 * coords_[2] * halo_size;
                }

                Point& currentPoint = Geometry_.pointsBox[global_coordX][global_coordY][global_coordZ];
                if (Utility_.isPointInBoundary(currentPoint, InletBoundary)) {
                    InletBoundaryNode_.push_back(lattice_[i][j][k]);
                }

                if (Utility_.isPointInBoundary(currentPoint, OutletBoundary)) {
                    OutletBoundaryNode_.push_back(lattice_[i][j][k]);
                }

                if (Utility_.isPointInBoundary(currentPoint, WallBoundary)) {
                    WallBoundaryNode_.push_back(lattice_[i][j][k]);
                }
            }
        }
    }
}


void SubDomain::WriteLatticeToVTK(const std::string& filename, const std::vector<std::shared_ptr<Node>>& lattice) {
    std::ofstream vtkFile(filename);

    if (!vtkFile.is_open()) {
        std::cerr << "Error: Unable to open VTK file for writing." << std::endl;
        return;
    }

    vtkFile << "# vtk DataFile Version 3.0" << std::endl;
    vtkFile << "Lattice Data" << std::endl;
    vtkFile << "ASCII" << std::endl;
    vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    vtkFile << "POINTS " << lattice.size() << " int" << std::endl;
    for (const auto& nodePtr : lattice) {
        vtkFile << nodePtr->x_position << " " << nodePtr->y_position << " " << nodePtr->z_position << std::endl;
    }

    vtkFile << "CELLS " << lattice.size() << " " << 2 * lattice.size() << std::endl;
    for (size_t i = 0; i < lattice.size(); ++i) {
        vtkFile << "1 " << i << std::endl;
    }


    vtkFile << "CELL_TYPES " << lattice.size() << std::endl;
    for (size_t i = 0; i < lattice.size(); ++i) {
        vtkFile << "1" << std::endl;
    }

    vtkFile.close();
    std::cout << "VTK file '" << filename << "' written successfully." << std::endl;
}


void SubDomain::WriteNotNullLatToVTKFile(const std::string& filename) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return;
    }


    outfile << "# vtk DataFile Version 3.0" << std::endl;
    outfile << "Lattice Nodes" << std::endl;
    outfile << "ASCII" << std::endl;
    outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    int numNonNullNodes = 0;
    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {
                if (lattice_[i][j][k] != nullptr) {
                    numNonNullNodes++;
                }
            }
        }
    }

    outfile << "POINTS " << numNonNullNodes << " int" << std::endl;
    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {
                if (lattice_[i][j][k] != nullptr) {
                    outfile << i << " " << j << " " << k << std::endl;
                }
            }
        }
    }

    outfile << "CELLS " << numNonNullNodes << " " << numNonNullNodes * 2 << std::endl;
    for (int cellId = 0; cellId < numNonNullNodes; cellId++) {
        outfile << "1 " << cellId << std::endl;
    }

    outfile << "CELL_TYPES " << numNonNullNodes << std::endl;
    for (int cellId = 0; cellId < numNonNullNodes; cellId++) {
        outfile << "1" << std::endl;
    }

    outfile.close();
}


const std::shared_ptr<VelocitySet> SubDomain::GetVelSet()
{
    return velSet_;
}



