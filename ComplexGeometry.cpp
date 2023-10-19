#include "ComplexGeometry.h"
#include <queue>
#include <unordered_set>
#include <map>


Geometry::Geometry(const std::string& FileName,std::vector<double>& InOutCoord, double dx, int mdepth){

    stlFileName = FileName;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    this->InOutCoord = InOutCoord;
    spacing = dx;
    maxDepth = mdepth;
    GetBoundingbox();
    DistributeTriangles();
    generatePointsInBoundingBox();
    Writer.writePointsToFile(pointsBox, "BoundingBox.vtk");
    octree.bBox = boundingbox;
    buildOctree(octree, maxDepth);
    labelPoints();
    Writer.writeLabeledPointsToFile(pointsBox, "LabeledPoints.vtk");
    std::vector<Point> boundaryPoints = findBoundaryPoints();
    Writer.writePointsToFile(boundaryPoints, "boundaryPoints.vtk");
    getInOutBP();
    excludeInletOutlet();
}


Geometry::~Geometry() {


}

void Geometry::GetBoundingbox() {

    if (myrank == 0) {
        globalTriangles = Reader.readSTLASCII(stlFileName);
        globalNumtris = globalTriangles.size();
        Point minPoint = { std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() };
        Point maxPoint = { -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity() };

        for (int t = 0; t < globalTriangles.size(); ++t) {
            const auto& triangle = globalTriangles[t];
            {
                minPoint.x = std::min(minPoint.x, std::min(triangle.x1, std::min(triangle.x2, triangle.x3)));
                minPoint.y = std::min(minPoint.y, std::min(triangle.y1, std::min(triangle.y2, triangle.y3)));
                minPoint.z = std::min(minPoint.z, std::min(triangle.z1, std::min(triangle.z2, triangle.z3)));
                maxPoint.x = std::max(maxPoint.x, std::max(triangle.x1, std::max(triangle.x2, triangle.x3)));
                maxPoint.y = std::max(maxPoint.y, std::max(triangle.y1, std::max(triangle.y2, triangle.y3)));
                maxPoint.z = std::max(maxPoint.z, std::max(triangle.z1, std::max(triangle.z2, triangle.z3)));
            }
        }

        boundingbox = { minPoint.x,maxPoint.x, minPoint.y,maxPoint.y, minPoint.z, maxPoint.z };
    }
    MPI_Bcast(&globalNumtris, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Datatype mpiBoundingBoxType = mpi_caller.createBoundingBoxMPIType();
    MPI_Bcast(&boundingbox, 1, mpiBoundingBoxType, 0, MPI_COMM_WORLD);


    //std::cout << std::endl;

   // std::cout << "Processor " << myrank << ": BoundingBox: "
   //    << "Min(" << boundingbox.min_x << ", " << boundingbox.min_y << ", " << boundingbox.min_z << ") "
   //    << "Max(" << boundingbox.max_x << ", " << boundingbox.max_y << ", " << boundingbox.max_z << ")"
   //    << std::endl;

    int nx = std::ceil((boundingbox.max_x - boundingbox.min_x) / spacing);
    int ny = std::ceil((boundingbox.max_y - boundingbox.min_y) / spacing);
    int nz = std::ceil((boundingbox.max_z - boundingbox.min_z) / spacing);

    GeomSize = { nx + 1, ny + 1, nz + 1 };
}

void Geometry::DistributeTriangles() {
    int localNumTriangles = globalNumtris / num_procs;
    int remainingTriangles = globalNumtris % num_procs;

    std::vector<int> sendcounts(num_procs, localNumTriangles);

    for (int i = 0; i < remainingTriangles; ++i) {
        sendcounts[i]++;
    }

    std::vector<int> displacements(num_procs, 0);

    for (int i = 1; i < num_procs; ++i) {
        displacements[i] = displacements[i - 1] + sendcounts[i - 1];
    }

    //std::cout << "Processor " << myrank << ": globaltris.size(): " << globalNumtris << std::endl;
    //std::cout << "Processor " << myrank << ": localNumTriangles: " << sendcounts[myrank] << std::endl;

    localTriangles.resize(sendcounts[myrank]);
    MPI_Datatype mpiTriangleType = mpi_caller.createTriangleMPIType();
    MPI_Scatterv(globalTriangles.data(), sendcounts.data(), displacements.data(), mpiTriangleType,
        localTriangles.data(), sendcounts[myrank], mpiTriangleType, 0, MPI_COMM_WORLD);

    //std::cout << "Processor " << myrank << ": localTriangles.size(): " << localTriangles.size() << std::endl;
}


bool Geometry::isPointInsideBoundingBox(double x, double y, double z, const BoundingBox& box) {
    return (x >= box.min_x && x <= box.max_x &&
        y >= box.min_y && y <= box.max_y &&
        z >= box.min_z && z <= box.max_z);
}

std::vector<Triangle> Geometry::findTrianglesInSubBox(const BoundingBox& subBox) {
    std::vector<Triangle> trianglesInBox;

    for (auto& triangle : localTriangles) {
        bool point1Inside = isPointInsideBoundingBox(triangle.x1, triangle.y1, triangle.z1, subBox);
        bool point2Inside = isPointInsideBoundingBox(triangle.x2, triangle.y2, triangle.z2, subBox);
        bool point3Inside = isPointInsideBoundingBox(triangle.x3, triangle.y3, triangle.z3, subBox);

        if (point1Inside || point2Inside || point3Inside) {
            trianglesInBox.push_back(triangle);
        }
    }

    return trianglesInBox;
}


void Geometry::buildOctree(Octree& octree, int depth) {

    octree.tri = findTrianglesInSubBox(octree.bBox);
   if (depth > 0) {
        double mid_x = 0.5 * (octree.bBox.min_x + octree.bBox.max_x);
        double mid_y = 0.5 * (octree.bBox.min_y + octree.bBox.max_y);
        double mid_z = 0.5 * (octree.bBox.min_z + octree.bBox.max_z);

        octree.children.resize(8);

        octree.children[0].bBox = { octree.bBox.min_x, mid_x, octree.bBox.min_y, mid_y, octree.bBox.min_z, mid_z };
        octree.children[1].bBox = { mid_x, octree.bBox.max_x, octree.bBox.min_y, mid_y, octree.bBox.min_z, mid_z };
        octree.children[2].bBox = { octree.bBox.min_x, mid_x, octree.bBox.min_y, mid_y, mid_z, octree.bBox.max_z };
        octree.children[3].bBox = { mid_x, octree.bBox.max_x, octree.bBox.min_y, mid_y, mid_z, octree.bBox.max_z };
        octree.children[4].bBox = { octree.bBox.min_x, mid_x, mid_y, octree.bBox.max_y, octree.bBox.min_z, mid_z };
        octree.children[5].bBox = { mid_x, octree.bBox.max_x, mid_y, octree.bBox.max_y, octree.bBox.min_z, mid_z };
        octree.children[6].bBox = { octree.bBox.min_x, mid_x, mid_y, octree.bBox.max_y, mid_z, octree.bBox.max_z };
        octree.children[7].bBox = { mid_x, octree.bBox.max_x, mid_y, octree.bBox.max_y, mid_z, octree.bBox.max_z };


#pragma omp parallel for if (depth > 0)
        for (int i = 0; i < 8; ++i) {
            buildOctree(octree.children[i], depth - 1);
        }
    }
}

std::vector<Octree> Geometry::findDeepestIntersectedLeaves(Line& line, Octree& octree) {
    double t1=0, t2=0;
    std::vector<Octree> intersectedLeaves;

    if (!DoesLineBoxIntersect(octree.bBox,line,t1,t2)) {
        return intersectedLeaves;  
    }

    if (octree.children.empty()) {
        intersectedLeaves.push_back(octree);
        return intersectedLeaves; 
    }

    for (auto& child : octree.children) {
        auto childIntersectedLeaves = findDeepestIntersectedLeaves(line, child);
        intersectedLeaves.insert(intersectedLeaves.end(), childIntersectedLeaves.begin(), childIntersectedLeaves.end());
    }

    return intersectedLeaves; 
}

void Geometry::generatePointsInBoundingBox() {

    int nx = std::ceil((boundingbox.max_x - boundingbox.min_x) / spacing);
    int ny = std::ceil((boundingbox.max_y - boundingbox.min_y) / spacing);
    int nz = std::ceil((boundingbox.max_z - boundingbox.min_z) / spacing);

    GeomSize = { nx + 1, ny + 1, nz + 1 };

    pointsBox.resize(nx + 1);
    for (int i = 0; i <= nx; ++i) {
        pointsBox[i].resize(ny + 1);
        for (int j = 0; j <= ny; ++j) {
            pointsBox[i][j].resize(nz + 1);
        }
    }

    double x, y, z;

    for (int k = 0; k <= nz; ++k) {
        for (int j = 0; j <= ny; ++j) {
            for (int i = 0; i <= nx; ++i) {
                x = boundingbox.min_x + i * spacing;
                y = boundingbox.min_y + j * spacing;
                z = boundingbox.min_z + k * spacing;
                Point p{ x, y, z };
                pointsBox[i][j][k] = p;
            }
        }
    }
}

bool Geometry::DoesLineBoxIntersect(BoundingBox& box, Line& line, double t0, double t1) {
    double tmin, tmax, tymin, tymax, tzmin, tzmax;

    int sign[3];
    sign[0] = (line.direction().x < 0);
    sign[1] = (line.direction().y < 0);
    sign[2] = (line.direction().z < 0);

    tmin = (box.min_x - line.start.x) / line.direction().x;
    tmax = (box.max_x - line.start.x) / line.direction().x;
    tymin = (box.min_y - line.start.y) / line.direction().y;
    tymax = (box.max_y - line.start.y) / line.direction().y;

    if ((tmin > tymax) || (tymin > tmax))
        return false;

    if (tymin > tmin)
        tmin = tymin;

    if (tymax < tmax)
        tmax = tymax;

    tzmin = (box.min_z - line.start.z) / line.direction().z;
    tzmax = (box.max_z - line.start.z) / line.direction().z;

    if ((tmin > tzmax) || (tzmin > tmax))
        return false;

    if (tzmin > tmin)
        tmin = tzmin;

    if (tzmax < tmax)
        tmax = tzmax;

    return ((tmin < t1) && (tmax > t0));
}



bool Geometry::doLinesIntersect(const Line& line1, const Line& line2) {

    double x1 = line1.start.x;
    double y1 = line1.start.y;
    double z1 = line1.start.z;
    double x2 = line1.end.x;
    double y2 = line1.end.y;
    double z2 = line1.end.z;

    double x3 = line2.start.x;
    double y3 = line2.start.y;
    double z3 = line2.start.z;
    double x4 = line2.end.x;
    double y4 = line2.end.y;
    double z4 = line2.end.z;

    double ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / ((y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1));
    double ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / ((y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1));

    const double epsilon = 1e-7;
    if (ua >= -epsilon && ua <= 1 + epsilon && ub >= -epsilon && ub <= 1 + epsilon) {
        return true;
    }

    return false;
}


bool Geometry::isPointInsideTriangle(double x, double y, const Triangle& triangle) {

    double x1 = triangle.x1;
    double y1 = triangle.y1;
    double x2 = triangle.x2;
    double y2 = triangle.y2;
    double x3 = triangle.x3;
    double y3 = triangle.y3;

    double d1 = ((x - x2) * (y1 - y2) - (y - y2) * (x1 - x2));
    double d2 = ((x - x3) * (y2 - y3) - (y - y3) * (x2 - x3));
    double d3 = ((x - x1) * (y3 - y1) - (y - y1) * (x3 - x1));

    bool hasNegative = (d1 < 0) || (d2 < 0) || (d3 < 0);
    bool hasPositive = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(hasNegative && hasPositive);
}


int Geometry::countIntersectinsTriangle(Line& line,std::vector<Triangle> &posintersectedtri) {
    const double threshold = 1e-7;
    int count = 0;

        for (const auto& triangle : posintersectedtri) {

            Line edge1, edge2, edge3;
            edge1.start = { triangle.x1, triangle.y1, triangle.z1 };
            edge1.end = { triangle.x2, triangle.y2, triangle.z2 };
            edge2.start = { triangle.x2, triangle.y2, triangle.z2 };
            edge2.end = { triangle.x3, triangle.y3, triangle.z3 };
            edge3.start = { triangle.x3, triangle.y3, triangle.z3 };
            edge3.end = { triangle.x1, triangle.y1, triangle.z1 };

            if (doLinesIntersect(line, edge1) || doLinesIntersect(line, edge2) || doLinesIntersect(line, edge3)) {
                count++;
                continue;
            }

            double px = line.start.x;
            double py = line.start.y;
            double pz = line.start.z;

            double nx = triangle.nx;
            double ny = triangle.ny;
            double nz = triangle.nz;

            double qx = triangle.x1;
            double qy = triangle.y1;
            double qz = triangle.z1;

            double t = ((qx - px) * nx + (qy - py) * ny + (qz - pz) * nz) /
                ((line.end.x - line.start.x) * nx + (line.end.y - line.start.y) * ny + (line.end.z - line.start.z) * nz);

            if (t >= -threshold && t <= 1 + threshold) {
                double ix = px + t * (line.end.x - line.start.x);
                double iy = py + t * (line.end.y - line.start.y);
                double iz = pz + t * (line.end.z - line.start.z);

                if (isPointInsideTriangle(ix, iy, triangle) && t >= -threshold && t <= 1 + threshold) {
                    count++;
                }
            }
        }

    return count;
}

void Geometry::adjustLineCoordinates(Line& line, BoundingBox& boundingbox, double rayOffset, int directionSwitch) {
    switch (directionSwitch) {
    case 0:  // positive z-direction
        if (line.start.x == boundingbox.min_x) {
            line.start.x += rayOffset;
            line.end.x += rayOffset;
        }

        if (line.start.x == boundingbox.max_x) {
            line.start.x -= rayOffset;
            line.end.x -= rayOffset;
        }

        if (line.start.y == boundingbox.min_y) {
            line.start.y += rayOffset;
            line.end.y += rayOffset;
        }

        if (line.start.y == boundingbox.max_y) {
            line.start.y -= rayOffset;
            line.end.y -= rayOffset;
        }

        break;

    case 1:  // negative z-direction
        if (line.start.x == boundingbox.min_x) {
            line.start.x += rayOffset;
            line.end.x += rayOffset;
        }

        if (line.start.x == boundingbox.max_x) {
            line.start.x -= rayOffset;
            line.end.x -= rayOffset;
        }

        if (line.start.y == boundingbox.min_y) {
            line.start.y += rayOffset;
            line.end.y += rayOffset;
        }

        if (line.start.y == boundingbox.max_y) {
            line.start.y -= rayOffset;
            line.end.y -= rayOffset;
        }

        break;

    case 2:  // positive x-direction
    case 3:  // negative x-direction
        if (line.start.z == boundingbox.min_z) {
            line.start.z += rayOffset;
            line.end.z += rayOffset;
        }

        if (line.start.z == boundingbox.max_z) {
            line.start.z -= rayOffset;
            line.end.z -= rayOffset;
        }

        if (line.start.y == boundingbox.min_y) {
            line.start.y += rayOffset;
            line.end.y += rayOffset;
        }

        if (line.start.y == boundingbox.max_y) {
            line.start.y -= rayOffset;
            line.end.y -= rayOffset;
        }
        break;

    case 4:  // positive y-direction
    case 5:  // negative y-direction
        if (line.start.x == boundingbox.min_x) {
            line.start.x += rayOffset;
            line.end.x += rayOffset;
        }

        if (line.start.x == boundingbox.max_x) {
            line.start.x -= rayOffset;
            line.end.x -= rayOffset;
        }

        if (line.start.z == boundingbox.min_z) {
            line.start.z += rayOffset;
            line.end.z += rayOffset;
        }

        if (line.start.z == boundingbox.max_z) {
            line.start.z -= rayOffset;
            line.end.z -= rayOffset;
        }
        break;
    }
}


void Geometry::labelPoints() {
    const double rayOffset = 0.01;
    std::vector<Triangle> possibleIntersectedTriangles;
    int directionSwitch;
    int local_intersected_triangles;
    int total_intersected_triangles;
    int intersectingCount;
    rayLength = boundingbox.max_z - boundingbox.min_z;
    for (int i = 0; i < pointsBox.size(); i++) {
        for (int j = 0; j < pointsBox[i].size(); j++) {
            for (int k = 0; k < pointsBox[i][j].size(); k++) {
                Line line;
                line.start.x = pointsBox[i][j][k].x;
                line.start.y = pointsBox[i][j][k].y;
                line.start.z = pointsBox[i][j][k].z;

                directionSwitch = 0;
                local_intersected_triangles=0;
                total_intersected_triangles=0;
                bool foundNonZeroIntersections = false;

                do {
                    switch (directionSwitch) {
                    case 0:  // positive z-direction
                        rayLength  = boundingbox.max_z - boundingbox.min_z;
                        line.end.x = pointsBox[i][j][k].x;
                        line.end.y = pointsBox[i][j][k].y;
                        line.end.z = pointsBox[i][j][k].z + rayLength + 0.1;
                        break;
                    case 1:  // negative z-direction
                        rayLength = boundingbox.max_z - boundingbox.min_z;
                        line.end.x = pointsBox[i][j][k].x;
                        line.end.y = pointsBox[i][j][k].y;
                        line.end.z = pointsBox[i][j][k].z - (rayLength + 0.1);
                        break;
                    case 2:  // positive x-direction
                        rayLength = boundingbox.max_x - boundingbox.min_x;
                        line.end.x = pointsBox[i][j][k].x + rayLength + 0.1;
                        line.end.y = pointsBox[i][j][k].y;
                        line.end.z = pointsBox[i][j][k].z;
                        break;
                    case 3:  // negative x-direction
                        rayLength = boundingbox.max_x - boundingbox.min_x;
                        line.end.x = pointsBox[i][j][k].x - (rayLength + 0.1);
                        line.end.y = pointsBox[i][j][k].y;
                        line.end.z = pointsBox[i][j][k].z;
                        break;
                    case 4:  // positive y-direction
                        rayLength = boundingbox.max_y - boundingbox.min_y;
                        line.end.x = pointsBox[i][j][k].x;
                        line.end.y = pointsBox[i][j][k].y + rayLength + 0.1;
                        line.end.z = pointsBox[i][j][k].z;
                        break;
                    case 5:  // negative y-direction
                        rayLength = boundingbox.max_y - boundingbox.min_y;
                        line.end.x = pointsBox[i][j][k].x;
                        line.end.y = pointsBox[i][j][k].y - (rayLength + 0.1);
                        line.end.z = pointsBox[i][j][k].z;
                        break;
                    }

                    adjustLineCoordinates(line, boundingbox, rayOffset, directionSwitch);

                    std::vector<Octree> intersectedLeaves = findDeepestIntersectedLeaves(line, octree);
                    possibleIntersectedTriangles.clear();

                    for (auto& leaf : intersectedLeaves) {
                        possibleIntersectedTriangles.insert(
                            possibleIntersectedTriangles.end(),
                            leaf.tri.begin(), leaf.tri.end());
                    }

                    local_intersected_triangles = countIntersectinsTriangle(line,localTriangles);
                    MPI_Allreduce(&local_intersected_triangles, &total_intersected_triangles, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                    if (total_intersected_triangles != 0) {
                        foundNonZeroIntersections = true;
                        break;
                    }
                    else {
                        directionSwitch++;
                        if (directionSwitch > 5) {
                            break;
                        }
                    }
                    
                } while (true);

                if (total_intersected_triangles % 2 == 1) {
                    pointsBox[i][j][k].Type = NodeType::Fluid;
                }
                else {
                    pointsBox[i][j][k].Type = NodeType::NotFluid;
                }
            }
        }
    }
}



std::vector<Point> Geometry::findBoundaryPoints() {
    for (int i = 0; i < pointsBox.size(); i++) {
        for (int j = 0; j < pointsBox[i].size(); j++) {
            for (int k = 0; k < pointsBox[i][j].size(); k++) {
                Point& point = pointsBox[i][j][k];

                if (point.Type == NodeType::Fluid) {
                    bool isBoundaryPoint = false;
                    for (int dx = -1; dx <= 1; dx++) {
                        for (int dy = -1; dy <= 1; dy++) {
                            for (int dz = -1; dz <= 1; dz++) {
                                if (i + dx >= 0 && i + dx < pointsBox.size() &&
                                    j + dy >= 0 && j + dy < pointsBox[i].size() &&
                                    k + dz >= 0 && k + dz < pointsBox[i][j].size()) {
                                    Point& neighbor = pointsBox[i + dx][j + dy][k + dz];
                                    if (neighbor.Type == NodeType::NotFluid) {
                                        isBoundaryPoint = true;
                                        break;
                                    }
                                }
                            }
                            if (isBoundaryPoint) {
                                break;
                            }
                        }
                        if (isBoundaryPoint) {
                            break;
                        }
                    }
                    bool isAtBoundary = (i == 0 || i == pointsBox.size() - 1) ||
                        (j == 0 || j == pointsBox[i].size() - 1) ||
                        (k == 0 || k == pointsBox[i][j].size() - 1);

                    if (isAtBoundary && !isBoundaryPoint) {
                        isBoundaryPoint = true;
                    }

                    if (isBoundaryPoint) {
                        boundaryPoints.push_back(point);
                    }
                }
            }
        }
    }

    return boundaryPoints;
}


std::vector<Point> Geometry::excludeInletOutlet() {
    std::vector<Point> remainingPoints;

    for (const Point& p : boundaryPoints) {
        if (std::find(InletboundaryPoints.begin(), InletboundaryPoints.end(), p) == InletboundaryPoints.end() &&
            std::find(OutletboundaryPoints.begin(), OutletboundaryPoints.end(), p) == OutletboundaryPoints.end()) {
            remainingPoints.push_back(p);
        }
    }

    return remainingPoints;
}

void Geometry::getInOutBP() {

    long double num1 = 1e-17L;
    for (int i = 0; i < boundaryPoints.size(); i++) {
        if (boundaryPoints[i].x - InOutCoord[0]<num1) {
            InletboundaryPoints.push_back(boundaryPoints[i]);
        }
        else if (boundaryPoints[i].x == InOutCoord[1]) {
            OutletboundaryPoints.push_back(boundaryPoints[i]);

        }
    }

    WallboundaryPoints = excludeInletOutlet();
}

