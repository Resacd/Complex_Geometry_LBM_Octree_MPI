#include "Boundary.h"
#include <cmath>
#include <iostream>

Boundary::Boundary() {



}

Boundary::~Boundary() {


}

double Boundary::Density() {

	return 1;
}


void Boundary::Apply_Equilibrium_Bc(std::shared_ptr<VelocitySet> velSet, Node& lat, std::vector<double>& vel, int dir)
{
        double rho = 1.0;  
        double du = vel[0] * velSet->GetDirection(dir)[0] +
            vel[1] * velSet->GetDirection(dir)[1] +
            vel[2] * velSet->GetDirection(dir)[2];
        double u_sqr = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];
        double equilibriumDist = velSet->GetWeight(dir) * rho *
            (1 + 3 * du + 9.0 / 2.0 * du * du - 3.0 / 2.0 * u_sqr);

        lat.m_distributions[dir] = equilibriumDist;
}


void Boundary::Apply_ZouHe_Bc(std::shared_ptr<VelocitySet> velSet, Node& lat, std::vector<double>& vel, int dir)
{
    if (lat.m_neighbours[dir] == nullptr) {
        double density = Density();
        double du = 0.0;
        double u_sqr = 0.0;
        int opp_dir = velSet->OppositeDirection(dir);

        du = vel[0] * velSet->GetDirection(dir)[0] +
            vel[1] * velSet->GetDirection(dir)[1] +
            vel[2] * velSet->GetDirection(dir)[2];
        u_sqr = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];

        double equilibriumDist = velSet->GetWeight(dir) * density *
            (1 + 3 * du + 9.0 / 2.0 * du * du - 3.0 / 2.0 * u_sqr);

        double equilibriumoppDist = velSet->GetWeight(opp_dir) * density *
            (1 + 3 * du + 9.0 / 2.0 * du * du - 3.0 / 2.0 * u_sqr);


        lat.m_distributions[opp_dir] = (lat.m_distributions[dir] - equilibriumDist) + equilibriumoppDist;

    }
}



void Boundary::Apply_BounceBack(Node& lat, int dir,int opp_dir)
{

    lat.m_distributions[opp_dir] = lat.m_newDistributions[dir];
}

void Boundary::ApplyFreeSlipBC(std::shared_ptr<VelocitySet> velSet, Node& lat, int dir, int opp_dir)
{
    if (lat.m_neighbours[dir] != nullptr) {
        std::swap(lat.m_distributions[dir], lat.m_neighbours[dir]->m_distributions[opp_dir]);
    }
    else {

    }
}


//void Boundary::ApplyPressureDrivenBC(std::shared_ptr<VelocitySet> velSet, Node& lat, int dir, double pressure)
//{
//    int opp_dir = velSet->OppositeDirection(dir);
//    double equilibriumDist = pressure / (velSet->GetWeight(dir) * lat.Density());
//
//    double distributionDiff = equilibriumDist - lat.m_distributions[dir];
//    double bounceBackDist = equilibriumDist - distributionDiff;
//
//    lat.m_distributions[dir] += distributionDiff;
//
//    if (lat.m_neighbours[dir] != nullptr) {
//        lat.m_neighbours[dir]->m_distributions[opp_dir] = bounceBackDist;
//    }
//}


void Boundary::ApplyConstPressure(std::shared_ptr<VelocitySet> velSet, Node& lat, int dir)
{

    if (lat.m_neighbours[dir] == nullptr) {

        int opp_dir = velSet->OppositeDirection(dir);

        if (lat.m_neighbours[opp_dir] != nullptr) {

            lat.m_distributions[opp_dir] = lat.m_neighbours[opp_dir]->m_distributions[opp_dir];
        }
    }
}



