#pragma once
#include "Node.h"

class Boundary
{

private:
	std::vector<double> BoundaryVel_;

public:
	Boundary();
	~Boundary();
	double Density();
	void Apply_Equilibrium_Bc(std::shared_ptr<VelocitySet> velSet, Node& lat, std::vector<double>& vel, int dir);
	void Apply_ZouHe_Bc(std::shared_ptr<VelocitySet> velSet, Node &lat, std::vector<double>& vel, int dir);
	void Apply_BounceBack(Node& lat, int dir, int opp_dir);
	void ApplyFreeSlipBC(std::shared_ptr<VelocitySet> velSet, Node& lat, int dir, int opp_dir);
	void ApplyConstPressure(std::shared_ptr<VelocitySet> velSet, Node& lat, int dir);
	void ApplyPressureDrivenBC(std::shared_ptr<VelocitySet> velSet, Node& lat, int dir, double pressure);
};

