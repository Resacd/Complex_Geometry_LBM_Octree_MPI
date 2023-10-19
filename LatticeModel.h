#pragma once

#include "VelocitySet.h"

class LatticeModel : public VelocitySet
{
public:
	LatticeModel();
	~LatticeModel();

	int OppositeDirection(int dir) const override;
};

