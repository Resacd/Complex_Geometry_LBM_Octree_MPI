#pragma once
#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <iostream>
#include <vector>
#include <unordered_set>
#include "Struct.h"
#include <tuple>
#include <map>
#include <set>

class Utility {

public:
    Utility();
    ~Utility();
    bool isPointInBoundary(const Point& targetPoint, const std::vector<Point>& boundaryPoints);
};



#endif // !_UTILITY_
