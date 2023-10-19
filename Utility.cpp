#pragma once
#include "Utility.h"

Utility::Utility() {


}

Utility::~Utility() {


}

bool Utility::isPointInBoundary(const Point& targetPoint, const std::vector<Point>& boundaryPoints) {
    for (const Point& boundaryPoint : boundaryPoints) {
        if (targetPoint == boundaryPoint) {
            return true; 
        }
    }
    return false; 
}