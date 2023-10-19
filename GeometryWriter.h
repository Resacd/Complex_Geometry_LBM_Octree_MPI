#pragma once
#pragma once
#ifndef _WRITER_H_
#define _WRITER_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <array>
#include "Struct.h"
#include <unordered_map>

class Writer {

public:

    Writer();
    ~Writer();
    void writePointoffToFile(const std::vector<Point>& points, const std::string& filename);
    void writePointsToFile(const std::vector<Point>& pointsBox, const std::string& filename);
    void writeBoxToFile(const std::vector<std::vector<std::vector<Point>>>& points, const std::string& filename);
    void writePointsToFile(const std::vector<std::vector<std::vector<Point>>>& pointsBox, const std::string& filename);
    void writeLabeledPointsToFile(const std::vector<std::vector<std::vector<Point>>>& pointsBox, const std::string& filename);
    void writeLineToVTK(const Line& line, const std::string& filename);
};	

#endif 