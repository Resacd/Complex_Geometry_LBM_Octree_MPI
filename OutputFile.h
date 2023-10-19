#pragma once

#include "Node.h"
#include <fstream>
#include <string>
#include <mpi.h>
#include "Utility.h"

class OutputFile
{
private:

	std::mutex m_fileMutex;
	int m_outputIdx;
	std::mutex m_dataMutex;
	Utility Utility_;

public:

	OutputFile(int outputIdx);
	~OutputFile();


	void WriteVTK(std::vector<std::vector<std::vector<std::shared_ptr<Node>>>>& lattice,
		std::shared_ptr<VelocitySet> velSet, int size_x, int size_y, int size_z,
		int iter);

	void WriteCSV(std::vector<std::vector<std::vector<std::shared_ptr<Node>>>>& lattice,
		std::shared_ptr<VelocitySet> velSet, int start_x, int start_y, int start_z,
		int end_x, int end_y, int end_z);
};

