#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <ctype.h>
#include "ConeOptimization.h"


void printHelp(const std::string& programName)
{
	std::cout << "Help: "
		<< programName
		<< " input_obj_file "
		<< "output_cones_file "
		<< "[--normBound=sigma] "
		<< std::endl;
}

bool doesArgExist(const std::string& arg, const std::string& searchStr)
{
	return arg.find(searchStr) != std::string::npos;
}

bool parseArg(const std::string& arg, const std::string& searchStr, std::string& value)
{
	if (doesArgExist(arg, searchStr)) {
		value = arg.substr(arg.find_first_of(searchStr[searchStr.size() - 1]) + 1);
		return true;
	}
	return false;
}

void parseArgs(int argc, const char* argv[], std::string& objPath, std::string& conesPath, double & distortion)
{
	if (argc < 3) {
		// input and/or output path not specified
		printHelp(argv[0]);
		exit(EXIT_FAILURE);

	}
	else {
		// parse arguments
		objPath = argv[1];
		conesPath = argv[2];
		distortion = std::stod(argv[3]);
	}
}

void saveCones(const std::vector<double>& conesK, std::string conesPath, Mesh& mesh, double eps = 1e-9)
{
	std::ofstream conesFile(conesPath);
	if (conesFile.fail())
	{
		std::cout << "Open " << conesPath << "failed\n";
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < conesK.size(); ++i)
	{
		if ((conesK[i] > -eps && conesK[i] < eps)||mesh.is_boundary(mesh.vertex_handle(i))) continue;
		conesFile << i + 1 << " " << conesK[i] << std::endl;
	}
	conesFile.close();
}
void saveInfo(std::string conesPath, double costTime, double u, int conenum)
{
	std::ofstream conesFile(conesPath);
	if (conesFile.fail())
	{
		std::cout << "Open " << conesPath;
		exit(EXIT_FAILURE);
	}
	conesFile << "Cost time : " << costTime << "\n";
	conesFile << "distortion : " << u << "\n";
	conesFile << "Cone Number : " << conenum << "\n";
	conesFile.close();
}

int main(int argc, const char *argv[])
{
	std::string objPath = "";
	std::string conesPath = "";
	double distortion;
	parseArgs(argc, argv, objPath, conesPath, distortion);

	Mesh mesh;
	std::cout << "load mesh from " << objPath << std::endl;
	if (!MeshTools::ReadMesh(mesh, objPath))
	{
		std::cout << "load failed!\n";
		exit(EXIT_FAILURE);
	}

	clock_t start, end;
	ConeOptimization ConeOpt;
	start = clock();
	ConeOpt.Initialization(mesh, distortion);
	ConeOpt.Optimization();
	end = clock();
	double costTime = (double)(end - start) / CLOCKS_PER_SEC;
	std::cout << "Cost time : " << costTime << "\n";
	saveCones(ConeOpt.kc, conesPath+"-cones.txt", mesh);
	return 0;
}