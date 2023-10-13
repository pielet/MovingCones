#pragma once
#include "MeshLaplace.h"
#include "Types.h"
#include <iostream>
#include <ctime>
class Eigen_Solver
{
public:
	Eigen_Solver(Mesh&mesh, std::vector<int>& ConeStatus, std::vector<int>& ConeList,
		std::vector<double>& cur_0, std::vector<double>& area_0, int& conenum);
	Eigen_Solver();
	~Eigen_Solver();
	Eigen::SimplicialLDLT<ColMajorSparseMatrix> solverB1;
	Eigen::SimplicialLDLT<ColMajorSparseMatrix> solverP1;
	int vn;
	bool mesh_boundary_status = false;
	ColMajorSparseMatrix P, tildeB0, m_lap, P_inner, m_lap_inner;
	ColMajorSparseMatrix tildeB0_inner;
	VectorX m_M, m_KT;
	VectorX K_temp;
	VectorX u_temp;
	VectorX p_temp;

	std::vector<double> opt_k;
	std::vector<int> opt_cone;
	double opt_distortion = INFINITY;


	void calc_distortion(std::vector<double>& _k, double& distortion);
	std::vector<bool> v_b;
	void Initialization(Mesh &_mesh);
	void ResetInnerMatrix(std::vector<int> & ConeStatus);
	// Update v_i, A'=A + C*C^t
	void Solve(std::vector<double>& _u, std::vector<double>& _p, std::vector<double>& _k, std::vector<int>& Cone_Status);
	void Solve(std::vector<double>& _u, std::vector<double>& _p, std::vector<double>& _k, std::vector<int>& ConeStatus, double & _distortion);
	void Solve_Dirichlet(std::vector<double>& _u, std::vector<double>& _p, std::vector<double>& _k, std::vector<int>& Cone_Status, double& area_distortion);
	void SolveP(std::vector<double>& _u, std::vector<double>& _p, std::vector<double>& _k, std::vector<int>& Cone_Status);

private:

};

