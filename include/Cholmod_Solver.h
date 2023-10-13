#pragma once
#include "cholmod.h"
#include "MeshLaplace.h"
#include <iostream>
#include <stdlib.h>
struct Tri
{
	int i;
	int j;
	double x;

	Tri(int _i, int _j, double _x) :i(_i), j(_j), x(_x)
	{}
};
class Cholmod_Solver
{
public:
	Cholmod_Solver(Mesh&_mesh, std::vector<int>& ConeStatus, std::vector<int>& ConeList, 
		std::vector<double>& cur_0, std::vector<double>& area_0, int& conenum, bool isinit);
	~Cholmod_Solver();
	clock_t start, end;

	int vn;
	int *ci, *cj, *cid, *cjd;
	double *cx, *cxd;
	cholmod_common common, *cm;
	cholmod_triplet *ctri, *ctri_d;

	cholmod_sparse *choL, *choB, *choP, *choPI, *choD_inv, *cho_inner_P;
	cholmod_sparse *choB_unsym;
	double *bx;
	int *bp, *bi;
	double sigma;

	cholmod_sparse *Update_C, *Update_CP;
	cholmod_dense *chob0, *chox;
	double *b0x;
	cholmod_factor *chofactorL;
	cholmod_dense *cu, *cp, *ck;
	double *uu, *pp, *kk;

	std::vector<double> m_M, m_KT;
	std::vector<int> v_to_in;
	bool boundary_status = false;

	std::vector<Tri> trips;
	void Initialization(Mesh &mesh);
	void Initialization_with_Boundary(Mesh &mesh);
	void Initialize_Cone(std::vector<int> & ConeStatus, std::vector<double>& cur_0, std::vector<double>& area_0);
	// Update v_i, A'=A + C*C^t
	void Update_Factorization(int id, std::vector<int>& ConeStatus);
	void Solve(std::vector<double>& _u, std::vector<double>& _p, std::vector<double>& _k, std::vector<int>& ConeStatus, double & _distortion);
	void Solve_with_Boundary(std::vector<double>& _u, std::vector<double>& _p, std::vector<double>& _k, std::vector<int>& ConeStatus, double & _distortion);
	void Solve_Zero_Boundary_Condition(std::vector<double>& _u, std::vector<double>& _p, std::vector<double>& _k, std::vector<int>& ConeStatus, double & _distortion);
	void First_Factorization(std::vector<int>& ConeList);
	double one[2], zero[2], negone[2];
	int var_count;
	std::vector<int> v_cur_idx;

	std::vector<double> opt_k;
	std::vector<int> opt_cone;
	double opt_distortion = INFINITY;
private:

};

