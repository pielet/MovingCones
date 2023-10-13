#pragma once
#include "Cholmod_Solver.h"
#include "Eigen_Solver.h"
#include <fstream>
#include <iostream>
#include <deque>

struct OneRing
{
	int id;
	int is_cone;
	double value;
	double u;
	OpenMesh::Vec3d vec;
	OneRing(int id, int is_cone, OpenMesh::Vec3d vec, double u, double value)
		:id(id), is_cone(is_cone), vec(vec), u(u), value(value) {};
};
struct MovePair
{
	int src;
	int tar;
	double value;
	int list_id;
	MovePair(int _src, int _tar, double _value, int _list_id)
		:src(_src), tar(_tar), value(_value), list_id(_list_id){};
};
struct GraphList
{
	int max_id;
	double max_u;
	double range_size;
	int cone_num;
	std::vector<int> v_list;
	std::vector<int> cone_list;
	bool sign;
	
	GraphList(int max_id, double max_u, double range_size, int cone_num, std::vector<int> v_list, std::vector<int> cone_list, bool sign)
		:max_id(max_id), max_u(max_u), range_size(range_size), cone_num(cone_num), v_list(v_list), cone_list(cone_list), sign(sign) {};
};


class ConeOptimization
{
public:
	ConeOptimization();
	~ConeOptimization();

	// ------------   Initial
	void Initialization(Mesh& _mesh, double distortion);

	void Optimization();
	
	// ------------   Cone Adding
	bool is_add_zero = false;
	void Cone_Adding();
	void Cone_Add_tiny_bounds();

	// ------------   Cone Moving
	int Calc_Gradient(int id, double &value);
	std::vector<OneRing> one_ring;


	// ------------   Info
	Mesh mesh;

	int vn, cone_num = 1;
	double area_distortion, distortion_old, input_distortion;

	std::vector<double> uc, pc, kc, Area_V;
	std::vector<int> ConeStatus, ConeList, old_cone;
	void ResetConeId();

	std::vector<double> ori_cur;
	std::vector<int> v_cur_idx;

	// ------------   Solver 
	int solver_choose = 0;
	void ChangeSolver();
	Cholmod_Solver *cholmod_s;
	Eigen_Solver *eigen_s;

	// ------------   Parameter
	double move_para = 0.5;
	int iter_num_move = 1e4; 
	int change_solver_number = 1500;
};

