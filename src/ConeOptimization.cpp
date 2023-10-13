#include "ConeOptimization.h"


ConeOptimization::ConeOptimization()
{
}

ConeOptimization::~ConeOptimization()
{
}

void ConeOptimization::Initialization(Mesh & _mesh, double distortion)
{
	mesh = _mesh;
	vn = mesh.n_vertices();
	uc.resize(vn);
	pc.resize(vn);
	kc.resize(vn);

	cholmod_s = new Cholmod_Solver(mesh, ConeStatus, ConeList, ori_cur, Area_V, cone_num, false);
	input_distortion = distortion;
	eigen_s = new Eigen_Solver();
}

int ConeOptimization::Calc_Gradient(int id, double &value)
{
	one_ring.clear();
	auto vh = mesh.vertex_handle(id);
	for (auto vv : mesh.vv_range(vh))
	{
		OneRing temp_one_ring(vv.idx(), ConeStatus[vv.idx()], mesh.point(vv) - mesh.point(vh)
			, -(uc[vv.idx()] - uc[vh.idx()])*(pc[vv.idx()] - pc[vh.idx()]) / pow((mesh.point(vv) - mesh.point(vh)).norm(), 2), 0.0);
		one_ring.push_back(temp_one_ring);
	}
	for (size_t i = 0; i < one_ring.size(); i++)
	{
		for (size_t j = 0; j < one_ring.size(); j++)
		{
			one_ring[j].value += (one_ring[i].vec | one_ring[j].vec)*one_ring[i].u / (one_ring[j].vec.norm()*one_ring[i].vec.norm());
		}
	}
	std::sort(one_ring.begin(), one_ring.end(), [](OneRing x, OneRing y) -> bool {
		return x.value > y.value;
	});
	for (size_t i = 0; i < one_ring.size(); i++)
	{
		if (one_ring[i].is_cone > 0 && (!mesh.is_boundary(mesh.vertex_handle(one_ring[i].id))))
		{
			double factor = (one_ring[i].value - one_ring[one_ring.size() - 1].value) / one_ring[i].value;
			if (factor < move_para)
			{
				return -1;
				value = 0;
			}
			else
			{
				value = one_ring[i].value;
				return one_ring[i].id;
			}
		}
	}
}


void ConeOptimization::Optimization()
{
	for (size_t iteration = 0; iteration < iter_num_move; iteration++)
	{
		ResetConeId();
		if (cone_num > change_solver_number && (solver_choose==0))
		{
			ChangeSolver();
		}
		if (solver_choose == 0)
		{
			cholmod_s->Solve(uc, pc, kc, ConeStatus, area_distortion);
		}
		else
		{
			eigen_s->Solve_Dirichlet(uc, pc, kc, ConeStatus, area_distortion);
		}
		int movecount = 0;
		std::vector<MovePair> move_list;
		for (size_t i = 0; i < ConeList.size(); i++)
		{
			int id = ConeList[i];
			int opt_id = -1;
			double value = 0;
			opt_id = Calc_Gradient(id, value);
			if (opt_id > -1)
			{
				MovePair temp(id, opt_id, value, i);
				move_list.push_back(temp);
			}
		}
		std::sort(move_list.begin(), move_list.end(), [](MovePair x, MovePair y) -> bool {
			return x.value > y.value;
		});
		for (size_t i = 0; i < move_list.size(); i++)
		{
			int src = move_list[i].src;
			int tar = move_list[i].tar;
			int lid = move_list[i].list_id;
			if (ConeStatus[tar] > 0)
			{
				if (solver_choose == 0)
				{
					cholmod_s->Update_Factorization(src, ConeStatus);
					ConeStatus[src] = 1;
					cholmod_s->Update_Factorization(tar, ConeStatus);
					ConeStatus[tar] = 0;
				}

				else
				{
					ConeStatus[src] = 1;
					ConeStatus[tar] = 0;
				}

				movecount++;
				ConeList[lid] = tar;
			}
			if (move_list[i].value < 5e-1)
			{
				break;
			}
		}
		if (iteration > 0 && (!(distortion_old > area_distortion)))
		{	
			if (distortion_old < input_distortion) break;		
			if (cone_num > change_solver_number)
			{
				Cone_Add_tiny_bounds();
			}
			else
			{
				Cone_Adding();
			}	
		}
		if ((iteration % 10 == 9) && (area_distortion > input_distortion))
		{
			if (cone_num > change_solver_number)
			{
				Cone_Add_tiny_bounds();
			}
			else
			{
				Cone_Adding();
			}	
		}
		distortion_old = area_distortion;
	}
	
	if (solver_choose == 0)
	{
		kc = cholmod_s->opt_k;
		ConeStatus = cholmod_s->opt_cone;
		area_distortion = cholmod_s->opt_distortion;
	}
	else
	{
		kc = eigen_s->opt_k;
		ConeStatus = eigen_s->opt_cone;
		area_distortion = eigen_s->opt_distortion;
	}
	std::cout << "Distortion : " << area_distortion << "  cone number : " << cone_num << std::endl;
}





void ConeOptimization::ChangeSolver()
{
	eigen_s->Initialization(mesh);
	solver_choose = 1;
	std::cout << " Solver changed !" << std::endl;
}


void ConeOptimization::ResetConeId()
{
	cone_num = 0;
	ConeList.clear();
	for (size_t i = 0; i < vn; i++)
	{
		if (ConeStatus[i] == 0)
		{
			cone_num++;
			ConeList.push_back(i);
		}
	}
}

void ConeOptimization::Cone_Adding()
{
	int add_num = int(cone_num*(area_distortion - input_distortion) / input_distortion);
	if (area_distortion < 0.25)
	{
		add_num = 1;
	}
	else if (area_distortion < 0.3)
	{
		add_num = 2;
	}
	else if (area_distortion < 0.4)
	{
		add_num = 3;
	}
	else if (area_distortion < 0.5)
	{
		add_num = 5;
	}
	else if (area_distortion < 0.7)
	{
		add_num = 7;
	}
	else if (area_distortion < 1)
	{
		add_num = 10;
	}
	else
	{
		add_num = 20;
	}
	if (cone_num < 20)
	{
		add_num = 1;
	}
	if (input_distortion < 1e-3)
	{
		add_num = 10;
	}

	int add_count = 0;
	std::vector<int> visit_status(vn, -1);     // -1: not visited, 0 : visited, 1 : visited unsigned, 2 : visited unsigned revisited 
	std::vector<bool> filter_status(vn, false);
	std::vector<GraphList> max_list;
	visit_status[0] = 1;
	std::deque<int> v_q;
	v_q.push_back(0);
	int list_pos = -1;
	double filter_ratio = area_distortion;

	if (fabs(uc[0]) > filter_ratio)
	{
		filter_status[0] = true;
		double max_temp_u = fabs(uc[0]);
		if (ConeStatus[0] == 0 || mesh.is_boundary(mesh.vertex_handle(0))) max_temp_u = 0;
		GraphList temp(0, max_temp_u, uc[0] * uc[0] * Area_V[0], 0, {}, {}, (uc[0] > 0));
		max_list.push_back(temp);
		list_pos++;
		if (1 - ConeStatus[0])
		{
			max_list[list_pos].cone_num++;
			max_list[list_pos].cone_list.push_back(0);
		}
		max_list[list_pos].v_list.push_back(0);
	}
	else
	{
		filter_status[0] = false;
	}
	while (!v_q.empty())
	{
		int front = v_q.front();
		v_q.pop_front();

		if (visit_status[front] == 2)
		{
			double max_temp_u = fabs(uc[front]);
			if (ConeStatus[front] == 0) max_temp_u = 0;
			GraphList temp(front, max_temp_u, uc[front] * uc[front] * Area_V[front], 0, {}, {}, (uc[front] > 0));
			max_list.push_back(temp);
			list_pos++;
			if (1 - ConeStatus[front])
			{
				max_list[list_pos].cone_num++;
				max_list[list_pos].cone_list.push_back(front);
			}
			max_list[list_pos].v_list.push_back(front);
			visit_status[front] = 1;
		}
		for (auto vv : mesh.vv_range(mesh.vertex_handle(front)))
		{
			if (visit_status[vv.idx()] < 0)
			{
				if (fabs(uc[vv.idx()]) > filter_ratio)
				{
					if (!filter_status[front])
					{
						filter_status[vv.idx()] = true;
						v_q.push_front(vv.idx());
						double max_temp_u = fabs(uc[vv.idx()]);
						if (ConeStatus[vv.idx()] == 0 || mesh.is_boundary(vv)) max_temp_u = 0;
						GraphList temp(vv.idx(), max_temp_u, uc[vv.idx()] * uc[vv.idx()] * Area_V[vv.idx()], 0, {}, {}, (uc[vv.idx()] > 0));
						max_list.push_back(temp);
						list_pos++;
						if (1 - ConeStatus[vv.idx()])
						{
							max_list[list_pos].cone_num++;
							max_list[list_pos].cone_list.push_back(vv.idx());
						}
						max_list[list_pos].v_list.push_back(vv.idx());
						visit_status[vv.idx()] = 1;
						break;
					}
					else if ((uc[vv.idx()] > 0) == max_list[list_pos].sign)
					{
						filter_status[vv.idx()] = true;
						v_q.push_front(vv.idx());
						max_list[list_pos].range_size += uc[vv.idx()] * uc[vv.idx()] * Area_V[vv.idx()];
						//max_list[list_pos].range_size += 1;
						if (fabs(uc[vv.idx()]) > max_list[list_pos].max_u && ConeStatus[vv.idx()] == 1 && (!mesh.is_boundary(vv)))
						{
							max_list[list_pos].max_u = fabs(uc[vv.idx()]);
							max_list[list_pos].max_id = vv.idx();
						}
						if (1 - ConeStatus[vv.idx()])
						{
							max_list[list_pos].cone_num++;
							max_list[list_pos].cone_list.push_back(vv.idx());
						}
						max_list[list_pos].v_list.push_back(vv.idx());
						visit_status[vv.idx()] = 1;
					}
					else
					{
						filter_status[vv.idx()] = true;
						visit_status[vv.idx()] = 2;
						v_q.push_back(vv.idx());
					}
				}
				else
				{
					v_q.push_back(vv.idx());
					visit_status[vv.idx()] = 1;
				}
			}
			else if (visit_status[vv.idx()] == 2)
			{
				if (!filter_status[front])
				{
					filter_status[vv.idx()] = true;
					v_q.push_front(vv.idx());
					double max_temp_u = fabs(uc[vv.idx()]);
					if (ConeStatus[vv.idx()] == 0 || mesh.is_boundary(vv)) max_temp_u = 0;
					GraphList temp(vv.idx(), max_temp_u, uc[vv.idx()] * uc[vv.idx()] * Area_V[vv.idx()], 0, {}, {}, (uc[vv.idx()] > 0));
					max_list.push_back(temp);
					list_pos++;
					if (1 - ConeStatus[vv.idx()])
					{
						max_list[list_pos].cone_num++;
						max_list[list_pos].cone_list.push_back(vv.idx());
					}
					max_list[list_pos].v_list.push_back(vv.idx());
					visit_status[vv.idx()] = 1;
					break;
				}
				if ((uc[vv.idx()] > 0) == max_list[list_pos].sign)
				{
					filter_status[vv.idx()] = true;
					v_q.push_front(vv.idx());
					max_list[list_pos].range_size += uc[vv.idx()] * uc[vv.idx()] * Area_V[vv.idx()];
					//max_list[list_pos].range_size += 1;
					if (fabs(uc[vv.idx()]) > max_list[list_pos].max_u && ConeStatus[vv.idx()] == 1 && (!mesh.is_boundary(vv)))
					{
						max_list[list_pos].max_u = fabs(uc[vv.idx()]);
						max_list[list_pos].max_id = vv.idx();
					}
					if (1 - ConeStatus[vv.idx()])
					{
						max_list[list_pos].cone_num++;
						max_list[list_pos].cone_list.push_back(vv.idx());
					}
					max_list[list_pos].v_list.push_back(vv.idx());
					visit_status[vv.idx()] = 1;
				}
			}
		}
	}

	std::vector<int> sort_id;
	for (int i = 0; i < max_list.size(); i++)
	{
		sort_id.push_back(i);
	}
	std::sort(sort_id.begin(), sort_id.end(), [&](const auto& v1_idx, const auto& v2_idx) {
		return max_list[v1_idx].range_size / (max_list[v1_idx].cone_num + 1) > max_list[v2_idx].range_size / (max_list[v2_idx].cone_num + 1);
		//return max_list[v1_idx].range_size  > max_list[v2_idx].range_size ;
	});

	std::vector<int> update_list;
	for (size_t i = 0; i < max_list.size(); i++)
	{
		if (max_list[sort_id[i]].cone_num < 1)
		{
			update_list.push_back(max_list[sort_id[i]].max_id);
			add_count++;
		}

		if (add_count == add_num) break;
	}
	if (update_list.size() < 1)
	{
		std::cout << " add wrong !" << std::endl;
		for (size_t i = 0; i < add_num; i++)
		{
			update_list.push_back(max_list[sort_id[i]].max_id);
		}
	}


	add_count = 0;
	for (size_t i = 0; i < update_list.size(); i++)
	{
		// check 2-ring
		bool is_update = true;
		for (auto vv : mesh.vv_range(mesh.vertex_handle(update_list[i])))
		{
			if (ConeStatus[vv.idx()] == 0) is_update = false;
			for (auto vvv : mesh.vv_range(vv))
			{
				if (ConeStatus[vv.idx()] == 0) is_update = false;
			}
		}
		if (is_add_zero) is_update = true;
		if (is_update)
		{
			if (solver_choose == 0)
			{
				cholmod_s->Update_Factorization(update_list[i], ConeStatus);
			}

			ConeStatus[update_list[i]] = 0;
			ConeList.push_back(update_list[i]);
			cone_num++;
			add_count++;
		}
	}
	is_add_zero = (add_count == 0);
}




void ConeOptimization::Cone_Add_tiny_bounds()
{
	int addnum = 50;
	int addcount = 0;
	v_cur_idx.clear();
	for (int i = 0; i < vn; i++)
	{
		v_cur_idx.push_back(i);
	}
	std::sort(v_cur_idx.begin(), v_cur_idx.end(), [&](const auto& v1_idx, const auto& v2_idx) {
		return abs(ori_cur[v1_idx]) > abs(ori_cur[v2_idx]);
	});
	for (size_t i = 0; i < vn; i++)
	{
		int id = v_cur_idx[i];
		if (ConeStatus[id] == 1)
		{
			ConeStatus[id] = 0;
			ConeList.push_back(id);
			addcount++;
			cone_num++;
		}
		if (addcount > addnum) break;
	}
}

