#include "Cholmod_Solver.h"
Cholmod_Solver::Cholmod_Solver(Mesh&mesh, std::vector<int>& ConeStatus, std::vector<int>& ConeList, 
	std::vector<double>& cur_0, std::vector<double>& area_0, int& conenum, bool isinit)
{
	// set boundary label
	vn = mesh.n_vertices();
	cur_0.resize(vn);
	area_0.resize(vn);
	opt_cone.resize(vn);
	opt_k.resize(vn);
	// start cholmod
	cm = &common;
	cholmod_start(cm);

	// ------------------   cholmod parameter --------------
	one[0] = 1;
	one[1] = 0;
	negone[0] = -1;
	negone[1] = 0;
	zero[0] = 0;
	zero[1] = 0;
	// ------------------   cholmod parameter --------------



	// ------------------  cholmod variables  --------------
	chob0 = cholmod_ones(vn, 1, CHOLMOD_REAL, cm);
	cu = cholmod_zeros(vn, 1, CHOLMOD_REAL, cm);
	cp = cholmod_zeros(vn, 1, CHOLMOD_REAL, cm);
	ck = cholmod_zeros(vn, 1, CHOLMOD_REAL, cm);
	b0x = (double *)chob0->x;
	uu = (double *)cu->x;
	pp = (double *)cp->x;
	kk = (double *)ck->x;
	// ------------------  cholmod variables  --------------


	// assemble Lap
	MeshTopo topo(mesh);
	vertArea(mesh.points(), topo, m_M);
	vertGauss(mesh.points(), topo, m_KT);
	std::vector<OpenMesh::Vec3d> f2a, f2cot;
	faceAngle(mesh.points(), topo, f2a);
	f2cot.assign(f2a.size(), OpenMesh::Vec3d(0));
	for (int i = 0; i < f2a.size(); ++i)
	{
		f2cot[i][0] = 0.5 / tan(f2a[i][0]);
		f2cot[i][1] = 0.5 / tan(f2a[i][1]);
		f2cot[i][2] = 0.5 / tan(f2a[i][2]);
	}

	trips.clear();
	for (int i = 0; i < f2a.size(); ++i)
	{
		trips.emplace_back(topo.f2v[i][0], topo.f2v[i][0], (f2cot[i][1] + f2cot[i][2])* 1.0);
		trips.emplace_back(topo.f2v[i][1], topo.f2v[i][1], (f2cot[i][0] + f2cot[i][2])* 1.0);
		trips.emplace_back(topo.f2v[i][2], topo.f2v[i][2], (f2cot[i][0] + f2cot[i][1])* 1.0);
		trips.emplace_back(topo.f2v[i][0], topo.f2v[i][1], -f2cot[i][2] * 1.0);
		trips.emplace_back(topo.f2v[i][1], topo.f2v[i][0], -f2cot[i][2] * 1.0);
		trips.emplace_back(topo.f2v[i][1], topo.f2v[i][2], -f2cot[i][0] * 1.0);
		trips.emplace_back(topo.f2v[i][2], topo.f2v[i][1], -f2cot[i][0] * 1.0);
		trips.emplace_back(topo.f2v[i][2], topo.f2v[i][0], -f2cot[i][1] * 1.0);
		trips.emplace_back(topo.f2v[i][0], topo.f2v[i][2], -f2cot[i][1] * 1.0);
	}
	//-----------------------   laplace matrix ----------------------------------  
	ctri = cholmod_allocate_triplet(vn, vn, trips.size(), 0, CHOLMOD_REAL, cm);
	ctri->nnz = trips.size();

	ci = (int *)ctri->i;
	cj = (int *)ctri->j;
	cx = (double *)ctri->x;

	for (size_t i = 0; i < trips.size(); i++)
	{
		ci[i] = trips[i].i;
		cj[i] = trips[i].j;
		cx[i] = trips[i].x;
	}
	choL = cholmod_triplet_to_sparse(ctri, 0, cm);
	//-----------------------   laplace matrix ----------------------------------  

	//-----------------------   area matrix ---------------------------------- 
	ctri = cholmod_allocate_triplet(vn, vn, vn, 0, CHOLMOD_REAL, cm);
	ctri->nnz = vn;
	ci = (int *)ctri->i;
	cj = (int *)ctri->j;
	cx = (double *)ctri->x;
	for (int i = 0; i < vn; i++)
	{
		b0x[i] = m_KT[i];
		ci[i] = i;
		cj[i] = i;
		cx[i] = pow(m_M[i], -1);
		cur_0[i] = m_KT[i];
		area_0[i] = m_M[i];
	}
	choD_inv = cholmod_triplet_to_sparse(ctri, 0, cm);
	//-----------------------   area matrix ---------------------------------- 


	//-----------------------   inner matrix ---------------------------------- 
	ctri = cholmod_allocate_triplet(vn, vn, vn, -1, CHOLMOD_REAL, cm);
	ctri_d = cholmod_allocate_triplet(vn, vn, vn, -1, CHOLMOD_REAL, cm);
	ctri->nnz = vn;
	ctri_d->nnz = vn;
	ci = (int *)ctri->i;
	cj = (int *)ctri->j;
	cx = (double *)ctri->x;
	cid = (int *)ctri_d->i;
	cjd = (int *)ctri_d->j;
	cxd = (double *)ctri_d->x;
	for (int i = 0; i < vn; i++)
	{
		ci[i] = i;
		cj[i] = i;
		cid[i] = i;
		cjd[i] = i;
		if (mesh.is_boundary(mesh.vertex_handle(i)))
		{
			boundary_status = true;
			cx[i] = 0;
			cxd[i] = 1;
			b0x[i] = 0;
		}
		else

		{
			cxd[i] = 0;
			cx[i] = 1;
		}
	}
	//-----------------------   inner matrix ---------------------------------- 

	if (isinit)
	{
		conenum = 0;
		ConeList.clear();
		for (int i = 0; i < vn; i++)
		{
			if (ConeStatus[i] == 0)
			{
				ConeList.push_back(i);
				conenum++;
			}
		}
	}
	else
	{
		//-----------------------   initialize cones ---------------------------------- 
		ConeStatus.resize(vn, 1);
		ConeList.clear();
		for (int i = 0; i < vn; i++)
		{
			v_cur_idx.push_back(i);
		}
		std::sort(v_cur_idx.begin(), v_cur_idx.end(), [&](const auto& v1_idx, const auto& v2_idx) {
			return abs(m_KT[v1_idx]) > abs(m_KT[v2_idx]);
		});
		conenum = 0;
		for (size_t i = 0; i < vn; ++i)
		{
			if (v_cur_idx[i] > 0)
			{
				if (!mesh.is_boundary(mesh.vertex_handle(v_cur_idx[i])))
				{
					auto vh = mesh.vertex_handle(v_cur_idx[i]);
					bool status = true;
					for (auto vv : mesh.vv_range(vh))
					{
						for (auto vvv : mesh.vv_range(vv))
						{
							for (auto vvvv : mesh.vv_range(vvv))
							{
								for (auto vvvvv : mesh.vv_range(vvvv))
								{
									if (ConeStatus[vvvvv.idx()] == 0)
									{
										status = false;
										break;
									}
								}
							}
						}
					}
					if (status)
					{
						ConeStatus[v_cur_idx[i]] = 0;
						ConeList.push_back(v_cur_idx[i]);
						conenum++;
						//b0x[v_cur_idx[i]] = 0;
						//std::cout << "init cone : " << v_cur_idx[i] + 1 << "  curvature : " << m_KT[v_cur_idx[i]] << std::endl;
					}
				}
				if (conenum > 2) break;       // default : 4,   2(4.20 dataset)
			}
		}
		//-----------------------   initialize cones ---------------------------------- 
	}

	






	//for (size_t i = 0; i < vn; ++i)
	//{
	//	if (v_cur_idx[i] > 0)
	//	{
	//		if (!mesh.is_boundary(mesh.vertex_handle(v_cur_idx[i])))
	//		{
	//			ConeStatus[v_cur_idx[i]] = 0;
	//			ConeList.push_back(v_cur_idx[i]);
	//			conenum++;
	//			//b0x[v_cur_idx[i]] = 0;
	//			std::cout << "init cone : " << v_cur_idx[i] + 1 << "  curvature : " << m_KT[v_cur_idx[i]] << std::endl;
	//		}
	//		if (conenum > 8000) break;       // default : 4,   2(4.20 dataset)
	//	}
	//}




	////-----------------------   random cones ---------------------------------- 
	//conenum = 0;
	//ConeList.clear();
	//ConeStatus.resize(vn, 1);
	//srand(time(nullptr));
	//for (size_t i = 0; i < vn; ++i)
	//{
	//	int rvn = rand()&vn;
	//	if (rvn > 0)
	//	{
	//		if (!mesh.is_boundary(mesh.vertex_handle(rvn)))
	//		{
	//			if (ConeStatus[rvn] == 1)
	//			{
	//				ConeStatus[rvn] = 0;
	//				ConeList.push_back(rvn);
	//				conenum++;
	//				std::cout << "init cone : " << rvn + 1 << "  curvature : " << m_KT[rvn] << std::endl;
	//			}
	//		}
	//		if (conenum > 19) break;       // default : 4
	//	}
	//}
	////-----------------------   random cones ---------------------------------- 




	//// tiny status direct factorize
	//for (size_t i = 0; i < vn; ++i)
	//{
	//	if (v_cur_idx[i] > 0)
	//	{
	//		if (!mesh.is_boundary(mesh.vertex_handle(v_cur_idx[i])))
	//		{
	//			if (ConeStatus[v_cur_idx[i]] == 1)
	//			{
	//				cx[v_cur_idx[i]] = 0;
	//				cxd[v_cur_idx[i]] = 1;
	//				b0x[v_cur_idx[i]] = 0;
	//				ConeStatus[v_cur_idx[i]] = 0;
	//				ConeList.push_back(v_cur_idx[i]);
	//				conenum++;
	//				//b0x[v_cur_idx[i]] = 0;
	//				std::cout << "init cone : " << v_cur_idx[i] + 1 << "  curvature : " << m_KT[v_cur_idx[i]] << std::endl;
	//			}
	//		}
	//		if (conenum > 20) break;       // default : 4
	//	}
	//}
	//ConeStatus[0] = 0;
	//ConeList.push_back(0);
	//conenum++;
	//cho_inner_P = cholmod_triplet_to_sparse(ctri, 0, cm);
	//choP = cholmod_triplet_to_sparse(ctri_d, 0, cm);

	//choB = cholmod_ssmult(choL, choD_inv, 0, 1, 1, cm);
	//choB = cholmod_ssmult(choB, choL, 1, 1, 1, cm);

	//choB = cholmod_ssmult(cho_inner_P, choB, 0, 1, 1, cm);
	//choB = cholmod_ssmult(choB, cho_inner_P, -1, 1, 1, cm);
	//choB = cholmod_add(choB, choP, one, one, 1, 1, cm);

	//choB_unsym = cholmod_ssmult(choL, choD_inv, 0, 1, 1, cm);
	//choB_unsym = cholmod_ssmult(choB_unsym, choL, 0, 1, 1, cm);

	//cholmod_print_sparse(choL, "L", cm);
	//cholmod_print_sparse(choB, "B", cm);
	//cholmod_print_sparse(choB_unsym, "B_unsymmetry", cm);

	//bp = (int *)choB_unsym->p;
	//bi = (int *)choB_unsym->i;
	//bx = (double *)choB_unsym->x;

	//chofactorL = cholmod_analyze(choB, cm);
	//cholmod_factorize(choB, chofactorL, cm);
	//// tiny status direct factorize


	//-----------------------   First factorize (if not boundary, fix first cone) ---------------------------------- 
	if (!boundary_status)
	{
		cx[0] = 0;
		cxd[0] = 1;
		b0x[0] = 0;
		cho_inner_P = cholmod_triplet_to_sparse(ctri, 0, cm);
		choP = cholmod_triplet_to_sparse(ctri_d, 0, cm);

		choB = cholmod_ssmult(choL, choD_inv, 0, 1, 1, cm);
		choB = cholmod_ssmult(choB, choL, 1, 1, 1, cm);

		choB = cholmod_ssmult(cho_inner_P, choB, 0, 1, 1, cm);
		choB = cholmod_ssmult(choB, cho_inner_P, -1, 1, 1, cm);
		choB = cholmod_add(choB, choP, one, one, 1, 1, cm);

		//choB = cholmod_ssmult(cho_inner_P, choL, 0, 1, 1, cm);
		//choB = cholmod_ssmult(choB, choD_inv, 0, 1, 1, cm);
		//choB = cholmod_ssmult(choB, choL, 0, 1, 1, cm);
		//choB = cholmod_ssmult(choB, cho_inner_P, 1, 1, 1, cm);
		//choB = cholmod_add(choB, choP, one, one, 1, 1, cm);
		choB_unsym = cholmod_ssmult(choL, choD_inv, 0, 1, 1, cm);
		choB_unsym = cholmod_ssmult(choB_unsym, choL, 0, 1, 1, cm);
	}
	else
	{
		cho_inner_P = cholmod_triplet_to_sparse(ctri, 0, cm);
		choP = cholmod_triplet_to_sparse(ctri_d, 0, cm);

		choD_inv = cholmod_ssmult(cho_inner_P, choD_inv, 0, 1, 1, cm);
		choD_inv = cholmod_ssmult(choD_inv, cho_inner_P, 0, 1, 1, cm);
		choL = cholmod_ssmult(choL, cho_inner_P, 0, 1, 1, cm);
		choL = cholmod_ssmult(cho_inner_P, choL, 0, 1, 1, cm);

		choB = cholmod_ssmult(choL, choD_inv, 0, 1, 1, cm);
		choB = cholmod_ssmult(choB, choL, 1, 1, 1, cm);
		choB = cholmod_add(choB, choP, one, one, 1, 1, cm);
		choB_unsym = cholmod_ssmult(choL, choD_inv, 0, 1, 1, cm);
		choB_unsym = cholmod_ssmult(choB_unsym, choL, 0, 1, 1, cm);
	}
	//cholmod_print_sparse(choL, "L", cm);
	//cholmod_print_sparse(choB, "B", cm);
	//cholmod_print_sparse(choB_unsym, "B_unsymmetry", cm);

	bp = (int *)choB_unsym->p;
	bi = (int *)choB_unsym->i;
	bx = (double *)choB_unsym->x;

	chofactorL = cholmod_analyze(choB, cm);
	cholmod_factorize(choB, chofactorL, cm);

	std::vector<int> cones(vn, 1);
	if (!boundary_status)
	{
		cones[0] = 0;
		for (size_t i = 0; i < ConeList.size(); i++)
		{
			if (ConeList[i] != 0)
			{
				Update_Factorization(ConeList[i], cones);
				cones[ConeList[i]] = 0;
			}
		}
		Update_Factorization(0, cones);
	}
	else
	{
		for (size_t i = 0; i < ConeList.size(); i++)
		{
			if (ConeList[i] != 0)
			{
				Update_Factorization(ConeList[i], cones);
				cones[ConeList[i]] = 0;
			}
		}
	}
	//-----------------------   First factorize (if not boundary, fix first cone) ---------------------------------- 

	//-----------------------   Initialization over ---------------------------------- 
}

Cholmod_Solver::~Cholmod_Solver()
{
}

void Cholmod_Solver::Initialization(Mesh &mesh)
{
	vn = mesh.n_vertices();
	// start cholmod
	cm = &common;
	cholmod_start(cm);

	// cholmod para
	one[0] = 1;
	one[1] = 0;
	negone[0] = -1;
	negone[1] = 0;
	zero[0] = 0;
	zero[1] = 0;

	// 
	chob0 = cholmod_ones(vn, 1, CHOLMOD_REAL, cm);
	cu = cholmod_zeros(vn, 1, CHOLMOD_REAL, cm);
	cp = cholmod_zeros(vn, 1, CHOLMOD_REAL, cm);
	ck = cholmod_zeros(vn, 1, CHOLMOD_REAL, cm);
	b0x = (double *)chob0->x;
	uu = (double *)cu->x;
	pp = (double *)cp->x;
	kk = (double *)ck->x;
	// assemble Lap
	MeshTopo topo(mesh);
	vertArea(mesh.points(), topo, m_M);
	vertGauss(mesh.points(), topo, m_KT);
	std::vector<OpenMesh::Vec3d> f2a, f2cot;
	faceAngle(mesh.points(), topo, f2a);
	f2cot.assign(f2a.size(), OpenMesh::Vec3d(0));
	for (int i = 0; i < f2a.size(); ++i)
	{
		f2cot[i][0] = 0.5 / tan(f2a[i][0]);
		f2cot[i][1] = 0.5 / tan(f2a[i][1]);
		f2cot[i][2] = 0.5 / tan(f2a[i][2]);
	}

	trips.clear();
	for (int i = 0; i < f2a.size(); ++i)
	{
		trips.emplace_back(topo.f2v[i][0], topo.f2v[i][0], (f2cot[i][1] + f2cot[i][2])* 1.0);
		trips.emplace_back(topo.f2v[i][1], topo.f2v[i][1], (f2cot[i][0] + f2cot[i][2])* 1.0);
		trips.emplace_back(topo.f2v[i][2], topo.f2v[i][2], (f2cot[i][0] + f2cot[i][1])* 1.0);
		trips.emplace_back(topo.f2v[i][0], topo.f2v[i][1], -f2cot[i][2] * 1.0);
		trips.emplace_back(topo.f2v[i][1], topo.f2v[i][0], -f2cot[i][2] * 1.0);
		trips.emplace_back(topo.f2v[i][1], topo.f2v[i][2], -f2cot[i][0] * 1.0);
		trips.emplace_back(topo.f2v[i][2], topo.f2v[i][1], -f2cot[i][0] * 1.0);
		trips.emplace_back(topo.f2v[i][2], topo.f2v[i][0], -f2cot[i][1] * 1.0);
		trips.emplace_back(topo.f2v[i][0], topo.f2v[i][2], -f2cot[i][1] * 1.0);
	}
	//-----------------------   laplace matrix ----------------------------------  
	ctri = cholmod_allocate_triplet(vn, vn, trips.size(), 0, CHOLMOD_REAL, cm);
	ctri->nnz = trips.size();

	ci = (int *)ctri->i;
	cj = (int *)ctri->j;
	cx = (double *)ctri->x;

	for (size_t i = 0; i < trips.size(); i++)
	{
		ci[i] = trips[i].i;
		cj[i] = trips[i].j;
		cx[i] = trips[i].x;
	}
	choL = cholmod_triplet_to_sparse(ctri, 0, cm);
	//-----------------------   laplace matrix ----------------------------------  


	//-----------------------   area matrix ---------------------------------- 
	ctri = cholmod_allocate_triplet(vn, vn, vn, 0, CHOLMOD_REAL, cm);
	ctri->nnz = vn;
	ci = (int *)ctri->i;
	cj = (int *)ctri->j;
	cx = (double *)ctri->x;
	for (int i = 0; i < vn; i++)
	{
		b0x[i] = m_KT[i];
		ci[i] = i;
		cj[i] = i;
		cx[i] = pow(m_M[i], -1);
	}
	choD_inv = cholmod_triplet_to_sparse(ctri, 0, cm);
	//-----------------------   area matrix ---------------------------------- 


	choB = cholmod_ssmult(choL, choD_inv, 0, 1, 1, cm);
	choB = cholmod_ssmult(choB, choL, 1, 1, 1, cm);
	choB_unsym = cholmod_ssmult(choL, choD_inv, 0, 1, 1, cm);
	choB_unsym = cholmod_ssmult(choB_unsym, choL, 0, 1, 1, cm);

	cholmod_print_sparse(choL, "L", cm);
	cholmod_print_sparse(choB, "B", cm);
	cholmod_print_sparse(choD_inv, "D", cm);
	cholmod_print_sparse(choB_unsym, "B_unsymmetry", cm);

	bp = (int *)choB_unsym->p;
	bi = (int *)choB_unsym->i;
	bx = (double *)choB_unsym->x;
}

void Cholmod_Solver::Initialization_with_Boundary(Mesh & mesh)
{
	std::vector<int> Cone_Status;
	// set boundary label
	vn = mesh.n_vertices();

	// start cholmod
	cm = &common;
	cholmod_start(cm);

	// cholmod para
	one[0] = 1;
	one[1] = 0;
	negone[0] = -1;
	negone[1] = 0;
	zero[0] = 0;
	zero[1] = 0;

	// 
	chob0 = cholmod_ones(vn, 1, CHOLMOD_REAL, cm);
	cu = cholmod_zeros(vn, 1, CHOLMOD_REAL, cm);
	cp = cholmod_zeros(vn, 1, CHOLMOD_REAL, cm);
	ck = cholmod_zeros(vn, 1, CHOLMOD_REAL, cm);
	b0x = (double *)chob0->x;
	uu = (double *)cu->x;
	pp = (double *)cp->x;
	kk = (double *)ck->x;
	// assemble Lap
	MeshTopo topo(mesh);
	vertArea(mesh.points(), topo, m_M);
	vertGauss(mesh.points(), topo, m_KT);
	std::vector<OpenMesh::Vec3d> f2a, f2cot;
	faceAngle(mesh.points(), topo, f2a);
	f2cot.assign(f2a.size(), OpenMesh::Vec3d(0));
	for (int i = 0; i < f2a.size(); ++i)
	{
		f2cot[i][0] = 0.5 / tan(f2a[i][0]);
		f2cot[i][1] = 0.5 / tan(f2a[i][1]);
		f2cot[i][2] = 0.5 / tan(f2a[i][2]);
	}

	trips.clear();
	for (int i = 0; i < f2a.size(); ++i)
	{
		trips.emplace_back(topo.f2v[i][0], topo.f2v[i][0], (f2cot[i][1] + f2cot[i][2])* 1.0);
		trips.emplace_back(topo.f2v[i][1], topo.f2v[i][1], (f2cot[i][0] + f2cot[i][2])* 1.0);
		trips.emplace_back(topo.f2v[i][2], topo.f2v[i][2], (f2cot[i][0] + f2cot[i][1])* 1.0);
		trips.emplace_back(topo.f2v[i][0], topo.f2v[i][1], -f2cot[i][2] * 1.0);
		trips.emplace_back(topo.f2v[i][1], topo.f2v[i][0], -f2cot[i][2] * 1.0);
		trips.emplace_back(topo.f2v[i][1], topo.f2v[i][2], -f2cot[i][0] * 1.0);
		trips.emplace_back(topo.f2v[i][2], topo.f2v[i][1], -f2cot[i][0] * 1.0);
		trips.emplace_back(topo.f2v[i][2], topo.f2v[i][0], -f2cot[i][1] * 1.0);
		trips.emplace_back(topo.f2v[i][0], topo.f2v[i][2], -f2cot[i][1] * 1.0);
	}
	//-----------------------   laplace matrix ----------------------------------  
	ctri = cholmod_allocate_triplet(vn, vn, trips.size(), 0, CHOLMOD_REAL, cm);
	ctri->nnz = trips.size();

	ci = (int *)ctri->i;
	cj = (int *)ctri->j;
	cx = (double *)ctri->x;

	for (size_t i = 0; i < trips.size(); i++)
	{
		ci[i] = trips[i].i;
		cj[i] = trips[i].j;
		cx[i] = trips[i].x;
	}
	choL = cholmod_triplet_to_sparse(ctri, 0, cm);
	//-----------------------   laplace matrix ----------------------------------  

	//-----------------------   area matrix ---------------------------------- 
	ctri = cholmod_allocate_triplet(vn, vn, vn, 0, CHOLMOD_REAL, cm);
	ctri->nnz = vn;
	ci = (int *)ctri->i;
	cj = (int *)ctri->j;
	cx = (double *)ctri->x;
	for (int i = 0; i < vn; i++)
	{
		b0x[i] = m_KT[i];
		ci[i] = i;
		cj[i] = i;
		cx[i] = pow(m_M[i], -1);
	}
	choD_inv = cholmod_triplet_to_sparse(ctri, 0, cm);
	//-----------------------   area matrix ---------------------------------- 


	//-----------------------   inner matrix ---------------------------------- 
	ctri = cholmod_allocate_triplet(vn, vn, vn, 0, CHOLMOD_REAL, cm);
	ctri_d = cholmod_allocate_triplet(vn, vn, vn, -1, CHOLMOD_REAL, cm);
	ctri->nnz = vn;
	ctri_d->nnz = vn;
	ci = (int *)ctri->i;
	cj = (int *)ctri->j;
	cx = (double *)ctri->x;
	cid = (int *)ctri_d->i;
	cjd = (int *)ctri_d->j;
	cxd = (double *)ctri_d->x;
	for (int i = 0; i < vn; i++)
	{
		ci[i] = i;
		cj[i] = i;
		cid[i] = i;
		cjd[i] = i;
		if (mesh.is_boundary(mesh.vertex_handle(i)))
		{
			boundary_status = true;
			cx[i] = 0;
			cxd[i] = 1;
			b0x[i] = 0;
		}
		else

		{
			cxd[i] = 0;
			cx[i] = 1;
		}
	}

	//-----------------------   inner matrix ---------------------------------- 

	if (!boundary_status)
	{
		cx[Cone_Status[0]] = 0;
		cxd[Cone_Status[0]] = 1;
		b0x[Cone_Status[0]] = 0;
		cho_inner_P = cholmod_triplet_to_sparse(ctri, 0, cm);
		choP = cholmod_triplet_to_sparse(ctri_d, 0, cm);

		choB = cholmod_ssmult(cho_inner_P, choL, 0, 1, 1, cm);
		choB = cholmod_ssmult(choB, choD_inv, 0, 1, 1, cm);
		choB = cholmod_ssmult(choB, choL, 0, 1, 1, cm);
		choB = cholmod_ssmult(choB, cho_inner_P, 1, 1, 1, cm);
		choB = cholmod_add(choB, choP, one, one, 1, 1, cm);
		choB_unsym = cholmod_ssmult(choL, choD_inv, 0, 1, 1, cm);
		choB_unsym = cholmod_ssmult(choB_unsym, choL, 0, 1, 1, cm);
	}
	else
	{
		cho_inner_P = cholmod_triplet_to_sparse(ctri, 0, cm);
		choP = cholmod_triplet_to_sparse(ctri_d, 0, cm);

		choD_inv = cholmod_ssmult(cho_inner_P, choD_inv, 0, 1, 1, cm);
		choD_inv = cholmod_ssmult(choD_inv, cho_inner_P, 0, 1, 1, cm);
		choL = cholmod_ssmult(choL, cho_inner_P, 0, 1, 1, cm);
		choL = cholmod_ssmult(cho_inner_P, choL, 0, 1, 1, cm);

		choB = cholmod_ssmult(choL, choD_inv, 0, 1, 1, cm);
		choB = cholmod_ssmult(choB, choL, 1, 1, 1, cm);
		choB = cholmod_add(choB, choP, one, one, 1, 1, cm);
		choB_unsym = cholmod_ssmult(choL, choD_inv, 0, 1, 1, cm);
		choB_unsym = cholmod_ssmult(choB_unsym, choL, 0, 1, 1, cm);
	}
	cholmod_print_sparse(choL, "L", cm);
	cholmod_print_sparse(choB, "B", cm);
	cholmod_print_sparse(choB_unsym, "B_unsymmetry", cm);

	bp = (int *)choB_unsym->p;
	bi = (int *)choB_unsym->i;
	bx = (double *)choB_unsym->x;
}



void Cholmod_Solver::Initialize_Cone(std::vector<int>& ConeStatus, std::vector<double>& cur_0, std::vector<double>& area_0)
{
	var_count = 0;
	cur_0.resize(vn);
	area_0.resize(vn);
	// first factorize
	ctri_d = cholmod_allocate_triplet(vn, vn, vn, -1, CHOLMOD_REAL, cm);
	ctri_d->nnz = vn;
	cid = (int *)ctri_d->i;
	cjd = (int *)ctri_d->j;
	cxd = (double *)ctri_d->x;
	ctri = cholmod_allocate_triplet(vn, vn, vn, -1, CHOLMOD_REAL, cm);
	ctri->nnz = vn;
	ci = (int *)ctri->i;
	cj = (int *)ctri->j;
	cx = (double *)ctri->x;
	for (size_t i = 0; i < vn; i++)
	{
		cid[i] = i;
		cjd[i] = i;
		ci[i] = i;
		cj[i] = i;
		if (ConeStatus[i] > 0)
		{
			cxd[i] = 1;
			cx[i] = 0;
		}
		else
		{
			cx[i] = 1;
			cxd[i] = 0;
			b0x[i] = 0;
			var_count++;
		}
		cur_0[i] = m_KT[i];
		area_0[i] = m_M[i];
	}
	choP = cholmod_triplet_to_sparse(ctri_d, 0, cm);
	choPI = cholmod_triplet_to_sparse(ctri, 0, cm);
	choB = cholmod_ssmult(choP, choB, 0, 1, 1, cm);
	choB = cholmod_ssmult(choB, choP, -1, 1, 1, cm);
	choB = cholmod_add(choB, choPI, one, one, 1, 1, cm);
	cholmod_print_sparse(choB, "B", cm);

	chofactorL = cholmod_analyze(choB, cm);
	cholmod_factorize(choB, chofactorL, cm);
}

void Cholmod_Solver::Update_Factorization(int id, std::vector<int>& ConeStatus)
{
	b0x[id] = (1 - ConeStatus[id]) * m_KT[id];
	if (ConeStatus[id] == 0)
	{
		var_count--;
	}
	else
	{
		var_count++;
	}
	// calc eigenvalue 
	int inner_size = 1;
	double delta = 0;
	double diag_value = 1.0 - bx[id], q_value;
	std::vector<int> change_status(bp[id + 1] - bp[id], 0);
	if (ConeStatus[id] == 0)
	{
		diag_value *= -1;
	}
	delta += pow(diag_value, 2);
	int k = 0;
	for (int row = bp[id]; row < bp[id + 1]; row++)
	{
		if ((bi[row] != id))
		{
			change_status[k] = -ConeStatus[bi[row]] * (ConeStatus[id] - abs(1 - ConeStatus[id]));
			delta += 4.0*pow(change_status[k] * bx[row], 2);
			if (change_status[k] != 0)
			{
				inner_size++;
			}
		}
		k++;
	}
	// first sigma_1 > 0
	int k1, k2;
	sigma = 0.5*(diag_value + sqrt(delta));
	q_value = sqrt(sigma / (sigma*sigma + 0.25*(delta - diag_value * diag_value)));
	// assemble C
	ctri = cholmod_allocate_triplet(vn, 1, inner_size, 0, CHOLMOD_REAL, cm);
	ctri->nnz = inner_size;
	ci = (int *)ctri->i;
	cj = (int *)ctri->j;
	cx = (double *)ctri->x;
	k1 = 0;
	k2 = 0;
	for (int row = bp[id]; row < bp[id + 1]; row++)
	{
		if (bi[row] == id)
		{
			ci[k1] = bi[row];
			cj[k1] = 0;
			cx[k1] = sigma * q_value;
			k1++;
		}
		else if (change_status[k2] != 0)
		{
			ci[k1] = bi[row];
			cj[k1] = 0;
			cx[k1] = change_status[k2] * bx[row] * q_value;
			k1++;
		}
		k2++;
	}
	Update_C = cholmod_triplet_to_sparse(ctri, 0, cm);
	Update_CP = cholmod_submatrix(Update_C, (int*)chofactorL->Perm, chofactorL->n, NULL, -1, 1, 1, cm);
	// update
	cholmod_updown(1, Update_CP, chofactorL, cm);

	// second sigma_2 < 0
	sigma = 0.5*(-diag_value + sqrt(delta));
	q_value = sqrt(sigma / (sigma*sigma + 0.25*(delta - diag_value * diag_value)));
	// assemble C
	ctri = cholmod_allocate_triplet(vn, 1, inner_size, 0, CHOLMOD_REAL, cm);
	ctri->nnz = inner_size;
	ci = (int *)ctri->i;
	cj = (int *)ctri->j;
	cx = (double *)ctri->x;
	k1 = 0;
	k2 = 0;
	for (int row = bp[id]; row < bp[id + 1]; row++)
	{
		if (bi[row] == id)
		{
			ci[k1] = bi[row];
			cj[k1] = 0;
			cx[k1] = sigma * q_value;
			k1++;
		}
		else if (change_status[k2] != 0)
		{
			ci[k1] = bi[row];
			cj[k1] = 0;
			cx[k1] = -1.0*change_status[k2] * bx[row] * q_value;
			k1++;
		}
		k2++;
	}
	Update_C = cholmod_triplet_to_sparse(ctri, 0, cm);
	Update_CP = cholmod_submatrix(Update_C, (int*)chofactorL->Perm, chofactorL->n, NULL, -1, 1, 1, cm);
	// downdate
	cholmod_updown(0, Update_CP, chofactorL, cm);
}



void Cholmod_Solver::Solve(std::vector<double>& _u, std::vector<double>& _p, std::vector<double>& _k, std::vector<int>& ConeStatus, double & _distortion)
{
	// check b

	//std::cout << "Cholmod Solver cone number : " << var_count << std::endl;
	double distortion = 0;
	chox = cholmod_solve(CHOLMOD_A, chofactorL, chob0, cm);
	cholmod_sdmult(choL, 1, one, zero, chox, cu, cm);
	cholmod_dense* cuu;
	double *uuu;
	cuu = cholmod_copy_dense(cu, cm);
	uu = (double *)cu->x;
	uuu = (double *)cuu->x;
	kk = (double *)ck->x;
	for (size_t i = 0; i < vn; i++)
	{
		distortion += uu[i] * uu[i] / m_M[i];
		_u[i] = uu[i] / m_M[i];
		//uu[i] = -2.0 * uu[i] / m_M[i];
		uuu[i] = _u[i];
		kk[i] = m_KT[i];
	}
	//cp = cholmod_solve(CHOLMOD_A, chofactorL0, cu, cm);
	////cholmod_sdmult(choL, 1, one, zero, chox, cp, cm);
	cholmod_sdmult(choL, 1, negone, one, cuu, ck, cm);
	pp = (double *)chox->x;
	////
	for (size_t i = 0; i < vn; i++)
	{
		//pp[i] = pp[i] / m_M[i];
		_p[i] = -2.0*pp[i];
		_k[i] = kk[i];
	}
	//std::cout << "Cholmod Solver distortion : " << sqrt(distortion) << std::endl;
	_distortion = sqrt(distortion);

	// update
	if (_distortion < opt_distortion)
	{
		opt_distortion = _distortion;
		for (size_t i = 0; i < vn; i++)
		{
			opt_cone[i] = ConeStatus[i];
			opt_k[i] = _k[i];
		}
	}
}

void Cholmod_Solver::Solve_with_Boundary(std::vector<double>& _u, std::vector<double>& _p, std::vector<double>& _k, std::vector<int>& ConeStatus, double & _distortion)
{
	//std::cout << "Cholmod Solver cone number : " << var_count << std::endl;
	double distortion = 0;
	chox = cholmod_solve(CHOLMOD_A, chofactorL, chob0, cm);
	cholmod_sdmult(choL, 1, one, zero, chox, cu, cm);
	cholmod_dense* cuu;
	double *uuu;
	cuu = cholmod_copy_dense(cu, cm);
	uu = (double *)cu->x;
	uuu = (double *)cuu->x;
	kk = (double *)ck->x;
	for (size_t i = 0; i < vn; i++)
	{
		distortion += uu[i] * uu[i] / m_M[i];
		_u[i] = uu[i] / m_M[i];
		kk[i] = m_KT[i];
		uuu[i] = _u[i];
	}
	cholmod_sdmult(choL, 1, negone, one, cuu, ck, cm);
	////
	for (size_t i = 0; i < vn; i++)
	{
		_k[i] = kk[i];
	}
	std::cout << "Cholmod Solver distortion : " << sqrt(distortion) << std::endl;
	_distortion = sqrt(distortion);
}

void Cholmod_Solver::Solve_Zero_Boundary_Condition(std::vector<double>& _u, std::vector<double>& _p, std::vector<double>& _k, std::vector<int>& ConeStatus, double & _distortion)
{
	//std::cout << "Cholmod Solver cone number : " << var_count << std::endl;
	double distortion = 0;
	chox = cholmod_solve(CHOLMOD_A, chofactorL, chob0, cm);
	//cholmod_sdmult(choL, 1, one, zero, chox, cu, cm);
	cu = cholmod_copy_dense(chox, cm);
	uu = (double *)cu->x;
	kk = (double *)ck->x;
	for (size_t i = 0; i < vn; i++)
	{
		distortion += uu[i] * uu[i] * m_M[i];
		_u[i] = uu[i];
		if (ConeStatus[i] == 0)
		{
			uu[i] = 0;
		}
		else
		{
			uu[i] = -2.0 * uu[i] * m_M[i];
		}		
		//kk[i] = m_KT[i];
	}
	chox = cholmod_solve(CHOLMOD_A, chofactorL, cu, cm);
	cholmod_sdmult(choL, 1, one, zero, chox, cp, cm);
	pp = (double *)cp->x;
	//
	for (size_t i = 0; i < vn; i++)
	{
		_p[i] = pp[i];
	}
	std::cout << "Cholmod Solver distortion : " << sqrt(distortion) << std::endl;
	_distortion = sqrt(distortion);
}

void Cholmod_Solver::First_Factorization(std::vector<int>& ConeList)
{
	std::vector<int> cones(vn, 1);
	//cones[0] = 0;
	for (size_t i = 0; i < ConeList.size(); i++)
	{
		if (ConeList[i] != 0)
		{
			Update_Factorization(ConeList[i], cones);
			cones[ConeList[i]] = 0;
		}
	}
}
