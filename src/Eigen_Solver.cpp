#include "Eigen_Solver.h"
Eigen_Solver::Eigen_Solver(Mesh&mesh, std::vector<int>& ConeStatus, std::vector<int>& ConeList,
	std::vector<double>& cur_0, std::vector<double>& area_0, int& conenum)
{
	clock_t start, end;
	start = clock();
	vn = mesh.n_vertices();
	MeshTopo topo(mesh);
	std::vector<double> vA, vK;
	vertArea(mesh.points(), topo, vA);
	vertGauss(mesh.points(), topo, vK);
	std::vector<Eigen::Triplet<double>> trip;
	laplaceTriplet(mesh.points(), topo, trip);
	m_M = Eigen::Map<VectorX>(vA.data(), topo.vN);
	m_KT = Eigen::Map<VectorX>(vK.data(), topo.vN);
	m_M = m_M / m_M.sum();
	m_lap.resize(topo.vN, topo.vN);
	m_lap.setFromTriplets(trip.begin(), trip.end());
	tildeB0 = m_lap * m_M.cwiseInverse().asDiagonal() * m_lap;
	//solverP1.compute(m_lap*m_lap);
	std::vector<int> v_cur_idx;
	conenum = 0;
	cur_0.resize(vn);
	area_0.resize(vn);
	ConeStatus.resize(vn, 1);
	ConeList.clear();
	for (int i = 0; i < vn; i++)
	{
		v_cur_idx.push_back(i);
		cur_0[i] = m_KT[i];
		area_0[i] = m_M[i];
	}
	std::sort(v_cur_idx.begin(), v_cur_idx.end(), [&](const auto& v1_idx, const auto& v2_idx) {
		return abs(m_KT[v1_idx]) > abs(m_KT[v2_idx]);
	});
	for (size_t i = 0; i < vn; ++i)
	{
		if (v_cur_idx[i] > 0)
		{
			if (!mesh.is_boundary(mesh.vertex_handle(v_cur_idx[i])))
			{
				ConeStatus[v_cur_idx[i]] = 0;
				ConeList.push_back(v_cur_idx[i]);
				conenum++;
				//b0x[v_cur_idx[i]] = 0;
				//std::cout << "init cone : " << v_cur_idx[i] + 1 << "  curvature : " << m_KT[v_cur_idx[i]] << std::endl;
			}
			if (conenum > 5000) break;       // default : 4,   2(4.20 dataset)
		}
	}



	v_b.resize(vn, false);
	trip.clear();
	for (auto v : mesh.vertices())
	{
		if (mesh.is_boundary(v))
		{
			v_b[v.idx()] = true;
			mesh_boundary_status = true;
		}
		else
		{
			trip.emplace_back(v.idx(), v.idx(), 1);
		}
	}
	P_inner.resize(vn, vn);
	P_inner.setFromTriplets(trip.begin(), trip.end());
	m_lap_inner = P_inner * m_lap*P_inner;
	tildeB0_inner = m_lap * P_inner*m_M.cwiseInverse().asDiagonal() * m_lap;
}

Eigen_Solver::Eigen_Solver()
{
}

Eigen_Solver::~Eigen_Solver()
{
}

void Eigen_Solver::calc_distortion(std::vector<double>& _k, double& distortion)
{
	VectorX KK;
	KK.resize(vn);
	for (size_t i = 0; i < vn; i++)
	{
		KK[i] = _k[i];
	}
	VectorX uu = solverP1.solve(m_KT - KK);
	std::cout << "res :" << (m_lap*uu - m_KT + KK).norm() << std::endl;
	uu = uu - ((m_M.cwiseProduct(uu)).sum() / m_M.sum()) *VectorX::Ones(uu.size());
	distortion = sqrt(m_M.dot(uu.cwiseAbs2()));
}

void Eigen_Solver::Initialization(Mesh & mesh)
{
	clock_t start, end;
	start = clock();
	vn = mesh.n_vertices();
	MeshTopo topo(mesh);
	std::vector<double> vA, vK;
	vertArea(mesh.points(), topo, vA);
	vertGauss(mesh.points(), topo, vK);
	std::vector<Eigen::Triplet<double>> trip;
	laplaceTriplet(mesh.points(), topo, trip);
	m_M = Eigen::Map<VectorX>(vA.data(), topo.vN);
	m_KT = Eigen::Map<VectorX>(vK.data(), topo.vN);
	m_M = m_M / m_M.sum();
	m_lap.resize(topo.vN, topo.vN);
	m_lap.setFromTriplets(trip.begin(), trip.end());
	tildeB0 = m_lap * m_M.cwiseInverse().asDiagonal() * m_lap;
	solverP1.compute(m_lap*m_lap);

	v_b.resize(vn, false);
	trip.clear();
	for (auto v : mesh.vertices())
	{
		if (mesh.is_boundary(v))
		{
			v_b[v.idx()] = true;
			mesh_boundary_status = true;
		}
		else
		{
			trip.emplace_back(v.idx(), v.idx(), 1);
		}
	}
	P_inner.resize(vn, vn);
	P_inner.setFromTriplets(trip.begin(), trip.end());
	m_lap_inner = P_inner * m_lap*P_inner;
	tildeB0_inner = m_lap * P_inner*m_M.cwiseInverse().asDiagonal() * m_lap;
	//solverP1.compute(tildeB0);
	opt_cone.resize(vn);
	opt_k.resize(vn);

	solverP1.compute(m_lap);
}

void Eigen_Solver::ResetInnerMatrix(std::vector<int>& ConeStatus)
{
	std::vector<double> vb;
	std::vector<Eigen::Triplet<double>> trip;
	int row_id = 0;
	for (int i = 0; i < vn; i++)
	{
		if (ConeStatus[i] == 0 || (v_b[i])) continue;
		else
		{
			trip.emplace_back(row_id, i, 1);
			vb.push_back(m_KT[i]);
			row_id++;
		}
	}
	P.resize(row_id, vn);
	P.setFromTriplets(trip.begin(), trip.end());
	//std::cout << "Eigen Solver, cone number : " << (vn- row_id) << std::endl;
}

void Eigen_Solver::Solve(std::vector<double>& _u, std::vector<double>& _p, std::vector<double>& _k, std::vector<int>& ConeStatus)
{
	ResetInnerMatrix(ConeStatus);
	solverB1.compute(P * tildeB0*P.transpose());
	//solverB1.compute(P * m_lap*P.transpose());
	//u_temp = P.transpose()*solverB1.solve(P * m_KT);
	u_temp = m_M.cwiseInverse().asDiagonal()*m_lap*P.transpose()*solverB1.solve(P * m_KT);
	p_temp = m_M.cwiseInverse().asDiagonal()*m_lap*P.transpose()*solverB1.solve(-2.0*P*m_M.cwiseSqrt().asDiagonal()*u_temp);
	K_temp = -m_lap * u_temp + m_KT;
	for (size_t i = 0; i < vn; i++)
	{
		_u[i] = u_temp[i];
		_p[i] = p_temp[i];
		_k[i] = K_temp[i];
	}
	std::cout << "Eigen Solver, distortion : " << sqrt(u_temp.transpose()*m_M.asDiagonal()*u_temp) << std::endl;
}

void Eigen_Solver::Solve(std::vector<double>& _u, std::vector<double>& _p, std::vector<double>& _k, std::vector<int>& ConeStatus, double & _distortion)
{
	ResetInnerMatrix(ConeStatus);
	solverB1.compute(P * tildeB0*P.transpose());
	//solverB1.compute(P * m_lap*P.transpose());
	//u_temp = P.transpose()*solverB1.solve(P * m_KT);
	u_temp = m_M.cwiseInverse().asDiagonal()*m_lap*P.transpose()*solverB1.solve(P * m_KT);
	p_temp = m_M.cwiseInverse().asDiagonal()*m_lap*P.transpose()*solverB1.solve(-2.0*P*m_M.cwiseSqrt().asDiagonal()*u_temp);
	K_temp = -m_lap * u_temp + m_KT;
	for (size_t i = 0; i < vn; i++)
	{
		_u[i] = u_temp[i];
		_p[i] = p_temp[i];
		_k[i] = K_temp[i];
	}
	_distortion = sqrt(u_temp.transpose()*m_M.asDiagonal()*u_temp);
	//std::cout << "Eigen Solver, distortion : " << sqrt(u_temp.transpose()*m_M.asDiagonal()*u_temp) << std::endl;
}

void Eigen_Solver::Solve_Dirichlet(std::vector<double>& _u, std::vector<double>& _p, std::vector<double>& _k, std::vector<int>& Cone_Status, double& area_distortion)
{
	ResetInnerMatrix(Cone_Status);
	solverB1.compute(P * tildeB0_inner*P.transpose());
	//solverB1.compute(P * m_lap*P.transpose());
	//u_temp = P.transpose()*solverB1.solve(P * m_KT);
	u_temp = m_M.cwiseInverse().asDiagonal()*m_lap_inner*P.transpose()*solverB1.solve(P * m_KT);

	for (size_t i = 0; i < vn; i++)
	{
		if (v_b[i])
		{
			//std::cout << u_temp[i] << std::endl;
		}
	}
	//p_temp = m_M.cwiseInverse().asDiagonal()*m_lap*P.transpose()*solverB1.solve(-2.0*P*m_M.cwiseSqrt().asDiagonal()*u_temp);
	p_temp = solverP1.solve(-2.0*m_lap*m_M.asDiagonal()*u_temp);
	K_temp = -m_lap_inner * u_temp + m_KT;
	for (size_t i = 0; i < vn; i++)
	{
		_u[i] = u_temp[i];
		_p[i] = p_temp[i];
		_k[i] = K_temp[i];
	}
	area_distortion = sqrt(u_temp.transpose()*m_M.asDiagonal()*u_temp);
	//std::cout << "Eigen Solver, distortion : " << area_distortion << std::endl;

	if (area_distortion < opt_distortion)
	{
		opt_distortion = area_distortion;
		for (size_t i = 0; i < vn; i++)
		{
			opt_cone[i] = Cone_Status[i];
			opt_k[i] = _k[i];
		}
	}
}

void Eigen_Solver::SolveP(std::vector<double>& _u, std::vector<double>& _p, std::vector<double>& _k, std::vector<int>& ConeStatus)
{
	u_temp.resize(vn);
	p_temp.resize(vn);
	//ResetInnerMatrix(ConeStatus);
	//solverP1.compute(P * m_lap*P.transpose());
	for (size_t i = 0; i < vn; i++)
	{
		u_temp[i] = _u[i];
		p_temp[i] = _p[i];
	}
	std::cout << " resdual 1 : " << (m_lap *p_temp + 2.0*m_M.asDiagonal()*u_temp).norm() << std::endl;


	p_temp = solverP1.solve(-2.0*m_lap*m_M.asDiagonal()*u_temp);
	VectorX resd = m_lap *p_temp + 2.0*m_M.asDiagonal()*u_temp;
	std::vector<double> resdv(vn);
	for (size_t i = 0; i < vn; i++)
	{
		resdv[i] = resd[i];
	}
	std::cout << " resdual 2 : " << (m_lap *p_temp + 2.0*m_M.asDiagonal()*u_temp).norm() << std::endl;
	//p_temp = P.transpose()*solverP1.solve(-2.0*m_M.asDiagonal()*u_temp);
	for (size_t i = 0; i < vn; i++)
	{
		_p[i] = p_temp[i];
	}
}


