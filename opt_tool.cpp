/*** 
 * ZICEN XIONG 
 * 2023-12-17 19:49:04
 */
# include "opt_tool.h"

struct Particle
{	//
	Eigen::Matrix<double, 1, Eigen::Dynamic> position;
	Eigen::Matrix<double, 1, Eigen::Dynamic> velocity;
	Eigen::Matrix<double, 1, Eigen::Dynamic> pBest;
	double fitness;
	double pBestFitness;
};

std::pair<Eigen::MatrixXd, double> CLPSO(const std::function<double(const Eigen::MatrixXd&)>& fhd, int Dimension, const Eigen::MatrixXd& Rmin, const Eigen::MatrixXd& Rmax, int Max_Gen, int Particle_Number)
{
	int ps = Particle_Number; // ��ʼ��Ⱥ����
	int D = Dimension;        // �ռ�ά��
	int me = Max_Gen;         // ����������

	if (Rmin.size() == 1)
	{
		Rmin.replicate(1, D);
		Rmax.replicate(1, D);
	}

	Eigen::MatrixXd mv = 0.2 * (Rmax - Rmin);
	Eigen::MatrixXd Rmin_M = Rmin.replicate(ps, 1); // ά��Ϊps*D��λ�������޺��ٶ�������
	Eigen::MatrixXd Rmax_M = Rmax.replicate(ps, 1);
	Eigen::MatrixXd Vmin_M = -mv.replicate(ps, 1); // ������С����ٶȣ�ֵΪ0.2���ı�����Χ
	Eigen::MatrixXd Vmax_M = -Vmin_M;

	Eigen::VectorXd w = Eigen::VectorXd::LinSpaced(me, 0.9, 0.9) - (Eigen::VectorXd::LinSpaced(me, 1, me) * (0.7 / me));
	double c1 = 0.8; // ����ѧϰ����
	double c2 = 1.49; // Ⱥ��ѧϰ����

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dist(0.0, 1.0);

	Eigen::MatrixXd pos = Rmin_M + (Rmax_M - Rmin_M).cwiseProduct(Eigen::MatrixXd::NullaryExpr(ps, D, [&]() { return dist(gen); }));
	Eigen::MatrixXd vel = Vmin_M + (Vmax_M - Vmin_M).cwiseProduct(Eigen::MatrixXd::NullaryExpr(ps, D, [&]() { return dist(gen); }));

	std::vector<Particle> particles(ps);
	for (int i = 0; i < ps; ++i)
	{
		particles[i].position = pos.row(i);
		particles[i].velocity = vel.row(i);
		particles[i].fitness = fhd(particles[i].position);
		particles[i].pBest = particles[i].position;
		particles[i].pBestFitness = particles[i].fitness;
	}

	double gbestval = particles[0].pBestFitness;
	int minIndex = 0;
	for (int i = 1; i < ps; ++i)
	{
		if (particles[i].pBestFitness < gbestval)
		{
			gbestval = particles[i].pBestFitness;
			minIndex = i;
		}
	}
	Eigen::MatrixXd gbest(1, D);
	gbest = particles[minIndex].pBest;


	for (int i = 1; i < me; ++i)
	{
		for (int k = 0; k < ps; ++k)
		{
			particles[k].velocity = w(i) * particles[k].velocity +
				c1 * (particles[k].pBest - particles[k].position).cwiseProduct(Eigen::MatrixXd::NullaryExpr(1, D, [&]() { return dist(gen); })) +
				c2 * (gbest - particles[k].position).cwiseProduct(Eigen::MatrixXd::NullaryExpr(1, D, [&]() { return dist(gen); }));

			particles[k].velocity = particles[k].velocity.cwiseMax(-mv).cwiseMin(mv);
			particles[k].position = particles[k].position + particles[k].velocity;

			particles[k].position = particles[k].position.cwiseMax(Rmin_M.row(k)).cwiseMin(Rmax_M.row(k));

			particles[k].fitness = fhd(particles[k].position);

			if (particles[k].fitness <= particles[k].pBestFitness)
			{
				particles[k].pBest = particles[k].position;
				particles[k].pBestFitness = particles[k].fitness;
			}

			if (particles[k].pBestFitness < gbestval)
			{
				gbest = particles[k].pBest;
				gbestval = particles[k].pBestFitness;
			}
		}

	}

	return std::make_pair(gbest, gbestval);

}

