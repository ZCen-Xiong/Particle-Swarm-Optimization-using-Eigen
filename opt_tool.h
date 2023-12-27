#pragma once

#pragma once
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include <functional>
#include <Eigen/Dense>
#include <random>

// struct Particle;
// std::pair<Eigen::MatrixXd, double> CLPSO(const std::function<double(const Eigen::MatrixXd&)>& fhd, int Dimension, const Eigen::MatrixXd& Rmin, const Eigen::MatrixXd& Rmax, int Max_Gen, int Particle_Number);
std::pair<Eigen::MatrixXd, double> CLPSO(int Dimension, const Eigen::MatrixXd& Rmin, const Eigen::MatrixXd& Rmax, int Max_Gen, int Particle_Number, const Eigen::VectorXd& x0, double max_TOF, const Eigen::VectorXd& Xf, double h);
