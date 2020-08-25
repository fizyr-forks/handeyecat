/***************************************
*
* LanXin TECH, All Rights Reserverd.
* Created at Thu May 16 10:27:58 2019
* Contributor: Ling Shi, Ph.D
* Email: lshi@robvision.cn
*
***************************************/

#include "common.h"

#include <Eigen/SVD>

#include <stdexcept>

namespace lanXin {

namespace {
	constexpr double epsilon = 1e-6;

	Eigen::Vector3f rodrigues2(Eigen::Matrix3f const & matrix) {
		Eigen::JacobiSVD<Eigen::Matrix3f> svd(matrix, Eigen::ComputeFullV | Eigen::ComputeFullU);
		Eigen::Matrix3f R = svd.matrixU() * svd.matrixV().transpose();

		double rx = R(2, 1) - R(1, 2);
		double ry = R(0, 2) - R(2, 0);
		double rz = R(1, 0) - R(0, 1);

		double s = sqrt((rx*rx + ry*ry + rz*rz)*0.25);
		double c = (R.trace() - 1) * 0.5;
		c = std::max(-1.0, std::min(1.0, c));

		double theta = acos(c);

		if (s < epsilon) {
			if (c > 0) {
				rx = ry = rz = 0;
			} else {
				rx = sqrt(std::max((R(0, 0) + 1) * 0.5, 0.0));
				ry = sqrt(std::max((R(1, 1) + 1) * 0.5, 0.0)) * (R(0, 1) < 0 ? -1.0 : 1.0);
				rz = sqrt(std::max((R(2, 2) + 1) * 0.5, 0.0)) * (R(0, 2) < 0 ? -1.0 : 1.0);

				if (fabs(rx) < fabs(ry) && fabs(rx) < fabs(rz) && (R(1, 2) > 0) != (ry*rz > 0)) {
					rz = -rz;
				}
				theta /= sqrt(rx*rx + ry*ry + rz*rz);
				rx *= theta;
				ry *= theta;
				rz *= theta;
			}
		} else {
			double vth = 1 / (2 * s);
			vth *= theta;
			rx *= vth; ry *= vth; rz *= vth;
		}
		return Eigen::Vector3f(rx, ry, rz);
	}

	Eigen::Matrix3f skew(Eigen::Vector3f v) {
		Eigen::Matrix3f rot = Eigen::Matrix3f::Zero();

		rot(0, 1) = -v(2);
		rot(0, 2) = v(1);
		rot(1, 2) = -v(0);

		rot(1, 0) = -rot(0, 1);
		rot(2, 0) = -rot(0, 2);
		rot(2, 1) = -rot(1, 2);

		return rot;
	}
}

Eigen::Isometry3f calibrateHandEye(
	std::vector<Eigen::Isometry3f> & vH_robot,
	std::vector<Eigen::Isometry3f> & vH_mark,
	HandEyeType t
) {
	const int n = std::min(vH_robot.size(), vH_mark.size());
	if (vH_robot.size() != vH_mark.size()) {
		throw std::runtime_error("calibrateHandEye requires equal number of robot and mark poses, got " + std::to_string(vH_robot.size()) + " and " + std::to_string(vH_mark.size()));
	}
	if (n < 3) {
		throw std::runtime_error("calibrateHandEye requires at-least 3 point pars, got only " + std::to_string(n) + " pairs");
	}

	std::vector<Eigen::Isometry3f> vA;
	std::vector<Eigen::Isometry3f> vB;
	vA.reserve(n);
	vB.reserve(n);

	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; ++j) {
			//if(i == 0 && j==i+1) continue;
			if (t == EyeToHand) {
				Eigen::Isometry3f A = vH_robot[j] * vH_robot[i].inverse();
				Eigen::Isometry3f B = vH_mark[j] * vH_mark[i].inverse();

				vA.push_back(A);
				vB.push_back(B);
			} else if (t == EyeInHand) {
				Eigen::Isometry3f A = vH_robot[j].inverse() * vH_robot[i];
				Eigen::Isometry3f B = vH_mark[j] * vH_mark[i].inverse();

				vA.push_back(A);
				vB.push_back(B);
			}
		}
	}

	return sovleAXequalXB(vA, vB);
}

Eigen::MatrixXf svdInverse(Eigen::MatrixXf  A) {
	Eigen::JacobiSVD<Eigen::MatrixXf> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV); // M = USV*
	float pinvtoler = 1.e-6; //tolerance
	int row = A.rows();
	int col = A.cols();
	int k = std::min(row, col);
	Eigen::MatrixXf singularValues_inv = svd.singularValues(); //奇异值
	Eigen::MatrixXf singularValues_inv_mat = Eigen::MatrixXf::Zero(col, row);

	for (long i = 0; i < k; ++i) {
		if (singularValues_inv(i) > pinvtoler) {
			singularValues_inv(i) = 1.0 / singularValues_inv(i);
		} else {
			singularValues_inv(i) = 0;
		}
	}

	for (long i = 0; i < k; ++i) {
		singularValues_inv_mat(i, i) = singularValues_inv(i);
	}

	return svd.matrixV() * singularValues_inv_mat * svd.matrixU().transpose(); //X=VS+U*
}

Eigen::Isometry3f sovleAXequalXB(std::vector<Eigen::Isometry3f>& vA, std::vector<Eigen::Isometry3f>& vB) {
	if (vA.size() != vB.size()) {
		throw std::runtime_error("sovleAXequalXB requires input vectors to have the same size, got " + std::to_string(vA.size()) + " and " + std::to_string(vB.size()));
	}

	const int n = vA.size();

	Eigen::MatrixXf A = Eigen::MatrixXf::Zero(n * 3, 3);
	Eigen::MatrixXf b = Eigen::MatrixXf::Zero(n * 3, 1);

	for (int i = 0; i < n; ++i) {
		Eigen::Matrix3f R_a = vA[i].linear();
		Eigen::Matrix3f R_b = vB[i].linear();

		Eigen::Vector3f rod_a = rodrigues2(R_a);
		Eigen::Vector3f rod_b = rodrigues2(R_b);

		float theta_a = rod_a.norm();
		float theta_b = rod_b.norm();

		rod_a /= theta_a;
		rod_b /= theta_b;

		Eigen::Vector3f P_a = 2 * sin(theta_a / 2) * rod_a;
		Eigen::Vector3f P_b = 2 * sin(theta_b / 2) * rod_b;

		A.middleRows(3 * i, 3) = skew(P_b + P_a);
		b.middleRows(3 * i, 3) = P_b - P_a;
	}

	// 3 by 3*n
	Eigen::MatrixXf pinA = svdInverse(A);

	// 3 by 1 = 3 by 3*n multi 3*n by 1
	Eigen::Vector3f H_ba_prime = pinA * b;

	Eigen::Vector3f H_ba = 2 * H_ba_prime / sqrt(1 + std::pow(H_ba_prime.norm(), 2));

	// 1 by 3
	Eigen::MatrixXf H_ba_Trs = H_ba.transpose();

	Eigen::Matrix3f R_ba = (1 - std::pow(H_ba.norm(), 2) / 2) * Eigen::Matrix3f::Identity()
		+ 0.5 * (H_ba * H_ba_Trs + sqrt(4 - std::pow(H_ba.norm(), 2)) * skew(H_ba));

	A.setZero();
	b.setZero();
	for (int i = 0; i < n; ++i) {
		Eigen::Matrix3f AA = vA[i].linear() - Eigen::Matrix3f::Identity();
		Eigen::Vector3f bb = R_ba * vB[i].translation() - vA[i].translation();
		//for (int row = 3 * i; row < 3 * i + 3; ++row)
		//{
		//	for (int col = 0; col < 3; ++col)
		//	{
		//		A(row, col) = AA(row - 3 * i, col);
		//	}

		//	b(row) = bb(row - 3 * i);
		//}
		A.middleRows(3 * i, 3) = AA;
		b.middleRows(3 * i, 3) = bb;
	}
	pinA = svdInverse(A);
	Eigen::Vector3f t_ba = pinA * b;

	Eigen::Isometry3f H = Eigen::Isometry3f::Identity();
	H.linear() = R_ba;
	H.translation() = t_ba;

	// check
	for (int i = 0; i < n; ++i) {
		// GeoTransform AX = vA[i] * H;
		// GeoTransform XB = H * vB[i];

		// Eigen::Vector3f angles1 = AX.linear().eulerAngles(0, 1, 2);
		// Eigen::Vector3f angles2 = XB.linear().eulerAngles(0, 1, 2);
		// cout << i << " Dist Error: " << (AX.translation() - XB.translation()).norm()
		// 	<< ". Rotation Error: " << (angles1 - angles2).transpose() << endl;
	}

	return H;
}

}
