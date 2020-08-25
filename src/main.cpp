/***************************************
*
* LanXin TECH, All Rights Reserverd.
* Created at Thu May 16 10:27:58 2019
* Contributor: Ling Shi, Ph.D
* Email: lshi@robvision.cn
*
***************************************/

#include "common.h"

#include <iostream>

int main() {
	Eigen::Vector3f O_world(0, 0, 0);

	// camera in world, unknown, to calibrate
	// transform camera to world, to estimate
	Eigen::Isometry3f H_c_in_w = lanXin::getTransM(Eigen::Vector3f(0.8, 0.12, 2.4), Eigen::Vector3f(M_PI*1.02, -0.12*M_PI, M_PI *0.51));

	std::cout << "camera in world:\n" << H_c_in_w.matrix() << "\n";

	// grid in end, unkonwn, but we don't not have to solve it
	Eigen::Isometry3f H_g_in_e = lanXin::getTransM(Eigen::Vector3f(0.10, 0.20, -0.401), Eigen::Vector3f(0.20, -0.20, 0.20));

	// Grid in world pose list
	std::vector<Eigen::Isometry3f> vH_e_in_w, vH_g_in_w;
	vH_e_in_w.push_back(lanXin::getTransM(Eigen::Vector3f(0.56, 0.4, 0.3), Eigen::Vector3f(0.2, 0.2, 0)));
	vH_e_in_w.push_back(lanXin::getTransM(Eigen::Vector3f(0.6, -0.14, 0.3), Eigen::Vector3f(-0.3, -0.17, 0)));
	vH_e_in_w.push_back(lanXin::getTransM(Eigen::Vector3f(0.98, 0.4, 0.40), Eigen::Vector3f(0.2, 0.15, 0.1)));
	vH_e_in_w.push_back(lanXin::getTransM(Eigen::Vector3f(0.85, -0.2, 0.32), Eigen::Vector3f(-0.19, 0.2, 0)));

	int n = vH_e_in_w.size();
	for (int i = 0; i < n; ++i) {
		vH_g_in_w.push_back(vH_e_in_w[i] * H_g_in_e);
	}

	// grid in camera = H_c_in_w^(-1) * H_g_in_w;
	std::vector<Eigen::Isometry3f> vH_g_in_c;
	for (auto it = vH_g_in_w.begin(); it!= vH_g_in_w.end(); ++it) {
		vH_g_in_c.push_back(H_c_in_w.inverse() * (*it));
	}

	// add random noise to the input data
	for (int i = 0; i < n; ++i) {
		// cout << rng.uniform(-0.01, 0.01) << endl;

		// Eigen::Vector3f v = Eigen::Vector3f(rng.uniform(-0.01, 0.01), rng.uniform(-0.01, 0.01), rng.uniform(-0.01, 0.01));

		// vH_g_in_c[i].translation() += v/102;

		//Eigen::Vector3f eulers = vH_g_in_c[i].linear().eulerAngles(2, 1, 0);
		// v = Eigen::Vector3f(rng.uniform(-0.01, 0.01), rng.uniform(-0.01, 0.01), rng.uniform(-0.01, 0.01)) /100;

		//cout <<"1 "<< vH_g_in_c[i].linear() << endl;
		//cout <<"2 "<<  fromEulers(eulers[2], eulers[1], eulers[0]) << endl;
		// vH_g_in_c[i].linear() = vH_g_in_c[i].linear() * fromEulers(v[0], v[1], v[2]);
	}


	// compare
	for (std::size_t i = 0; i < vH_g_in_c.size(); ++i) {
		Eigen::Isometry3f result = H_c_in_w * vH_g_in_c[i];
		Eigen::Isometry3f g_in_w = vH_g_in_w[i];
		std::cout << "Compare: --- At " << i << "\n" << result.matrix() - g_in_w.matrix() << "\n";
	}

	// compose A and B
	Eigen::Isometry3f H = lanXin::calibrateHandEye(vH_e_in_w, vH_g_in_c);
	std::cout << "H:\n" << H.matrix() << "\n";

	return 0;
}
