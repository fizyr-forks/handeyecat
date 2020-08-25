/***************************************
*
* LanXin TECH, All Rights Reserverd.
* Created at Thu May 16 10:27:58 2019
* Contributor: Ling Shi, Ph.D
* Email: lshi@robvision.cn
*
***************************************/

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <cmath>

namespace lanXin {

// hand eye calibration
// ---------------------------------------------

#define XEPS					1e-6

template<typename T>
bool non_zero(const T& val)
{
	return std::abs(val) > 1e-6;
}


// Rodrigues transformation
Eigen::Vector3f rodrigues2(const Eigen::Matrix3f& matrix);


inline Eigen::Matrix3f fromEulers(float rx, float ry = .0f, float rz = .0f)
{
	Eigen::AngleAxisf quat = Eigen::AngleAxisf(rx, Eigen::Vector3f::UnitX());
	if (non_zero(ry))
	{
		quat = Eigen::AngleAxisf(ry, Eigen::Vector3f::UnitY()) * quat;
	}
	if(non_zero(rz))
		quat = Eigen::AngleAxisf(rz, Eigen::Vector3f::UnitZ()) * quat;
	return quat.matrix();
}

inline Eigen::Isometry3f getTransM(Eigen::Vector3f t, Eigen::Vector3f eulers)
{
	Eigen::Isometry3f H;
	H.setIdentity();
	H.linear() = fromEulers(eulers[0], eulers[1], eulers[2]);
	H.translation() = t;
	return H;
}


// Solve AX = XB Problem
// cv::Mat calibrateHandEye(std::vector<cv::Mat> Hgij, std::vector<cv::Mat> Hcij);

enum HandEyeType
{
	EyeToHand,
	EyeInHand
};

// calibrate Hand to Eye
//@ vH_robot: robot pose (read from the robot)
//@ vH_mark: mark pose in camera (computed from the camera)
Eigen::Isometry3f calibrateHandEye(std::vector<Eigen::Isometry3f>& vH_robot, std::vector<Eigen::Isometry3f>& vH_mark, HandEyeType t = EyeToHand);

Eigen::Isometry3f sovleAXequalXB(std::vector<Eigen::Isometry3f>& vA, std::vector<Eigen::Isometry3f>& vB);

} /* End of namespace lanXin */
