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


// Solve AX = XB Problem
// cv::Mat calibrateHandEye(std::vector<cv::Mat> Hgij, std::vector<cv::Mat> Hcij);

enum HandEyeType {
	EyeToHand,
	EyeInHand
};

// calibrate Hand to Eye
//@ vH_robot: robot pose (read from the robot)
//@ vH_mark: mark pose in camera (computed from the camera)
Eigen::Isometry3f calibrateHandEye(
	std::vector<Eigen::Isometry3f> const & vH_robot,
	std::vector<Eigen::Isometry3f> const & vH_mark,
	HandEyeType t = EyeToHand
);

Eigen::Isometry3f sovleAXequalXB(
	std::vector<Eigen::Isometry3f> const & vA,
	std::vector<Eigen::Isometry3f> const & vB
);

}
