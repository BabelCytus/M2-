#ifndef __RigidTransform_H
#define __RigidTransform_H

#include <vector>

#include <Eigen/Dense>

using Indices = std::vector<unsigned int>;
using Points = std::vector<Eigen::Vector3d>;
using RigidTransform = std::tuple<Eigen::Matrix3d, Eigen::Vector3d>;

const double pi{ std::acos(-1.0) };

void Rotation(const Points& points, Eigen::Vector3d center, Eigen::Vector3d axis, double angle, Points& newPoints);
void Translation(const Points& points, const Eigen::Vector3d& t, Points& newPoints);
RigidTransform Procrustes(const Points& pointsA, const Points& pointsB);
std::tuple<RigidTransform, Indices, double> RANSAC(const Points& modelPoints, const Points& dataPoints, unsigned int sampleMinimalSize, unsigned int nIteration, double distanceThreshold, unsigned int nInliers);


#endif // !__RigidTransform_H