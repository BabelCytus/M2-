#ifndef __EchoSequence_H
#define __EchoSequence_H

#include <string>
#include <vector>

#include "EchoPair.h"

class EchoSequence
{
public:
	EchoSequence(const std::string& dirName, int numX, int numY, double fovX, double fovY, double ratioX, double ratioY, 
		double RANSACMinimalSampleSize = 0.2, double RANSACMinimalInliersNumber = 0.5, int RANSACIterationNumber = 1000, double RANSACDistanceThreshold = 0.01, int timeStep = 1);
	~EchoSequence();
	void setNumBlocks(int numX, int numY);
	void setDebugDisplay(bool display);
	void importRayleighDecorrelationCurves(const std::string& axialName, const std::string& lateralName, const std::string& elevationName);
	void computeDisplacements();
	void displayDisplacements();

protected:
	std::vector<EchoPair*> m_EchoPairs;
	std::vector<Eigen::Matrix3d> m_Rotations;
	std::vector<Eigen::Vector3d> m_Translations;
	std::vector<EchoPair::MatrixOutliers> m_OutlierMatrices;

private:
	double m_FOVX{}, m_FOVY{};
	int m_TimeStep{ 1 };
	EchoPair::MatrixCurve m_rhoAxialR, m_rhoLateralR, m_rhoElevationR;
	int m_RANSACMinimalSampleSize{ 10 };		// Minimal size of points in random samples (% of block numbers)
	int m_RANSACMinimalInliersNumber{ 30 };		// Minimal number of inliers to validate a model (% of block numbers)
	int m_RANSACIterationNumber{ 1000 };		// Number of RANSAC iterations
	double m_RANSACDistanceThreshold{ 0.1 };	// Initial distance threshold in mm for inlier selection

	bool isPowerOf2(int n);
	template<class T> 
	std::vector<T> getTokensFromString(std::string& str, int begin = 0, const std::string& delimiter = ";");
};

#endif // !__EchoSequence_H