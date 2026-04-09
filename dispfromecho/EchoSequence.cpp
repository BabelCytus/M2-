#include <filesystem>
#include <set>
#include <cassert>
#include <fstream>
#include <iostream>

#include "EchoSequence.h"

namespace fs = std::filesystem;

#define assertm(exp, msg) assert(((void)msg, exp))

EchoSequence::EchoSequence(const std::string& dirName, int numX, int numY, double fovX, double fovY, double ratioX, double ratioY, 
	double RANSACMinimalSampleSize, double RANSACMinimalInliersNumber, int RANSACIterationNumber, double RANSACDistanceThreshold, int timeStep)
	: m_RANSACMinimalSampleSize((int)(RANSACMinimalSampleSize * numX * numY)),
	  m_RANSACMinimalInliersNumber((int)(RANSACMinimalInliersNumber * numX * numY)), 
	  m_RANSACIterationNumber(RANSACIterationNumber),
	  m_RANSACDistanceThreshold(RANSACDistanceThreshold),
	  m_TimeStep(timeStep)
{
	std::set<std::string> sorted_by_name;

	for (auto& entry : fs::directory_iterator(dirName))
		sorted_by_name.insert(entry.path().string());

	ImageUInt8* img0RGB = new ImageUInt8(std::next(sorted_by_name.begin(), 0)->c_str()), *imgRGB;
	for (int n{timeStep}; n < sorted_by_name.size(); n += timeStep)
	{
		imgRGB = new ImageUInt8(std::next(sorted_by_name.begin(), n)->c_str());
		m_EchoPairs.push_back(new EchoPair(img0RGB, imgRGB, numX, numY, fovX, fovY, ratioX, ratioY));
		img0RGB = imgRGB;
	}
}

EchoSequence::~EchoSequence()
{
	m_EchoPairs[0]->deleteImg1();
	for (int n{ 0 }; n < m_EchoPairs.size() - 1; ++n)
		m_EchoPairs[n]->deleteImg2();
	m_EchoPairs.clear();
}

bool EchoSequence::isPowerOf2(int n)
{
	// & operation between n and n - 1
	int i{ n & (n - 1) };

	// check if n is a power of 2
	return i == 0;
}

void EchoSequence::setNumBlocks(int numX, int numY)
{
	// Test if block size is power of two
	ImageUInt8* img0{ m_EchoPairs[0]->getImg1() };

	double blockWidth{ static_cast<double>(img0->width()) / numX };
	assertm(blockWidth - floor(blockWidth) == 0, "Image width is not a multiple of the block number!");
	assertm(isPowerOf2(static_cast<int>(blockWidth)), "Block width is not a power of two!");

	double blockHeight{ static_cast<double>(img0->height()) / numY };
	assertm(blockHeight - floor(blockHeight) == 0, "Image height is not a multiple of the block number!");
	assertm(isPowerOf2(static_cast<int>(blockHeight)), "Block height is not a power of two!");

	for (auto echoPair : m_EchoPairs)
		echoPair->setNumBlocks(numX, numY);
}

void EchoSequence::setDebugDisplay(bool display)
{
	for (auto echoPair : m_EchoPairs)
		echoPair->setDebugDisplay(display);
}

void EchoSequence::computeDisplacements()
{
	const double pi{ std::acos(-1.0) };
	int n{};
	Eigen::Vector3d tMoy{ Eigen::Vector3d::Zero() };
	for (auto echoPair : m_EchoPairs)
	{
		auto [R, t] { echoPair->computeRigidTransform(m_rhoAxialR, m_rhoLateralR, m_rhoElevationR, 
			m_RANSACMinimalSampleSize, m_RANSACMinimalInliersNumber, m_RANSACIterationNumber, m_RANSACDistanceThreshold) };

		m_Rotations.push_back(R);
		Eigen::AngleAxisd axisAngle(R);
		std::cout << "Rotation between frame #" << n * m_TimeStep << " and frame #" << (n + 1) * m_TimeStep << "\n";
		std::cout << R << "\n";
		std::cout << "\taxis: " << axisAngle.axis().transpose() << "\n";
		std::cout << "\tangle: " << axisAngle.angle() * 180. / pi << "\n";

		m_Translations.push_back(t);
		std::cout << "Tranlation between frame #" << n * m_TimeStep << " and frame #" << (n + 1) * m_TimeStep << "\n";
		std::cout << t << "\n";

		if (R != Eigen::Matrix3d::Identity())
		{
			tMoy += t;
			++n;
		}
	}
	tMoy /= n;
	std::cout << "Average translation\n" << tMoy << "\n";
}

void EchoSequence::displayDisplacements()
{
	for (auto echoPair : m_EchoPairs)
		echoPair->getDisplay().display();
}

template<class T> 
std::vector<T> EchoSequence::getTokensFromString(std::string& str, int begin, const std::string& delimiter)
{
	std::vector<T> tokens;
	size_t pos = 0;
	std::string token;
	int nToken{};
	while ((pos = str.find(delimiter)) != std::string::npos) {
		token = str.substr(0, pos);
		if (nToken >= begin)
			tokens.push_back(static_cast<T>(std::stod(token)));
		str.erase(0, pos + delimiter.length());
		nToken++;
	}
	tokens.push_back(static_cast<T>(std::stod(str)));
	return tokens;
}

void EchoSequence::importRayleighDecorrelationCurves(const std::string& axialName, const std::string& lateralName, const std::string& elevationName)
{
	int numBlockX{ m_EchoPairs[0]->getNumBlockX() }, numBlockY{ m_EchoPairs[0]->getNumBlockY() };

	// Axial decorrelation curve
	std::ifstream axialCSVFile(axialName);
	std::string line;
	std::getline(axialCSVFile, line);
	std::getline(axialCSVFile, line);
	std::getline(axialCSVFile, line);
	auto axialCoordinates{ getTokensFromString<double>(line, 3) };
	std::getline(axialCSVFile, line);
	std::getline(axialCSVFile, line);
	m_rhoAxialR.resize(numBlockX, std::vector<EchoPair::Curve>(numBlockY));
	while (std::getline(axialCSVFile, line))
	{
		auto correlations{ getTokensFromString<double>(line, 1) };
		int i{ static_cast<int>(correlations[0]) }, j{ static_cast<int>(correlations[1]) };
		for (int k{}; k < correlations.size() - 2; ++k)
			m_rhoAxialR[i][j].push_back({ axialCoordinates[k], correlations[k + 2] });
	}

	// Lateral decorrelation curve
	std::ifstream lateralCSVFile(lateralName);
	std::getline(lateralCSVFile, line);
	std::getline(lateralCSVFile, line);
	std::getline(lateralCSVFile, line);
	auto lateralCoordinates{ getTokensFromString<double>(line, 3) };
	std::getline(lateralCSVFile, line);
	std::getline(lateralCSVFile, line);
	m_rhoLateralR.resize(numBlockX, std::vector<EchoPair::Curve>(numBlockY));
	while (std::getline(lateralCSVFile, line))
	{
		auto correlations{ getTokensFromString<double>(line, 1) };
		int i{ static_cast<int>(correlations[0]) }, j{ static_cast<int>(correlations[1]) };
		for (int k{}; k < correlations.size() - 2; ++k)
			m_rhoLateralR[i][j].push_back({ lateralCoordinates[k], correlations[k + 2] });
	}

	// Elevation decorrelation curve
	std::ifstream elevationCSVFile(elevationName);
	std::getline(elevationCSVFile, line);
	std::getline(elevationCSVFile, line);
	std::getline(elevationCSVFile, line);
	auto elevationCoordinates{ getTokensFromString<double>(line, 3) };
	std::getline(elevationCSVFile, line);
	std::getline(elevationCSVFile, line);
	m_rhoElevationR.resize(numBlockX, std::vector<EchoPair::Curve>(numBlockY));
	while (std::getline(elevationCSVFile, line))
	{
		auto correlations{ getTokensFromString<double>(line, 1) };
		int i{ static_cast<int>(correlations[0]) }, j{ static_cast<int>(correlations[1]) };
		for (int k{}; k < correlations.size() - 2; ++k)
			m_rhoElevationR[i][j].push_back({ elevationCoordinates[k], correlations[k + 2] });
	}
}