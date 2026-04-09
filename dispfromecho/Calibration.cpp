#include <fstream>

#include "Calibration.h"

Calibration::Calibration(const std::string& dirName, int numX, int numY, double fovX, double fovY, const std::string& probeName, double elevationStep, double ratioX, double ratioY, int timeStep)
	: EchoSequence(dirName, numX, numY, fovX, fovY, ratioX, ratioY, timeStep), 
		m_ElevationStep(elevationStep), m_ProbeName(probeName)
{

}

// Compute axial, lateral and elevation decorrelation curves for each block
// Axial and lateral decorrelations are computed by the average of autocorrelation 
//		for all elevation positions
// Elevation decorrelation is computed by cross-correlation 
//		between each elevation frame and frame at zero elevation
void Calibration::computeCalibrationDecorrelationCurves()
{
	int numBlockX{ m_EchoPairs[0]->getNumBlockX() }, numBlockY{ m_EchoPairs[0]->getNumBlockY() };

	// Axial and lateral decorrelation curves
	m_AxialDecorrelationCurves.resize(numBlockX, std::vector<EchoPair::Curve>(numBlockY));
	m_LateralDecorrelationCurves.resize(numBlockX, std::vector<EchoPair::Curve>(numBlockY));
	double pixelSizeX{ m_EchoPairs[0]->getPixelSizeX() }, pixelSizeY{ m_EchoPairs[0]->getPixelSizeY() };
	for (int i{}; i < numBlockX; ++i)
	{
		for (int j{}; j < numBlockY; ++j)
		{
			ImageDouble moyAutocorrelation{ m_EchoPairs[0]->getBlockAutocorrelation(i, j) };
			for (int k{ 1 }; k < m_EchoPairs.size(); ++k)
				moyAutocorrelation += m_EchoPairs[k]->getBlockAutocorrelation(i, j);
			moyAutocorrelation /= m_EchoPairs.size();
			//moyAutocorrelation.display();
			m_AxialDecorrelationCurves[i][j] = EchoPair::extractCorrelationCurve(moyAutocorrelation, pixelSizeY, EchoPair::axial);
			m_LateralDecorrelationCurves[i][j] = EchoPair::extractCorrelationCurve(moyAutocorrelation, pixelSizeY, EchoPair::lateral);
		}
	}

	// Elevation decorrelation curves
	m_ElevationDecorrelationCurves.resize(numBlockX, std::vector<EchoPair::Curve>(numBlockY));
	for (int i{}; i < numBlockX; ++i)
	{
		for (int j{}; j < numBlockY; ++j)
		{
			ImageDouble imgRef{ m_EchoPairs[0]->getImg1Block(i, j) };
			int width_2{ imgRef.width() / 2 }, height_2{ imgRef.height() / 2 };
			m_ElevationDecorrelationCurves[i][j].push_back(std::make_pair(0.0, 1.0));
			double z{ m_ElevationStep };
			for (int k{ 1 }; k < m_EchoPairs.size(); ++k)
			{
				ImageDouble img{ m_EchoPairs[k]->getImg1Block(i, j) };
				// If for some reason image pairs are the same, cross-correlation will fail
				if (imgRef == img)
				{
					m_ElevationDecorrelationCurves[i][j].push_back(std::make_pair(z, 1.0));
				}
				else
				{
					ImageDouble crossCorrelationImage{ normalized_correlation(imgRef, img) };
					double crossCorrelation{ crossCorrelationImage(width_2, height_2) };
					m_ElevationDecorrelationCurves[i][j].push_back(std::make_pair(z, crossCorrelation));
				}
				z += m_ElevationStep;
			}
		}
	}
}

void Calibration::exportCalibrationDecorrelationCurves()
{
	int numBlockX{ m_EchoPairs[0]->getNumBlockX() }, numBlockY{ m_EchoPairs[0]->getNumBlockY() };

	// Axial decorrelation curve
	std::string axialCSVFileName{ m_ProbeName + "_axial_"
		+ std::to_string(numBlockX) + "x" + std::to_string(numBlockY) + ".csv"};
	std::ofstream axialCSVFile(axialCSVFileName);
	axialCSVFile	<< "Axial decorrelation curve for " << m_ProbeName << " probe "
					<< " with (" + std::to_string(numBlockX) + ", " + std::to_string(numBlockX) + ") blocks\n";
	axialCSVFile << "\nAxial coordinate in mm;;";
	for (const auto& decorrelation : m_AxialDecorrelationCurves[0][0])
		axialCSVFile << ";" << decorrelation.first;
	axialCSVFile << "\n\nBlock;X;Y\n";
	for (int i{}; i < numBlockX; ++i)
	{
		for (int j{}; j < numBlockY; ++j)
		{
			axialCSVFile << ";" << i << ";" << j;
			for (const auto& decorrelation : m_AxialDecorrelationCurves[i][j])
				axialCSVFile << ";" << decorrelation.second;
			axialCSVFile << "\n";
		}
	}
	axialCSVFile.close();

	// Lateral decorrelation curves
	std::string lateralCSVFileName{ m_ProbeName + "_lateral_"
		+ std::to_string(numBlockX) + "x" + std::to_string(numBlockY) + ".csv"};
	std::ofstream lateralCSVFile(lateralCSVFileName);
	lateralCSVFile << "Lateral decorrelation curve for " << m_ProbeName << " probe "
		<< " with (" + std::to_string(numBlockX) + ", " + std::to_string(numBlockX) + ") blocks\n";
	lateralCSVFile << "\nLateral coordinate in mm;;";
	for (const auto& decorrelation : m_LateralDecorrelationCurves[0][0])
		lateralCSVFile << ";" << decorrelation.first;
	lateralCSVFile << "\n\nBlock;X;Y\n";
	for (int i{}; i < numBlockX; ++i)
	{
		for (int j{}; j < numBlockY; ++j)
		{
			lateralCSVFile << ";" << i << ";" << j;
			for (const auto& decorrelation : m_LateralDecorrelationCurves[i][j])
				lateralCSVFile << ";" << decorrelation.second;
			lateralCSVFile << "\n";
		}
	}
	lateralCSVFile.close();

	// Elevation decorrelation curves
	std::string elevationCSVFileName{ m_ProbeName + "_elevation_"
		+ std::to_string(numBlockX) + "x" + std::to_string(numBlockY) + ".csv"};
	std::ofstream elevationCSVFile(elevationCSVFileName);
	elevationCSVFile << "Elevation decorrelation curve for " << m_ProbeName << " probe "
		<< " with (" + std::to_string(numBlockX) + ", " + std::to_string(numBlockX) + ") blocks\n";
	elevationCSVFile << "\nElevation coordinate in mm;;";
	for (const auto& decorrelation : m_ElevationDecorrelationCurves[0][0])
		elevationCSVFile << ";" << decorrelation.first;
	elevationCSVFile << "\n\nBlock;X;Y\n";
	for (int i{}; i < numBlockX; ++i)
	{
		for (int j{}; j < numBlockY; ++j)
		{
			elevationCSVFile << ";" << i << ";" << j;
			for (const auto& decorrelation : m_ElevationDecorrelationCurves[i][j])
				elevationCSVFile << ";" << decorrelation.second;
			elevationCSVFile << "\n";
		}
	}
	elevationCSVFile.close();
}