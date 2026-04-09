#include <sstream>
#include <iomanip>

#include "Correlation.h"
#include "RigidTransform.h"

#include "EchoPair.h"

EchoPair::EchoPair(ImageUInt8 *img1, ImageUInt8 *img2, int numX, int numY, double fovX, double fovY, double ratioX, double ratioY)
	: m_Img1(img1), m_Img2(img2), m_Axial2ElevationResolutionRatio(ratioX), m_Lateral2ElevationResolutionRatio(ratioY)
{
	setNumBlocks(numX, numY);
	setPixelSize(fovX, fovY);
	blockDecomposition();
}

EchoPair::~EchoPair()
{
	if (m_Img1)
	{
		delete m_Img1;
		m_Img1 = nullptr;
	}
	if (m_Img2)
	{
		delete m_Img2;
		m_Img2 = nullptr;
	}
}

void EchoPair::setPixelSize(double fovX, double fovY)
{
	m_PixelSizeX = fovX / m_Img1->width();
	m_PixelSizeY = fovY / m_Img1->height();
}

void EchoPair::setNumBlocks(int numX, int numY)
{
	m_NumBlockX = numX;
	m_NumBlockY = numY;
	m_BlockWidth = m_Img1->width() / numX;
	m_BlockHeight = m_Img2->height() / numY;

	deleteBlocks();
	blockDecomposition();
}

void EchoPair::blockDecomposition()
{
	m_BlockMatrix1.resize(m_NumBlockX, std::vector<ImageDouble>(m_NumBlockY));
	m_BlockMatrix2.resize(m_NumBlockX, std::vector<ImageDouble>(m_NumBlockY));
	for (int i{}; i < m_NumBlockX; ++i)
		for (int j{}; j < m_NumBlockY; ++j)
		{
			int x0{ i * m_BlockWidth }, x1{ x0 + m_BlockWidth - 1 };
			int y0{ j * m_BlockHeight }, y1{ y0 + m_BlockHeight - 1};
			m_BlockMatrix1[i][j] = ImageDouble(m_Img1->get_channel(0).get_crop(x0, y0, x1, y1));
			m_BlockMatrix2[i][j] = ImageDouble(m_Img2->get_channel(0).get_crop(x0, y0, x1, y1));
		}
}

void EchoPair::deleteBlocks()
{
	if (m_BlockMatrix1.size() != 0)
	{
		for (auto blockColumn : m_BlockMatrix1)
		{
			blockColumn.clear();
		}
		m_BlockMatrix1.clear();
	}

	if (m_BlockMatrix2.size() != 0)
	{
		for (auto blockColumn : m_BlockMatrix2)
		{
			blockColumn.clear();
		}
		m_BlockMatrix2.clear();
	}
}

EchoPair::Curve EchoPair::extractCorrelationCurve(const ImageDouble& correlationImage, double pixelSize, CorrelationType type)
{
	Curve Curve;
	if (type == lateral)
	{
		double x{};
		for (int i{ correlationImage.width() / 2 }; i < correlationImage.width(); ++i)
		{
			Curve.push_back(std::make_pair(x, correlationImage(i, correlationImage.height() / 2)));
			x += pixelSize;
		}
	}
	else if (type == axial)
	{
		double y{};
		for (int j{ correlationImage.height() / 2 }; j < correlationImage.height(); ++j)
		{
			Curve.push_back(std::make_pair(y, correlationImage(correlationImage.width() / 2, j)));
			y += pixelSize;
		}
	}
	return Curve;
}

ImageDouble EchoPair::getTranslatedBlock(const Vec2& translation, int i, int j)
{
	int x0{ i * m_BlockWidth + (int)translation[0] }, x1{x0 + m_BlockWidth - 1};
	int y0{ j * m_BlockHeight + (int)translation[1] }, y1{y0 + m_BlockHeight - 1};
	return ImageDouble(m_Img1->get_channel(0).get_crop(x0, y0, x1, y1));
}

ImageDouble EchoPair::getBlockAutocorrelation(int i, int j)
{
	ImageDouble correlationImage{ normalized_correlation(m_BlockMatrix1[i][j], m_BlockMatrix1[i][j]) };
	return correlationImage;
}

void EchoPair::computeInplaneDisplacement()
{
	m_InplaneDisp.resize(m_NumBlockX, std::vector<Vec2>(m_NumBlockY));
	for (int i{}; i < m_NumBlockX; ++i)
	{
		for (int j{}; j < m_NumBlockY; ++j)
		{
			ImageDouble correlationImage{ normalized_correlation(m_BlockMatrix1[i][j], m_BlockMatrix2[i][j]) };
			translation(correlationImage, m_InplaneDisp[i][j][0], m_InplaneDisp[i][j][1]);
			//m_BlockMatrix1[i][j].display();
			//m_BlockMatrix2[i][j].display();
			//correlationImage.display();
		}
	}
}

EchoPair::Curve EchoPair::computeCoherentScatteringAmount(const Curve& rhoR, const Curve& rho)
{
	Curve kCurve;
	for (Curve::const_iterator itR{ rhoR.cbegin() }, it{ rho.cbegin() }; it < rho.cend(); ++itR, ++it)
	{
		// If decorrelation equals Rayleigh, k = 0
		if (it->second <= itR->second) kCurve.push_back({ it->first, 0.0 });
		else
		{
			auto k{ ((1 - it->second * itR->second)
					- std::sqrt((1 - it->second * it->second) * (1 - itR->second * itR->second)))
					/ (it->second - itR->second) };
			kCurve.push_back({ it->first, k });
		}
	}
	return kCurve;
}

double EchoPair::abscissaLinearInterpolateAt(const Curve& curve, double y)
{
	int i{};
	while (i < curve.size() && curve[i].second >= y) ++i;
	if (i == curve.size()) --i;

	// x is linearly interpolated between curve[i - 1] and curve[i]
	return curve[i - 1].first
		+ (curve[i].first - curve[i - 1].first)
		* (y - curve[i - 1].second) / (curve[i].second - curve[i - 1].second);
}

double EchoPair::ordinateLinearInterpolateAt(const Curve& curve, double x)
{
	int i{};
	while (i < curve.size() && curve[i].first <= x) ++i;
	if (i == curve.size()) --i;

	// y is linearly interpolated between curve[i - 1] and curve[i]
	return curve[i - 1].second
		+ (curve[i].second - curve[i - 1].second)
		* (x - curve[i - 1].first) / (curve[i].first - curve[i - 1].first);
}

EchoPair::Curve EchoPair::averageBetween2CurvesWithScaling(const Curve& curve1, double ratio1, const Curve& curve2, double ratio2)
{
	Curve average;
	// The average is computed by interpolating the longest curve 
	// at the positions of the shortest one
	if (curve1.back().first * ratio1 < curve2.back().first * ratio2)
	{
		for (const auto& point1 : curve1)
		{
			auto value2{ ordinateLinearInterpolateAt(curve2, point1.first * ratio1 / ratio2) };
			average.push_back({ point1.first * ratio1, 0.5 * (point1.second + value2) });
		}
	}
	else
	{
		for (const auto& point2 : curve2)
		{
			auto value1{ ordinateLinearInterpolateAt(curve1, point2.first * ratio2 / ratio1) };
			average.push_back({ point2.first * ratio2, 0.5 * (value1 + point2.second) });
		}

	}
	return average;
}

EchoPair::Curve EchoPair::resampleCurve(const Curve& curve, const Curve& curveRef)
{
	Curve resampled;
	for (const auto& pointRef : curveRef)
	{
		if (pointRef.first > curve.back().first) break;
		auto value{ ordinateLinearInterpolateAt(curve, pointRef.first) };
		resampled.push_back({ pointRef.first, value });
	}
	return resampled;
}

EchoPair::Curve EchoPair::computeElevationDecorrelationCurve(const Curve& k, const Curve& rhoR)
{
	Curve rhoElevationCurve;
	for (Curve::const_iterator itK{ k.cbegin() }, itRhoR{ rhoR.cbegin() }; itK < k.cend(); ++itK, ++itRhoR)
	{
		double tmp{ 1 - itK->second * itK->second };
		double rhoElevation{ (tmp * itRhoR->second + 2 * itK->second) / (tmp + 2 * itK->second * itRhoR->second) };
		rhoElevationCurve.push_back({ itK->first, rhoElevation });
	}
	return rhoElevationCurve;
}

void EchoPair::computeOutplaneDisplacement(const MatrixCurve& rhoAxialR, const MatrixCurve& rhoLateralR, const MatrixCurve& rhoElevationR)
{
	m_OutplaneDisp.resize(m_NumBlockX, std::vector<double>(m_NumBlockY));

	// All image blocks
	for (int i{}; i < m_NumBlockX; ++i)
	{
		for (int j{}; j < m_NumBlockY; ++j)
		{
			// Compute first image autocorrelation
			auto autocorrelation{ getBlockAutocorrelation(i, j) };
			if (m_DebugDisplay) autocorrelation.display("First block autocorrelation");

			// Compute axial correlation curve
			auto rhoAxial{ extractCorrelationCurve(autocorrelation, m_PixelSizeX, axial) };
			if (m_DebugDisplay) displayGraph(rhoAxial, "Axial correlation curve");

			// Compute axial coherent scattering amount
			auto kAxial{ computeCoherentScatteringAmount(rhoAxialR[i][j], rhoAxial) };
			if (m_DebugDisplay) displayGraph(rhoAxialR[i][j], "Rayleigh axial correlation curve");
			if (m_DebugDisplay) displayGraph(kAxial, "Axial coherent scattering amount");

			// Compute lateral correlation curve
			auto rhoLateral{ extractCorrelationCurve(autocorrelation, m_PixelSizeY, lateral) };
			if (m_DebugDisplay) displayGraph(rhoLateral, "Lateral correlation curve");

			// Compute lateral coherent scattering
			auto kLateral{ computeCoherentScatteringAmount(rhoLateralR[i][j], rhoLateral) };
			if (m_DebugDisplay) displayGraph(rhoLateralR[i][j], "Rayleigh lateral correlation curve");
			if (m_DebugDisplay) displayGraph(kLateral, "Lateral coherent scattering amount");

			// Compute elevation coherent scattering (average of kAxial and kLateral)
			auto kElevation{ averageBetween2CurvesWithScaling(kAxial, m_Axial2ElevationResolutionRatio, kLateral, m_Lateral2ElevationResolutionRatio) };
			if (m_DebugDisplay) displayGraph(kElevation, "Elevation coherent scattering amount");

			// Compute block cross-correlation after second block translation
			auto translatedBlock{ getTranslatedBlock(m_InplaneDisp[i][j], i, j) };
			if (m_DebugDisplay) translatedBlock.display("First block after translation");
			if (m_DebugDisplay) m_BlockMatrix2[i][j].display("Second block");
			
			auto crossCorrelation{ normalized_correlation(m_BlockMatrix2[i][j], translatedBlock) };
			if (m_DebugDisplay) crossCorrelation.display("Cross correlation between first translated block and second block");
			double tx{}, ty{};
			auto correlation{ translation(crossCorrelation, tx, ty) };

			// Resample elevation coherent scattering
			auto kElevationResampled{ resampleCurve(kElevation, rhoElevationR[i][j]) };
			if (m_DebugDisplay) displayGraph(kElevation, "Resampled elevation coherent scattering amount");

			// Compute elevation decorrelation curve
			auto rhoElevation{ computeElevationDecorrelationCurve(kElevationResampled, rhoElevationR[i][j])};
			if (m_DebugDisplay) displayGraph(rhoElevationR[i][j], "Rayleigh elevation correlation curve");
			if (m_DebugDisplay) displayGraph(rhoElevation, "Elevation correlation curve");

			// Get elevation corresponding 
			m_OutplaneDisp[i][j] = abscissaLinearInterpolateAt(rhoElevation, correlation);
		}
	}
}

RigidTransform EchoPair::computeRigidTransform(EchoPair::MatrixCurve rhoAxialR, EchoPair::MatrixCurve rhoLateralR, EchoPair::MatrixCurve rhoElevationR,
	int minimalSampleSize, int minimalInliersNumber, int iterationNumber, double distanceThreshold)
{
	// Compute inplane displacement using maximal cross-correlation
	computeInplaneDisplacement();

	// Compute outplane component using signal decorrelation
	computeOutplaneDisplacement(rhoAxialR, rhoLateralR, rhoElevationR);

	// Compute the corresponding rigid transform
	Points points1, points2;
	for (int i{}; i < m_NumBlockX; ++i)
	{
		double x1{ static_cast<double>(i * m_BlockWidth + m_BlockWidth / 2) * m_PixelSizeX };
		for (int j{}; j < m_NumBlockY; ++j)
		{
			double y1{ static_cast<double>(j * m_BlockHeight + m_BlockHeight / 2) * m_PixelSizeY };
			double z1{};
			points1.push_back(Eigen::Vector3d(x1, y1, z1));
			double x2{ x1 + m_InplaneDisp[i][j][0] * m_PixelSizeX };
			double y2{ y1 + m_InplaneDisp[i][j][1] * m_PixelSizeY };
			double z2{ m_OutplaneDisp[i][j]};
			points2.push_back(Eigen::Vector3d(x2, y2, z2));
		}
	}

	//auto [RAB, tAB] { Procrustes(points1, points2) };
	auto [rigidTransform, inliers, error] { RANSAC(points1, points2, minimalSampleSize, iterationNumber, distanceThreshold, minimalInliersNumber) };

	m_Outliers.resize(m_NumBlockX, std::vector<bool>(m_NumBlockY, true));
	for (const auto& inlier : inliers)
	{
		auto i{ inlier / m_NumBlockY };
		auto j{ inlier % m_NumBlockY };
		m_Outliers[i][j] = false;
	}

	return rigidTransform;
}

ImageUInt8 EchoPair::getDisplay()
{
	ImageUInt8 display(*m_Img1);
	
	// Draw arrows for inplane displacements
	const unsigned char red[]{ 255, 0, 0 };
	const unsigned char green[]{ 0, 255, 0 };
	for (int i{}; i < m_NumBlockX; ++i)
	{
		for (int j{}; j < m_NumBlockY; ++j)
		{
			int x0{ i * m_BlockWidth + m_BlockWidth / 2 }, x1{ x0 + (int)m_InplaneDisp[i][j][0] };
			int y0{ j * m_BlockHeight + m_BlockHeight / 2 }, y1{ y0 + (int)m_InplaneDisp[i][j][1] };
			if (m_Outliers[i][j])
				display.draw_arrow(x0, y0, x1, y1, red);
			else
				display.draw_arrow(x0, y0, x1, y1, green);
		}
	}

	// Draw grid at block limits
	const unsigned char black[]{ 0, 255, 0 };
	display.draw_grid((float)m_BlockWidth, (float)m_BlockHeight, 0.0f, 0.0f, false, false, black, 0.5);

	// Display text for elevation component
	const unsigned char cyan[]{ 0, 255, 255 };
	for (int i{}; i < m_NumBlockX; ++i)
	{
		for (int j{}; j < m_NumBlockY; ++j)
		{
			int x0{ i * m_BlockWidth + m_BlockWidth / 4 };
			int y0{ j * m_BlockHeight + m_BlockHeight / 2 };
			std::stringstream stream;
			stream << std::fixed << std::setprecision(2) << m_OutplaneDisp[i][j];
			display.draw_text(x0, y0, stream.str().c_str(), cyan, 0, 1, 8);
		}
	}

	return display;
}

void EchoPair::displayGraph(const Curve& curve, const std::string& title, const std::string& legendX, const std::string& legendY, PlotType type)
{
	const float x0{ static_cast<float>(curve.front().first) };
	const float x1{ static_cast<float>(curve.back().first) };
	const unsigned int plotType{ static_cast<unsigned int>(type) };
	const unsigned int vertexType{ 1 };
	const unsigned int n{ static_cast<unsigned int>(curve.size()) };

	ImageDouble values(1, n, 1, 1, 0);
	for (int j{}; j < curve.size(); ++j)
		values(0, j) = curve[j].second;

	values.display_graph(title.c_str(), plotType, vertexType, legendX.c_str(), x0, x1, legendY.c_str());
}