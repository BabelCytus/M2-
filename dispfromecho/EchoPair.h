#ifndef __EchoPair_H
#define __EchoPair_H

#include <vector>
#include <array>

#include <Eigen/Dense>

#include "Correlation.h"
#include "RigidTransform.h"
#include "CImg.h"


using ImageUInt8 = cimg_library::CImg<unsigned char>;
using MatrixBlockDouble = std::vector<std::vector<ImageDouble>>;
using MatrixDouble = std::vector<std::vector<double>>;
using Vec2 = std::array<double, 2>;
using MatrixVec2 = std::vector<std::vector<Vec2>>;

class EchoPair
{

public:
	enum CorrelationType { axial, lateral, elevation };
	enum PlotType { poinst, segments, splines, bars };
	using Point = std::pair<double, double>;
	using Curve = std::vector<Point>;
	using MatrixOutliers = std::vector<std::vector<bool>>;
	using MatrixCurve = std::vector<std::vector<Curve>>;

	EchoPair(ImageUInt8 *img1, ImageUInt8 *img2, int numX, int numY, double fovX, double fovY, double ratioX, double ratioY);
	~EchoPair();
	ImageUInt8* getImg1() { return m_Img1; }
	ImageUInt8* getImg2() { return m_Img2; }
	void deleteImg1() { if (m_Img1) { delete m_Img1; m_Img1 = nullptr; } }
	void deleteImg2() { if (m_Img2) { delete m_Img2; m_Img2 = nullptr; } }
	void setNumBlocks(int numX, int numY);
	void setDebugDisplay(bool display) { m_DebugDisplay = display; }
	int getNumBlockX() { return m_NumBlockX; }
	int getNumBlockY() { return m_NumBlockY; }
	double getPixelSizeX() { return m_PixelSizeX; }
	double getPixelSizeY() { return m_PixelSizeY; }
	const ImageDouble& getImg1Block(int i, int j) const { return m_BlockMatrix1[i][j]; }
	static Curve extractCorrelationCurve(const ImageDouble& correlationImage, double pixelSize, CorrelationType type);
	ImageDouble getBlockAutocorrelation(int i, int j);
	RigidTransform computeRigidTransform(MatrixCurve rhoAxialR, MatrixCurve rhoLateralR, MatrixCurve rhoElevationR, int minimalSampleSize, int minimalInliersNumber, int iterationNumber, double distanceThreshold);
	ImageUInt8 getDisplay();

private:
	ImageUInt8* m_Img1{ nullptr }, * m_Img2{ nullptr };
	MatrixBlockDouble m_BlockMatrix1, m_BlockMatrix2;
	MatrixVec2 m_InplaneDisp;
	MatrixDouble m_OutplaneDisp;
	Eigen::Matrix3d m_Rotation;
	Eigen::Vector3d m_Translation;
	MatrixOutliers m_Outliers;
	int m_NumBlockX{}, m_NumBlockY{};
	int m_BlockWidth{}, m_BlockHeight{};
	double m_PixelSizeX{}, m_PixelSizeY{};
	double m_Axial2ElevationResolutionRatio{}, m_Lateral2ElevationResolutionRatio{};
	bool m_DebugDisplay{ false };
	MatrixCurve m_rhoAxial, m_rhoLateral, m_rhoElevation;
	MatrixCurve m_kAxial, m_kLateral, m_kElevation;
	void setPixelSize(double fovX, double fovY);
	void blockDecomposition();
	void deleteBlocks();
	ImageDouble getTranslatedBlock(const Vec2& tranlation, int i, int j);
	void computeInplaneDisplacement();
	Curve computeCoherentScatteringAmount(const Curve& rhoR, const Curve& rho);
	double ordinateLinearInterpolateAt(const Curve& curve, double x);
	double abscissaLinearInterpolateAt(const Curve& curve, double y);
	Curve averageBetween2CurvesWithScaling(const Curve& curve1, double ratio1, const Curve& curve2, double ratio2);
	Curve resampleCurve(const Curve& curve, const Curve& curveRef);
	Curve computeElevationDecorrelationCurve(const Curve& kElevation, const Curve& rhoElevation);
	void computeOutplaneDisplacement(const MatrixCurve& rhoAxialR, const MatrixCurve& rhoLateralR, const MatrixCurve& rhoElevationR);
	void displayGraph(const Curve& curve, const std::string& title, const std::string& legendX = "X axis", const std::string& legendY = "Y axis", PlotType type = segments);
};

#endif // !__EchoPair_H