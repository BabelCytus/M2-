#ifndef __Calibration_H
#define __Calibration_H

#include "EchoSequence.h"

class Calibration : public EchoSequence
{
public:
	Calibration(const std::string& dirName, int numX, int numY, double fovX, double fovY, const std::string& probeName, double elevationStep, double ratioX = 1.0, double ratioY = 1.0, int timeStep = 1);
	void computeCalibrationDecorrelationCurves();
	void exportCalibrationDecorrelationCurves();

private:
	std::string m_ProbeName{ "" };
	double m_ElevationStep{};
	EchoPair::MatrixCurve m_AxialDecorrelationCurves, m_LateralDecorrelationCurves, m_ElevationDecorrelationCurves;
};

#endif // !__Calibration_H
