#include "Calibration.h"

int main(void)
{
	const std::string srcPath{ "K:/Rob-echO/dispfromecho" };
	const bool calibration{ false };

	const int numBlockX{ 8 }, numBlockY{ 8 };	// Number of blocks in horizontal and vertical directions
	const double fovX{ 25.6 }, fovY{ 25.6 };	// Field of view in horizontal and vertical directions (in mm)
	
	// Calibration with phantom images supposed to be be pure fully developed speckle (Rayleigh distribution)
	if (calibration)
	{
		const std::string calibrationPath{ srcPath + "/Images/Acq3/Fantome/Elevation" };
		Calibration calibration(calibrationPath, numBlockX, numBlockY, fovX, fovY, "Vascular", 0.05);
		calibration.computeCalibrationDecorrelationCurves();
		//calibration.computeResolutionRatios();
		calibration.exportCalibrationDecorrelationCurves();
	}
	else
	{
		//const std::string imagePath{ srcPath + "/Images/Acq3/Fantome/Lateral" };
		const std::string imagePath{ srcPath + "/Images/Acq3/Fantome/Rotation" };
		//const std::string imagePath{ srcPath + "/Images/Acq3/Cuisse/Elevation" };
		//const std::string imagePath{ srcPath + "/Images/Acq3/Bras/Rotation" };

		// Ratio between autocorrelation widths in elevation direction vs axial and lateral directions
		// Could be measured automatically in the calibration step
		const double axial2ElevationRatio{ 4 }, lateral2ElevationRatio{ axial2ElevationRatio * 0.1 / 0.17 };
		double RANSACMinimalSampleSize{ 0.2 };
		double RANSACMinimalInliersNumber{ 0.4 }; 
		int RANSACIterationNumber{ 1000 };
		double RANSACDistanceThreshold{ 0.1 };
		const int timeStep{ 1 };
		EchoSequence seq(imagePath, numBlockX, numBlockY, fovX, fovY, axial2ElevationRatio, lateral2ElevationRatio,
			RANSACMinimalSampleSize,RANSACMinimalInliersNumber, RANSACIterationNumber, RANSACDistanceThreshold, timeStep);
		seq.setDebugDisplay(false);

		// New block decomposition
		//const int newNumBlockX{ 16 };
		//const int newNumBlockY{ 16 };
		//seq.setNumBlocks(newNumBlockX, newNumBlockX);

		// Import calibration curves (choose the files generated with the right number of blocks)
		seq.importRayleighDecorrelationCurves("Vascular_axial_8x8.csv", "Vascular_lateral_8x8.csv", "Vascular_elevation_8x8.csv");
		//seq.importRayleighDecorrelationCurves("Vascular_axial_16x16.csv", "Vascular_lateral_16x16.csv", "Vascular_elevation_16x16.csv");

		// Estimate all displacements between image pairs
		seq.computeDisplacements();

		// Display displacements between image pairs
		seq.displayDisplacements();
	}

	return 0;
}