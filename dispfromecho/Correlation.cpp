#include <iostream>

#include "Correlation.h"

typedef struct
{
	int x, y;
}Point;

using namespace cimg_library;

void covariance(const ImageDouble& input1, const ImageDouble& input2, ImageDouble& output)
{
	// Zero padding if images are not of the same size
	unsigned int max_dimx = (input1.width() > input2.width()) ? input1.width() : input2.width();
	unsigned int max_dimy = (input1.height() > input2.height()) ? input1.height() : input2.height();
	CImgList<double> list_input1(2, max_dimx, max_dimy);
	if (input1.width() == max_dimx && input1.height() == max_dimy)
		list_input1[0] = input1;
	else
	{
		list_input1[0].fill(0);
		cimg_forXY(input1, x, y) { list_input1[0](x, y) = input1(x, y); }
	}
	list_input1[1].fill(0);
	CImgList<double> list_input2(2, max_dimx, max_dimy);
	if (input2.width() == max_dimx && input2.height() == max_dimy)
		list_input2[0] = input2;
	else
	{
		list_input2[0].fill(0);
		cimg_forXY(input2, x, y) { list_input2[0](x, y) = input2(x, y); }
	}
	list_input2[1].fill(0);

	CImgList<double> spectrum_input1 = list_input1.get_FFT();
	CImgList<double> spectrum_input2 = list_input2.get_FFT();
	CImgList<double> spectrum_ri(2, max_dimx, max_dimy);
	for (int i = 0; i < max_dimx; i++)
		for (int j = 0; j < max_dimy; j++)
		{
			double real = spectrum_input1[0](i, j) * spectrum_input2[0](i, j) + spectrum_input1[1](i, j) * spectrum_input2[1](i, j);
			double imaginary = spectrum_input1[1](i, j) * spectrum_input2[0](i, j) - spectrum_input2[1](i, j) * spectrum_input1[0](i, j);
			spectrum_ri[0](i, j) = real;
			spectrum_ri[1](i, j) = imaginary;
		}
	spectrum_ri[0](0, 0) = spectrum_input1[0](0, 0) * spectrum_input2[0](0, 0);
	spectrum_ri[1](0, 0) = spectrum_input1[1](0, 0) * spectrum_input2[1](0, 0);

	CImgList<double> list_output = spectrum_ri.get_FFT(true);

	int width = list_output[0].width(), dimx_2 = width / 2;
	int height = list_output[1].height(), dimy_2 = height / 2;
	output = CImg<double>(width, height);
	cimg_forXY(output, i, j)
	{
		int x = (i + dimx_2) % width;
		int y = (j + dimy_2) % height;
		output(x, y) = list_output[0](i, j);
	}
}

void correlation(const ImageDouble& input1, const ImageDouble& input2, ImageDouble& output)
{
	CImg<double> input_m1(input1);
	input_m1 -= input1.mean();
	CImg<double> input_m2(input2);
	input_m2 -= input2.mean();

	return covariance(input_m1, input_m2, output);
}

ImageDouble normalized_correlation(const ImageDouble& input1, const ImageDouble& input2)
{
	ImageDouble covBig, covBig2;
	double sigmaSmall;
	if (input1.width() > input2.width() && input1.height() > input2.height())
	{
		CImg<double> U(input2.width(), input2.height());
		U.fill(1.0);
		//U.display("U");
		covariance(input1.get_sqr(), U, covBig2);
		covBig2 *= input2.width() * input2.height();
		covariance(input1, U, covBig);
		sigmaSmall = std::sqrt(input2.variance());
	}
	else
	{
		CImg<double> U(input1.width(), input1.height());
		U.fill(1.0);
		covariance(input2.get_sqr(), U, covBig2);
		covBig2 *= input1.width() * input1.height();
		covariance(input2, U, covBig);
		sigmaSmall = std::sqrt(input1.variance());
	}

	CImg<double> den = covBig2 - covBig.get_sqr();
	den.sqrt();
	den *= sigmaSmall;

	CImg<double> num;
	correlation(input1, input2, num);

	return num.get_div(den);
}

double translation(ImageDouble& input, double& tx, double& ty)
{
	double max = 0;
	tx = ty = 0;
	for (int i = 0; i<input.width(); i++)
		for (int j = 0; j<input.height(); j++)
		{
			double valtemp = input(i, j);
			if (valtemp > max)
			{
				max = valtemp;
				tx = i;
				ty = j;
			}
		}
	tx -= input.width() / 2;
	ty -= input.height() / 2;
	return max;
}
