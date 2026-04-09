#ifndef __Correlation_H
#define __Correlation_H

#include "CImg.h"

using ImageDouble = cimg_library::CImg<double>;

void covariance(const ImageDouble& input1, const ImageDouble& input2, ImageDouble& output);
void correlation(const ImageDouble& input1, const ImageDouble& input2, ImageDouble& output);
ImageDouble normalized_correlation(const ImageDouble& input1, const ImageDouble& input2);
double translation(ImageDouble& input, double& tx, double& ty);

#endif // !__Correlation_H
