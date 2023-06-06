
#ifndef ITERATION
#define ITERATION

#include <opencv2/opencv.hpp>
#include <math.h>
#include <thread>
#include <numeric>
#include "PointsWithValue.hpp"
extern "C" {
#include "interpol.h"
}
/*--------------------------------------------------------------------------*/
extern std::vector<double> iteration(const cv::Mat &img, float *fptr_img1, const unsigned int &img1_rows, const unsigned int &img1_cols, const unsigned int &i, const unsigned int &j, const std::vector<double> &P0, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold, const unsigned int &ShapeFunction, const unsigned int &xStart, const unsigned int &yStart);
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Point> get_valid_Neighbours(const cv::Mat &M_valid_points, const cv::Mat &Computed_Points, const unsigned int &x, const unsigned int &y, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned &offset);
/*--------------------------------------------------------------------------*/
extern cv::Mat get_valid_points(const cv::Mat &img, const unsigned int &SubsetLength, const unsigned int &offset);
/*--------------------------------------------------------------------------*/
//static double getDerivativeValue(float *fptr_img1, const unsigned int &cols, const unsigned int &rows, const double &x, const double &y, const unsigned int &SplineDegree, const unsigned int &direction);
/*--------------------------------------------------------------------------*/
#endif
