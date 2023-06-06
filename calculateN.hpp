#ifndef DetN
#define DetN

#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <numeric>
# include <vector>
# include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <future>

#include "SnellsLaw.hpp"
#include "PositionDirection.hpp"

//extern "C" {
//#include "coeff.h"
//#include "interpol.h"
//}
using namespace cv;
/*--------------------------------------------------------------------------*/
extern double calculateNorm(const double &a, const double &b, const double &c);
extern double computeLc(const double &c, const double &normPD, const double &L_m, const double &L_s, const double &L_g, const double &L_t);
/*--------------------------------------------------------------------------*/
extern void CalibrationFigures(const cv::Mat &GridX, const cv::Mat &GridY, const cv::Mat &Dx, const cv::Mat &Dy, const cv::Mat &CorrelationCoefficient, const double &focal_length, const std::vector<double> &Lengths, const double &Distance_From_Pixels_To_Meters, const double &n_0, const double &n_1, const double &n, const double &nref, const std::string &path, const double &corr_cut_off);
extern std::vector<double> Calibration(const cv::Mat &GridX, const cv::Mat &GridY, const cv::Mat &Dx, const cv::Mat &Dy, const cv::Mat &CorrelationCoefficient, const double &focal_length, const std::vector<double> &Lengths, const double &Distance_From_Pixels_To_Meters, const double &n_0, const double &n_1, const double &n, const double &nref, const std::string &path, const double &corr_cut_off, const int &directionsToInclude);
/*--------------------------------------------------------------------------*/
extern cv::Mat CalculateN(const cv::Mat &GridX, const cv::Mat &GridY, const double &meanGridX, const double &meanGridY, const cv::Mat &Dx, const cv::Mat &Dy, const double &focal_length, const std::vector<double> &Lengths, const double &Distance_From_Pixels_To_Meters, const std::vector<double> &PlaneDefinition, const double &n_0, const double &n_1, const double &nref, const int &directionsToInclude);
/*--------------------------------------------------------------------------*/
#endif
