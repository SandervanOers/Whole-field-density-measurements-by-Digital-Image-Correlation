#ifndef PIXEL_TRANSLATION
#define PIXEL_TRANSLATION

#include <opencv2/opencv.hpp>
#include <random>
#include <thread>
#include <chrono>
#include "nonlineariteration.hpp"
#include "InputVariables.hpp"
/*--------------------------------------------------------------------------*/
struct Points_With_Location_And_Data {
    double DX, DY;
    unsigned int I, J;
	Points_With_Location_And_Data(double k, double l, unsigned int i, unsigned int j) : DX(k), DY(l), I(i), J(j) {}
};
/*--------------------------------------------------------------------------*/
//std::vector<double> calculatePixelTranslationRandom_SinglePoint(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &Indexi, const unsigned int &Indexj, const InputVariables &inputvariables);
std::vector<double> calculatePixelTranslationRandom_SinglePoint(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &SubsetLength, const unsigned int &Indexi, const unsigned int &Indexj, const unsigned int &MaxPixelYVertical, const unsigned int &xStart, const unsigned int &yStart);
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculateInitialGuess_Iteration2(const cv::Mat &Reference, const cv::Mat &Deformed, float *fptr_img1, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &horx_ROI, const unsigned int &very_ROI, const unsigned int &offset, const unsigned int &Number_Of_Threads, const unsigned int &MaxPixelYVertical, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold, const unsigned int &ShapeFunction, const double &minimum_corrcoeff);
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculateInitialGuess_Iteration(const cv::Mat &Reference, const cv::Mat &Deformed, float *fptr_img1, const InputVariables &inputvariables);
/*--------------------------------------------------------------------------*/
//void calculateInitialGuess_Thread_Iteration(const cv::Mat &Reference, const cv::Mat &Deformed, float *fptr_img1, cv::Mat &DispX, cv::Mat &DispY, cv::Mat &Ux, cv::Mat &Vx, cv::Mat &Uy, cv::Mat &Vy, cv::Mat &Uxy, cv::Mat &Vxy, cv::Mat &Uxx, cv::Mat &Vxx, cv::Mat &Uyy, cv::Mat &Vyy, cv::Mat &CorrelationCoefficient, cv::Mat &Computed_Points, const InputVariables &inputvariables, const unsigned int &xl, const unsigned int &xr, const unsigned int &yl, const unsigned int &yr);
void calculateInitialGuess_Thread_Iteration(const unsigned int &Number_Of_Threads, const cv::Mat &Reference, const cv::Mat &Deformed, float *fptr_img1, cv::Mat &DispX, cv::Mat &DispY, cv::Mat &Ux, cv::Mat &Vx, cv::Mat &Uy, cv::Mat &Vy, cv::Mat &Uxy, cv::Mat &Vxy, cv::Mat &Uxx, cv::Mat &Vxx, cv::Mat &Uyy, cv::Mat &Vyy, cv::Mat &CorrelationCoefficient, cv::Mat &Computed_Points, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &offset, const unsigned int &xl, const unsigned int &xr, const unsigned int &yl, const unsigned int &yr, const unsigned int &MaxPixelYVertical, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold, const unsigned int &ShapeFunction, const double &minimum_corrcoeff, const unsigned int &xStart, const unsigned int &yStart);
/*--------------------------------------------------------------------------*/
#endif
