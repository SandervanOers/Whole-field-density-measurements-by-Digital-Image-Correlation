#ifndef InOut
#define InOut

//#include <stdio.h>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <fstream>
#include <thread>
#include <sys/stat.h>
#include "InputVariables.hpp"
#include "csvload.hpp"

#include <opencv2/opencv.hpp>
using namespace cv;
/*--------------------------------------------------------------------------*/
//extern int checkinput(unsigned int argc, char *argv[], std::string &pathname_string, unsigned int &SplineDegree, unsigned int &SubsetLength, unsigned int &GridLength, unsigned int &ShapeFunction, unsigned int &propagationfunction, unsigned int &ordering, unsigned int &xStart, unsigned int &xEnd, unsigned int &yStart, unsigned int &yEnd, unsigned int &offset, unsigned int &Number_Of_Threads, unsigned int &MaxPixelYVertical, double &abs_tolerance_threshold, double &rel_tolerance_threshold, double &minimum_corrcoeff_IG);
/*--------------------------------------------------------------------------*/
extern int readinput(unsigned int argc, char *argv[], InputVariables &inputvariables);
/*--------------------------------------------------------------------------*/
extern int checkImageDataFromFile(const cv::String &path, cv::Mat &img, cv::Mat &img1, const unsigned int &xStart, const unsigned int &xEnd, const unsigned int &yStart, const unsigned int &yEnd, const unsigned int &SubsetLength, const unsigned int &offset, unsigned int &xStart_ROI, unsigned int &yStart_ROI, unsigned int &horx_ROI, unsigned int &very_ROI, const unsigned int &ordering);
/*--------------------------------------------------------------------------*/
extern int readImageDataFromFile(cv::Mat &img, cv::Mat &img1, InputVariables &inputvariables);
/*--------------------------------------------------------------------------*/
#endif
