#ifndef DImage
#define DImage

#include <stdio.h>
#include <opencv2/opencv.hpp>

#include <math.h>
#include <numeric>
#include <algorithm>
#include <sstream>
#include <random>
#include <functional>
#include <sys/types.h>
#include <sys/stat.h>
#include <thread>
#include <future>
//#include <iostream>
//#include <fstream>

#include "InputOut.hpp"
#include "DIC.hpp"
#include "pixeltranslation.hpp"
#include "nonlineariteration.hpp"
#include "PointsWithValue.hpp"
#include "PositionDirection.hpp"
#include "InputVariables.hpp"
#include "ExperimentalSetupVariables.hpp"
#include "CalibrationValues.hpp"
#include "calculateN.hpp"
#include "csvload.hpp"

extern "C" {
#include "coeff.h"
#include "interpol.h"
}
using namespace cv;
/*--------------------------------------------------------------------------*/
static void store_matrix(std::string path, std::string filename, cv::Mat Matrix_To_Be_Stored);
/*--------------------------------------------------------------------------*/
static double calculateMean(const cv::Mat &Mat);
/*--------------------------------------------------------------------------*/
#endif
