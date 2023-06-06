#ifndef H_DIC
#define H_DIC

#include "pixeltranslation.hpp"
#include "nonlineariteration.hpp"
#include "PointsWithValue.hpp"
#include "InputVariables.hpp"
#include "InputOut.hpp"

extern "C" {
#include "coeff.h"
#include "interpol.h"
}
/*--------------------------------------------------------------------------*/
extern void DIC(const cv::Mat &img, const cv::Mat &img1, const InputVariables &inputvariables);
/*--------------------------------------------------------------------------*/
#endif
