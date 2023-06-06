#ifndef CSVImage
#define CSVImage
#include <opencv2/opencv.hpp>
#include <fstream>
/*--------------------------------------------------------------------------*/
extern cv::Mat load_matrix(std::string path, std::string filename, const int &skiplines);
/*--------------------------------------------------------------------------*/
#endif
