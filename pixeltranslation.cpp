# include "pixeltranslation.hpp"
/*--------------------------------------------------------------------------*/
std::vector<double> calculatePixelTranslationRandom_SinglePoint(const cv::Mat &Reference, const cv::Mat &Deformed, const unsigned int &SubsetLength, const unsigned int &Indexi, const unsigned int &Indexj, const unsigned int &MaxPixelYVertical, const unsigned int &xStart, const unsigned int &yStart)
{
	cv::Mat temp = Reference(cv::Range(Indexj-SubsetLength/2,Indexj+SubsetLength/2+1), cv::Range(Indexi-SubsetLength/2, Indexi+SubsetLength/2+1));

	double DispX=0, DispY=0, CorrelationCoefficient=-1;

	// check wheter the reference image is (nearly) uniform.
	cv::Mat     mean;
	cv::Mat     stddev;
	cv::meanStdDev (temp, mean, stddev );
	double stddev_pxl = stddev.data[0];
	if (stddev_pxl>=1)
	{
		unsigned int dyl, dyr;
		if (Indexj+yStart>MaxPixelYVertical)
		{
			dyl = Indexj+yStart-MaxPixelYVertical;
		}
		else
		{
			dyl = 0;
		}
		if (Indexj+MaxPixelYVertical+1+yStart < static_cast<unsigned int>(Deformed.rows) )
		{
			dyr = Indexj+MaxPixelYVertical+1+yStart;
		}
		else
		{
			dyr = Deformed.rows;
		}
		cv::Mat Deformed2 = Deformed(cv::Range(dyl,dyr), cv::Range::all());

		cv::Mat result;
		double minVal; double maxVal;
		cv::Point minLoc; cv::Point maxLoc;
		cv::Point matchLoc;
		std::cout << "here" << std::endl;
		cv::matchTemplate(Deformed2, temp, result, cv::TM_CCOEFF_NORMED );  //CV_TM_CCOEFF_NORMED
		std::cout << "here 1s" << std::endl;
		// Localizing the best match with minMaxLoc
		cv::minMaxLoc( result, &minVal, &maxVal, &minLoc, &maxLoc, cv::Mat() );
		matchLoc = maxLoc;
		//DispY = (double)matchLoc.y + (SubsetLength/2) + dyl -(double)Indexj;
		//DispX = (double)matchLoc.x + SubsetLength/2 - (double)Indexi;
		DispY = (double)matchLoc.y + (SubsetLength/2) + dyl -(double)Indexj - (double)yStart;
		DispX = (double)matchLoc.x + SubsetLength/2 - (double)Indexi - (double)xStart;
		CorrelationCoefficient = maxVal;
	}
    std::vector<double> ReturnVector;
    ReturnVector.push_back(DispX);
    ReturnVector.push_back(DispY);
    ReturnVector.push_back(CorrelationCoefficient);
    return ReturnVector;
}
/*--------------------------------------------------------------------------*/
void calculateInitialGuess_Thread_Iteration(const unsigned int &Number_Of_Threads, const cv::Mat &Reference, const cv::Mat &Deformed, float *fptr_img1, cv::Mat &DispX, cv::Mat &DispY, cv::Mat &Ux, cv::Mat &Vx, cv::Mat &Uy, cv::Mat &Vy, cv::Mat &Uxy, cv::Mat &Vxy, cv::Mat &Uxx, cv::Mat &Vxx, cv::Mat &Uyy, cv::Mat &Vyy, cv::Mat &CorrelationCoefficient, cv::Mat &Computed_Points, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &offset, const unsigned int &xl, const unsigned int &xr, const unsigned int &yl, const unsigned int &yr, const unsigned int &MaxPixelYVertical, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold, const unsigned int &ShapeFunction, const double &minimum_corrcoeff_IG, const unsigned int &xStart, const unsigned int &yStart)
{
	thread_local std::uniform_int_distribution<unsigned> ux(xl, xr);
  thread_local std::uniform_int_distribution<unsigned> uy(yl, yr);
	std::random_device rdx{};
    std::random_device rdy{};
    std::default_random_engine ex{rdx()};
    std::default_random_engine ey{rdy()};
    auto dicex = std::bind(ux, ex);
    auto dicey = std::bind(uy, ey);

	std::map<double, Points_With_Location_And_Data> List_Of_Initial_Guesses;
    for (unsigned int k = 0; k < ceil(100/Number_Of_Threads); k++)
    {
    // Pick X_positions and Y_positions randomly on Grid (from GridLength)
    // Get random i in range (1,horx_ROI/GridLength) => unsigned int Indexi = offset+SubsetLength/2 + i * GridLength; => X_positions(k) = Indexi;
    // Get random j in range (1,very_ROI/GridLength) => unsigned int Indexj = offset+SubsetLength/2 + j * GridLength; => Y_posiitons(k) = Indexj;
    // Do calculatePixelTranslation for these randomly chosen points
		unsigned int i = dicex();
        unsigned int j = dicey();
        unsigned int Indexi = offset+SubsetLength/2 + i * GridLength;
        unsigned int Indexj = offset+SubsetLength/2 + j * GridLength;
        std::vector <double> Displ;
        Displ = calculatePixelTranslationRandom_SinglePoint(Reference, Deformed, SubsetLength, Indexi, Indexj, MaxPixelYVertical, xStart, yStart);

		List_Of_Initial_Guesses.insert(std::pair<double, Points_With_Location_And_Data>(Displ[2],Points_With_Location_And_Data(Displ[0], Displ[1], i, j) ));
	}

	bool better_solution = 0;
	while (!List_Of_Initial_Guesses.empty() && !better_solution)
	{
		double CC = (--List_Of_Initial_Guesses.end())->first;
		double DX = (--List_Of_Initial_Guesses.end())->second.DX;
		double DY = (--List_Of_Initial_Guesses.end())->second.DY;
		unsigned int I = (--List_Of_Initial_Guesses.end())->second.I;
		unsigned int J = (--List_Of_Initial_Guesses.end())->second.J;
		// Remove Last Entry from map
		std::map<double, Points_With_Location_And_Data>:: iterator it = List_Of_Initial_Guesses.end();
		it--;
		List_Of_Initial_Guesses.erase(it);
		// Use the Initial Guess to Calculate the Iterated Solution
		std::vector<double> InitialCondition = {DX, DY, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

		std::cout << "DX = " << DX << ", DY = " << DY << ", CC = " << CC << std::endl;
		// TO DO: Need to rewrite nonlineariteration.cpp to use img1.cols and img1.rows instead of img.cols and img.rows since interpol.c needs width and height of image.
		std::vector<double> point1 = iteration(Reference, fptr_img1, Deformed.rows, Deformed.cols, I, J, InitialCondition, SplineDegree, SubsetLength, GridLength, abs_tolerance_threshold, rel_tolerance_threshold, ShapeFunction, xStart, yStart);
		std::cout << "DX = " << point1[0] << ", DY = " << point1[1] << ", CC = " << point1[12] << std::endl << std::endl;
		if (point1.back()>CC && point1.back()>minimum_corrcoeff_IG)
		{
			better_solution = 1;
			DispX.at<double>(J,I) = point1[0];
			DispY.at<double>(J,I)= point1[1];
			Ux.at<double>(J,I) = point1[2];
			Vx.at<double>(J,I) = point1[3];
			Uy.at<double>(J,I) = point1[4];
			Vy.at<double>(J,I) = point1[5];
			Uxy.at<double>(J,I) = point1[6];
			Vxy.at<double>(J,I) = point1[7];
			Uxx.at<double>(J,I) = point1[8];
			Vxx.at<double>(J,I)= point1[9];
			Uyy.at<double>(J,I) = point1[10];
			Vyy.at<double>(J,I) = point1[11];
			CorrelationCoefficient.at<double>(J,I) = point1[12];
			Computed_Points.at<uchar>(J,I) = 1;
		}
	}
}
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculateInitialGuess_Iteration2(const cv::Mat &Reference, const cv::Mat &Deformed, float *fptr_img1, const unsigned int &SplineDegree, const unsigned int &SubsetLength, const unsigned int &GridLength, const unsigned int &horx_ROI, const unsigned int &very_ROI, const unsigned int &offset, const unsigned int &Number_Of_Threads, const unsigned int &MaxPixelYVertical, const double &abs_tolerance_threshold, const double &rel_tolerance_threshold, const unsigned int &ShapeFunction, const double &minimum_corrcoeff)
{
	//std::cout << "Computing Initial Guess" << std::endl;
	cv::Mat DispX(very_ROI/GridLength+1, horx_ROI/GridLength+1, CV_64F, 0.0);
    cv::Mat DispY(very_ROI/GridLength+1, horx_ROI/GridLength+1, CV_64F, 0.0);
    cv::Mat CorrelationCoefficient(very_ROI/GridLength+1, horx_ROI/GridLength+1, CV_64F, 0.0);
    cv::Mat Computed_Points(DispX.size(), CV_8UC1, 0.0);
    cv::Mat Ux(DispX.size(), CV_64FC1, 0.0);
    cv::Mat Vx(DispX.size(), CV_64FC1, 0.0);
    cv::Mat Uy(DispX.size(), CV_64FC1, 0.0);
    cv::Mat Vy(DispX.size(), CV_64FC1, 0.0);
    cv::Mat Uxy(DispX.size(), CV_64FC1, 0.0);
    cv::Mat Vxy(DispX.size(), CV_64FC1, 0.0);
    cv::Mat Uxx(DispX.size(), CV_64FC1, 0.0);
    cv::Mat Vxx(DispX.size(), CV_64FC1, 0.0);
    cv::Mat Uyy(DispX.size(), CV_64FC1, 0.0);
    cv::Mat Vyy(DispX.size(), CV_64FC1, 0.0);

	//std::cout << "Reference.rows = " << Reference.rows << ", Reference.cols = " << Reference.cols << std::endl;
	//std::cout << "Deformed.rows = " << Deformed.rows << ", Deformed.cols = " << Deformed.cols << std::endl;
	//std::cout << "Horizontal Difference = "<< Deformed.rows - Reference.rows << std::endl;
	//std::cout << "Vertical Difference = "<< Deformed.cols - Reference.cols << std::endl;

	unsigned int HorizontalDifference_DeformedReference = Deformed.rows - Reference.rows;
	unsigned int VerticalDifference_DeformedReference = Deformed.cols - Reference.cols;

	std::vector<std::thread> threads;
	for (unsigned int l = 0; l < Number_Of_Threads; l++)
	{
		unsigned int xl = static_cast<unsigned int>(1+round(0.1*horx_ROI/GridLength));
		unsigned int xr = static_cast<unsigned int>(horx_ROI/GridLength-round(0.1*horx_ROI/GridLength));
		unsigned int yl = static_cast<unsigned int>(1+round(0.1*very_ROI/GridLength) + (double)l/Number_Of_Threads * (very_ROI/GridLength-round(0.1*very_ROI/GridLength)-(1+round(0.1*very_ROI/GridLength))));
		unsigned int yr = static_cast<unsigned int>(1+round(0.1*very_ROI/GridLength) + (l+1.0)/Number_Of_Threads * (very_ROI/GridLength-round(0.1*very_ROI/GridLength)-(1+round(0.1*very_ROI/GridLength)))-1);
		// Possible Error when Number_Of_Threads is large and GridLength is large: yl == yr or even yl>yr.
		threads.push_back(std::thread(calculateInitialGuess_Thread_Iteration, Number_Of_Threads, Reference, Deformed, fptr_img1, std::ref(DispX), std::ref(DispY), std::ref(Ux), std::ref(Vx), std::ref(Uy), std::ref(Vy), std::ref(Uxy), std::ref(Vxy), std::ref(Uxx), std::ref(Vxx), std::ref(Uyy), std::ref(Vyy), std::ref(CorrelationCoefficient), std::ref(Computed_Points), SplineDegree, SubsetLength, GridLength, offset, xl, xr, yl, yr, MaxPixelYVertical, abs_tolerance_threshold, rel_tolerance_threshold, ShapeFunction, minimum_corrcoeff, HorizontalDifference_DeformedReference, VerticalDifference_DeformedReference));
	}

	for (auto& th : threads)
		th.join();

	//std::cout << "Computation Initial Guess Completed " << std::endl;
	std::vector<cv::Mat> ReturnVector;
    ReturnVector.push_back(DispX);
    ReturnVector.push_back(DispY);
    ReturnVector.push_back(Ux);
    ReturnVector.push_back(Vx);
    ReturnVector.push_back(Uy);
    ReturnVector.push_back(Vy);
    ReturnVector.push_back(Uxy);
    ReturnVector.push_back(Vxy);
    ReturnVector.push_back(Uxx);
    ReturnVector.push_back(Vxx);
    ReturnVector.push_back(Uyy);
    ReturnVector.push_back(Vyy);
    ReturnVector.push_back(CorrelationCoefficient);
    ReturnVector.push_back(Computed_Points);
    return ReturnVector;
}
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculateInitialGuess_Iteration(const cv::Mat &Reference, const cv::Mat &Deformed, float *fptr_img1, const InputVariables &inputvariables)
{
	cv::Mat DispX(inputvariables.very_ROI/inputvariables.GridLength+1, inputvariables.horx_ROI/inputvariables.GridLength+1, CV_64F, 0.0);
  cv::Mat DispY(inputvariables.very_ROI/inputvariables.GridLength+1, inputvariables.horx_ROI/inputvariables.GridLength+1, CV_64F, 0.0);
  cv::Mat CorrelationCoefficient(inputvariables.very_ROI/inputvariables.GridLength+1, inputvariables.horx_ROI/inputvariables.GridLength+1, CV_64F, 0.0);
  cv::Mat Computed_Points(DispX.size(), CV_8UC1, 0.0);
  cv::Mat Ux(DispX.size(), CV_64FC1, 0.0);
  cv::Mat Vx(DispX.size(), CV_64FC1, 0.0);
  cv::Mat Uy(DispX.size(), CV_64FC1, 0.0);
  cv::Mat Vy(DispX.size(), CV_64FC1, 0.0);
  cv::Mat Uxy(DispX.size(), CV_64FC1, 0.0);
  cv::Mat Vxy(DispX.size(), CV_64FC1, 0.0);
  cv::Mat Uxx(DispX.size(), CV_64FC1, 0.0);
  cv::Mat Vxx(DispX.size(), CV_64FC1, 0.0);
  cv::Mat Uyy(DispX.size(), CV_64FC1, 0.0);
  cv::Mat Vyy(DispX.size(), CV_64FC1, 0.0);

	std::cout << "xStart = " << inputvariables.xStart << std::endl;
	std::cout << "yStart = " << inputvariables.yStart << std::endl;

	std::vector<std::thread> threads;
	for (unsigned int l = 0; l < inputvariables.Number_Of_Threads; l++)
	{
		unsigned int xl = static_cast<unsigned int>(1+round(0.1*inputvariables.horx_ROI/inputvariables.GridLength));
		unsigned int xr = static_cast<unsigned int>(inputvariables.horx_ROI/inputvariables.GridLength-round(0.1*inputvariables.horx_ROI/inputvariables.GridLength));
		unsigned int yl = static_cast<unsigned int>(1+round(0.1*inputvariables.very_ROI/inputvariables.GridLength) + (double)l/inputvariables.Number_Of_Threads * (inputvariables.very_ROI/inputvariables.GridLength-round(0.1*inputvariables.very_ROI/inputvariables.GridLength)-(1+round(0.1*inputvariables.very_ROI/inputvariables.GridLength))));
		unsigned int yr = static_cast<unsigned int>(1+round(0.1*inputvariables.very_ROI/inputvariables.GridLength) + (l+1.0)/inputvariables.Number_Of_Threads * (inputvariables.very_ROI/inputvariables.GridLength-round(0.1*inputvariables.very_ROI/inputvariables.GridLength)-(1+round(0.1*inputvariables.very_ROI/inputvariables.GridLength)))-1);
		// Possible Error when Number_Of_Threads is large and GridLength is large: yl == yr or even yl>yr.
		threads.push_back(std::thread(calculateInitialGuess_Thread_Iteration, inputvariables.Number_Of_Threads, Reference, Deformed, fptr_img1, std::ref(DispX), std::ref(DispY), std::ref(Ux), std::ref(Vx), std::ref(Uy), std::ref(Vy), std::ref(Uxy), std::ref(Vxy), std::ref(Uxx), std::ref(Vxx), std::ref(Uyy), std::ref(Vyy), std::ref(CorrelationCoefficient), std::ref(Computed_Points), inputvariables.SplineDegree, inputvariables.SubsetLength, inputvariables.GridLength, inputvariables.offset, xl, xr, yl, yr, inputvariables.MaxPixelYVertical, inputvariables.abs_tolerance_threshold, inputvariables.rel_tolerance_threshold, inputvariables.ShapeFunction, inputvariables.minimum_corrcoeff_IG, inputvariables.xDiff, inputvariables.yDiff ));//inputvariables.xStart, inputvariables.yStart));

	}

	for (auto& th : threads)
		th.join();

	//std::cout << "Computation Initial Guess Completed " << std::endl;
	std::vector<cv::Mat> ReturnVector;
    ReturnVector.push_back(DispX);
    ReturnVector.push_back(DispY);
    ReturnVector.push_back(Ux);
    ReturnVector.push_back(Vx);
    ReturnVector.push_back(Uy);
    ReturnVector.push_back(Vy);
    ReturnVector.push_back(Uxy);
    ReturnVector.push_back(Vxy);
    ReturnVector.push_back(Uxx);
    ReturnVector.push_back(Vxx);
    ReturnVector.push_back(Uyy);
    ReturnVector.push_back(Vyy);
    ReturnVector.push_back(CorrelationCoefficient);
    ReturnVector.push_back(Computed_Points);
    return ReturnVector;
}
