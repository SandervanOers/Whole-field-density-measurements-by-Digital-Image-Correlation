#include "DIC.hpp"
/*--------------------------------------------------------------------------*/
static bool sort_by_C_value (const  Points_With_Value &lhs, const Points_With_Value &rhs)
{
    return lhs.C_value > rhs.C_value;
}
/*--------------------------------------------------------------------------*/
static void compute_Save_GridX_Y(const cv::Size &Size, const InputVariables &inputvariables);
/*--------------------------------------------------------------------------*/
void store_matrix(std::string path, std::string filename, cv::Mat Matrix_To_Be_Stored)
{
    std::ofstream myfile;
    myfile.open(path+"/"+filename+".csv");
    myfile << format(Matrix_To_Be_Stored,cv::Formatter::FMT_MATLAB);
    myfile.close();
}
/*--------------------------------------------------------------------------*/
extern void DIC(const cv::Mat &img, const cv::Mat &img1, const InputVariables &inputvariables)
{
  cv::Mat M_valid_points = get_valid_points(img, inputvariables.SubsetLength, inputvariables.offset);
  /*--------------------------------------------------------------------------*/
  // Calculate interpolation coefficients for g
  cv::Mat prova_img1= img1.clone();
  float *fptr_img1 = prova_img1.ptr<float>(0);
  SamplesToCoefficients(fptr_img1, img1.cols, img1.rows, inputvariables.SplineDegree);
	std::cout << "\033[1;32mB-splines coefficients calculated\033[0m\n" << std::endl;
	/*--------------------------------------------------------------------------*/

	std::vector<cv::Mat> IG = calculateInitialGuess_Iteration(img, img1, fptr_img1, inputvariables);
	cv::Mat DispX = IG[0].clone();
	cv::Mat DispY = IG[1].clone();
	cv::Mat Ux =  	IG[2].clone();
	cv::Mat Vx =   	IG[3].clone();
	cv::Mat Uy =    IG[4].clone();
	cv::Mat Vy =    IG[5].clone();
	cv::Mat Uxy =   IG[6].clone();
	cv::Mat Vxy =   IG[7].clone();
	cv::Mat Uxx =   IG[8].clone();
	cv::Mat Vxx =   IG[9].clone();
	cv::Mat Uyy =   IG[10].clone();
	cv::Mat Vyy =   IG[11].clone();
	cv::Mat CorrelationCoefficient =  IG[12].clone();
	cv::Mat Computed_Points =  IG[13].clone();

	compute_Save_GridX_Y(DispX.size(), inputvariables);
	std::cout << "\033[1;32mInitial Points Computed\033[0m\n" << std::endl;
	/*--------------------------------------------------------------------------*/

    std::vector<Points_With_Value> Locations_Best_Correlation;
    cv::Mat nonZeroCoordinates;
    cv::findNonZero(CorrelationCoefficient>0, nonZeroCoordinates);
    for (unsigned int i = 0; i < nonZeroCoordinates.total(); i++ )
	{
		Locations_Best_Correlation.push_back(Points_With_Value(CorrelationCoefficient.at<double>(nonZeroCoordinates.at<Point>(i)), nonZeroCoordinates.at<Point>(i)));
    }

    // Sort first elements
    sort(Locations_Best_Correlation.begin(), Locations_Best_Correlation.end(), sort_by_C_value);
    for (auto i = Locations_Best_Correlation.begin(); i < Locations_Best_Correlation.end(); i++)
	{
		std::cout << (*i).Loc << ": " << (*i).C_value << std::endl;
	}
	std::cout << std::endl << "\033[1;32mList of Points with Values Sorted\033[0m\n" << std::endl;
	/*--------------------------------------------------------------------------*/
    const bool plotting = 0;
    /*int ct = 0;
    std::string name = "CC_";
    std::string type = ".png";
    std::string folderName = "Output/";
    if (plotting)
    {
        std::stringstream ss;
        ss<<name<<(ct)<<type;
        std::string filename = ss.str();
        ss.str("");
        cv::Mat Copy_CCF = CorrelationCoefficient.clone();
        Copy_CCF.convertTo(Copy_CCF, CV_8UC1, 255.0);
        imwrite(filename, Copy_CCF);
    }*/
	/*--------------------------------------------------------------------------*/

	cv::Ptr<cv::Formatter> fmt = cv::Formatter::get(cv::Formatter::FMT_DEFAULT);
	fmt->set64fPrecision(2);
	fmt->set32fPrecision(2);

	int outputcount = 1;
	while(static_cast<unsigned int>(cv::sum(Computed_Points).val[0]) < Computed_Points.total())
    {
        // Get best point in queue
        cv::Point matchLoc = Locations_Best_Correlation[0].Loc;
        // Use Initial Guess from (neighbouring) initial point
        std::vector<double> InitialCondition = {DispX.at<double>(matchLoc), DispY.at<double>(matchLoc), Ux.at<double>(matchLoc), Vx.at<double>(matchLoc), Uy.at<double>(matchLoc), Vy.at<double>(matchLoc), Uxy.at<double>(matchLoc), Vxy.at<double>(matchLoc), Uxx.at<double>(matchLoc), Vxx.at<double>(matchLoc), Uyy.at<double>(matchLoc), Vyy.at<double>(matchLoc)};
        std::vector<cv::Point> Neighbours = get_valid_Neighbours(M_valid_points, Computed_Points, matchLoc.x, matchLoc.y, inputvariables.SubsetLength, inputvariables.GridLength, inputvariables.offset);
        for (auto i = Neighbours.begin(); i < Neighbours.end(); i++)
        {
            if (inputvariables.propagationfunction==1 &&  inputvariables.GridLength <  inputvariables.SubsetLength/2)
            {
                //std::cout << "InitialCondition old: " << InitialCondition[0] << ", " << InitialCondition[1] << std::endl;
                InitialCondition[0] += Ux.at<double>(matchLoc)*((*i).x-matchLoc.x)*(double)inputvariables.GridLength + Uy.at<double>(matchLoc)*((*i).y-matchLoc.y)*(double)inputvariables.GridLength;
                InitialCondition[1] += Vx.at<double>(matchLoc)*((*i).x-matchLoc.x)*(double)inputvariables.GridLength + Vy.at<double>(matchLoc)*((*i).y-matchLoc.y)*(double)inputvariables.GridLength;
                //std::cout << "InitialCondition new: " << InitialCondition[0] << ", " << InitialCondition[1] << std::endl;
                ///InitialCondition[2] += Uxx.at<double>(matchLoc)*((*i).x-matchLoc.x)*(double)GridLength;// + Uxy.at<double>(matchLoc)*((*i).y-matchLoc.y)*(double)GridLength;
                //InitialCondition[3] += Vxx.at<double>(matchLoc)*((*i).x-matchLoc.x)*(double)GridLength;// + Vxy.at<double>(matchLoc)*((*i).y-matchLoc.y)*(double)GridLength;
                //InitialCondition[4] += Uyy.at<double>(matchLoc)*((*i).y-matchLoc.y)*(double)GridLength;// + Uxy.at<double>(matchLoc)*((*i).x-matchLoc.x)*(double)GridLength;
                ///InitialCondition[5] += Vyy.at<double>(matchLoc)*((*i).y-matchLoc.y)*(double)GridLength;// + Vxy.at<double>(matchLoc)*((*i).x-matchLoc.x)*(double)GridLength;
            }
			//auto tr3 = std::chrono::high_resolution_clock::now();
			//std::cout << "before" << std::endl;
            std::vector<double> point2 = iteration(img, fptr_img1, img1.rows, img1.cols, (*i).x, (*i).y, InitialCondition, inputvariables.SplineDegree, inputvariables.SubsetLength, inputvariables.GridLength, inputvariables.abs_tolerance_threshold, inputvariables.rel_tolerance_threshold, inputvariables.ShapeFunction, inputvariables.xDiff, inputvariables.yDiff);//inputvariables.xStart, inputvariables.yStart);
			//std::cout << "after" << std::endl;
			//auto tr4 = std::chrono::high_resolution_clock::now();
			//std::chrono::duration<double> elapsed_seconds = tr4-tr3;
			//T += elapsed_seconds.count();
            DispX.at<double>(*i) = point2[0];
            DispY.at<double>(*i) = point2[1];
            Ux.at<double>(*i) = point2[2];
            Vx.at<double>(*i) = point2[3];
            Uy.at<double>(*i) = point2[4];
            Vy.at<double>(*i) = point2[5];
            Uxy.at<double>(*i) = point2[6];
            Vxy.at<double>(*i) = point2[7];
            Uxx.at<double>(*i) = point2[8];
            Vxx.at<double>(*i) = point2[9];
            Uyy.at<double>(*i) = point2[10];
            Vyy.at<double>(*i) = point2[11];
			//std::cout << "11" << std::endl;
            CorrelationCoefficient.at<double>(*i) = point2.back();
			//std::cout << "12: " << point2.back() << std::endl;
            Computed_Points.at<uchar>(*i) = 1;
			//std::cout << "13" << std::endl;

			//std::cout << fmt->format(CorrelationCoefficient) << std::endl;
			//std::cout << *i << std::endl << std::endl;
			//std::cout << std::fixed << std::setprecision(2) << CorrelationCoefficient << std::endl;
            Locations_Best_Correlation.push_back(Points_With_Value(point2.back(), *i));

        }
		//std::cout << "Before Loc"<< std::endl;
        Locations_Best_Correlation.erase(Locations_Best_Correlation.begin());
        sort(Locations_Best_Correlation.begin(), Locations_Best_Correlation.end(), sort_by_C_value);

    //for (auto i = Locations_Best_Correlation.begin(); i < Locations_Best_Correlation.end(); i++)
	//{
	//	std::cout << (*i).Loc << ": " << (*i).C_value << std::endl;
	//}
	//std::cout << "After Loc"<< std::endl;

        if (Neighbours.empty())
        {
            //std::cout << "no valid or uncomputed neighbours" << std::endl;
        }
        else
        {
            if (plotting)
            {/*
                ct++;
                std::stringstream ss;
                ss<<name<<(ct)<<type;
                std::string filename = ss.str();
                ss.str("");
                cv::Mat Copy_CCF = CorrelationCoefficient.clone();
                Copy_CCF.convertTo(Copy_CCF, CV_8UC1, 255.0);
                imwrite(filename, Copy_CCF);*/
            }
            else
            {
                double computed = cv::sum(Computed_Points).val[0]/static_cast<double>(Computed_Points.total())*100.0;
                if (computed>outputcount)
                {
                    outputcount++;
                    outputcount++;
                    outputcount++;
                    std::cout << "Points Computed: " << std::setprecision(2) << computed << "%" << std::endl;
                }
            }
        }
    }
	std::cout << std::endl << "\033[1;32mAll Points Computed\033[0m\n" << std::endl;

	std::cout << "Reliability = " << std::setprecision(5) << cv::sum(CorrelationCoefficient).val[0]/static_cast<double>(CorrelationCoefficient.total()) << std::endl;
	/*--------------------------------------------------------------------------*/
  store_matrix(inputvariables.path,"U0", DispX);
  store_matrix(inputvariables.path,"V0", DispY);
  store_matrix(inputvariables.path,"Ux", Ux);
  store_matrix(inputvariables.path,"Vx", Vx);
  store_matrix(inputvariables.path,"Uy", Uy);
  store_matrix(inputvariables.path,"Vy", Vy);
  store_matrix(inputvariables.path,"Uxy", Uxy);
  store_matrix(inputvariables.path,"Vxy", Vxy);
  store_matrix(inputvariables.path,"Uxx", Uxx);
  store_matrix(inputvariables.path,"Uyy", Uyy);
  store_matrix(inputvariables.path,"Vxx", Vxx);
  store_matrix(inputvariables.path,"Vyy", Vyy);
  store_matrix(inputvariables.path,"CorrelationCoefficient", CorrelationCoefficient);
	std::cout << std::endl << "\033[1;32mSaving Matrices Completed\033[0m\n" << std::endl;
	/*--------------------------------------------------------------------------*/
}
/*--------------------------------------------------------------------------*/
static void compute_Save_GridX_Y(const cv::Size &Size, const InputVariables &inputvariables)
{
  cv::Mat GridX(Size, CV_64F, 0.0);
  cv::Mat GridY(Size, CV_64F, 0.0);
	for (unsigned int i = 0; i < static_cast<unsigned int>(GridX.cols); i++)
	{
		for (unsigned int j = 0; j < static_cast<unsigned int>(GridX.rows); j++)
		{
			GridX.at<double>(j,i) = static_cast<float>(inputvariables.xStart_ROI) + static_cast<float>(i*inputvariables.GridLength) - static_cast<float>(inputvariables.x_p0);
      GridY.at<double>(j,i) = static_cast<float>(inputvariables.yStart_ROI) + static_cast<float>(j*inputvariables.GridLength) - static_cast<float>(inputvariables.y_p0);
		}
	}
  store_matrix(inputvariables.path,"GridX", GridX);
	store_matrix(inputvariables.path,"GridY", GridY);
}
/*--------------------------------------------------------------------------*/
