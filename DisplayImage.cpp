#include "DisplayImage.hpp"

int main(int argc, char** argv )
{
	auto tr1 = std::chrono::high_resolution_clock::now();

	InputVariables inputvariables;
	int returnvalue = readinput(argc, argv, inputvariables);
	if (returnvalue < 0)
	{
		return -1;
	}
	std::cout << std::endl << "\033[1;32mInput Parsed\033[0m\n" << std::endl;
	/*--------------------------------------------------------------------------*/
	if (inputvariables.DICNeeded == 1)
	{
		// Read Images from Disk
		cv::Mat img, img1;
		returnvalue = readImageDataFromFile(img, img1, inputvariables);
		std::cout << std::endl << "\033[1;32mImages Loaded\033[0m\n" << std::endl;
		/*--------------------------------------------------------------------------*/
		DIC(img, img1, inputvariables);
		auto tr2= std::chrono::high_resolution_clock::now();
		std::cout << "DIC took " << std::chrono::duration_cast<std::chrono::milliseconds>(tr2-tr1).count()
		<< " milliseconds = " << std::chrono::duration_cast<std::chrono::seconds>(tr2-tr1).count() << " seconds = " << std::chrono::duration_cast<std::chrono::minutes>(tr2-tr1).count() << " minutes\n"<<std::endl;
	}
	/*--------------------------------------------------------------------------*/
	cv::Mat DX = load_matrix(inputvariables.path, "U0", 1);
	cv::Mat DY = load_matrix(inputvariables.path, "V0", 1);
	cv::Mat GridX = load_matrix(inputvariables.path, "GridX", 1);
	cv::Mat GridY = load_matrix(inputvariables.path, "GridY", 1);
	inputvariables.x_p0 = 0.0;
	inputvariables.y_p0 = 0.0;
	GridX = GridX - calculateMean(GridX);
	GridY = GridY - calculateMean(GridY);

	cv::Mat CC = load_matrix(inputvariables.path, "CorrelationCoefficient", 1);
	std::cout << std::endl << "\033[1;32mLoading Matrices Completed\033[0m\n" << std::endl;
	/*--------------------------------------------------------------------------*/
	auto tr3= std::chrono::high_resolution_clock::now();
	ExperimentalSetupVariables experimentalsetupvariables;
	// Camera
	// Camera Lyon (sider)
	experimentalsetupvariables.focal_length = 25.0e-3;
	experimentalsetupvariables.Distance_From_Pixels_To_Meters = 3.45e-6;
	// Lengths Lyon
	experimentalsetupvariables.L_c = 1e4;
	experimentalsetupvariables.L_g = 0.6/100;
	experimentalsetupvariables.L_t = 0.13;
	experimentalsetupvariables.L_s = 0.0;

	std::vector<double> Lengths{experimentalsetupvariables.L_c, experimentalsetupvariables.L_g, experimentalsetupvariables.L_t, experimentalsetupvariables.L_s};
	double corr_cut_off = inputvariables.minimum_corrcoeff_IG;
		CalibrationValues calibrationValues;
	if (inputvariables.CalibrationNeeded==1)
	{
		/*--------------------------------------------------------------------------*/
		// CalibrationFigures(GridX, GridY, DX, DY, CC, experimentalsetupvariables.focal_length, Lengths, experimentalsetupvariables.Distance_From_Pixels_To_Meters, experimentalsetupvariables.n_0, experimentalsetupvariables.n_1, experimentalsetupvariables.n, inputvariables.path, corr_cut_off);
		/*--------------------------------------------------------------------------*/
		// calibration of water with known refractive index nref to air (n_0)
		experimentalsetupvariables.n = 1.337;//inputvariables.nref;
		experimentalsetupvariables.n_ref = experimentalsetupvariables.n_0;
	  std::vector<double> CalibrationNumbers = Calibration(GridX, GridY, DX, DY, CC, experimentalsetupvariables.focal_length, Lengths, experimentalsetupvariables.Distance_From_Pixels_To_Meters, experimentalsetupvariables.n_0, experimentalsetupvariables.n_1, experimentalsetupvariables.n, experimentalsetupvariables.n_ref, inputvariables.path, corr_cut_off, inputvariables.directionsToInclude);
		calibrationValues.a = CalibrationNumbers[0];
		calibrationValues.b = CalibrationNumbers[1];
		calibrationValues.c = CalibrationNumbers[2];
		calibrationValues.L_m = CalibrationNumbers[3];
		calibrationValues.meanGridX = CalibrationNumbers[4];
		calibrationValues.meanGridY = CalibrationNumbers[5];
		std::cout << std::endl << "\033[1;32mCalibration Completed\033[0m\n" << std::endl;

		for (auto i = CalibrationNumbers.begin(); i != CalibrationNumbers.end(); ++i)
			std::cout << *i << ' ';

		std::cout << std::endl<< std::endl<< std::endl;
		auto tr4= std::chrono::high_resolution_clock::now();
		std::cout << "Calibration took " << std::chrono::duration_cast<std::chrono::milliseconds>(tr4-tr3).count()
		<< " milliseconds = " << std::chrono::duration_cast<std::chrono::seconds>(tr4-tr3).count() << " seconds = " << std::chrono::duration_cast<std::chrono::minutes>(tr4-tr3).count() << " minutes"<<std::endl;
	}
	else
	{
		// Load Numbers
		cv::Mat CalibrationNumbers = load_matrix(inputvariables.path, "CalibrationFigures", 1);
		std::cout << "CalibrationNumbers = " << CalibrationNumbers << std::endl;
		calibrationValues.a = CalibrationNumbers.at<double>(0);
		calibrationValues.b = CalibrationNumbers.at<double>(1);
		calibrationValues.c = CalibrationNumbers.at<double>(2);
		calibrationValues.L_m = CalibrationNumbers.at<double>(3);
		calibrationValues.meanGridX = CalibrationNumbers.at<double>(4);
		calibrationValues.meanGridY = CalibrationNumbers.at<double>(5);
		std::cout << std::endl << "\033[1;32mReading Calibration Completed\033[0m\n" << std::endl;
	}
	/*--------------------------------------------------------------------------*/
	if (inputvariables.CalculateRefractionIndex == 1)
	{
		std::cout << std::endl << "\033[1;32mComputing Refraction Index\033[0m\n" << std::endl;

		// calculation of unknown refractive index n to known known refractive index nref
		experimentalsetupvariables.n_ref = inputvariables.nref;

		double normPD = calculateNorm(calibrationValues.a, calibrationValues.b, calibrationValues.c);
		double d = -calibrationValues.c/normPD*calibrationValues.L_m;
		std::vector<double> PlaneDefinition{0, 0, 0, 0};
		PlaneDefinition[0] = calibrationValues.a/normPD;
		PlaneDefinition[1] = calibrationValues.b/normPD;
		PlaneDefinition[2] = calibrationValues.c/normPD;
		PlaneDefinition[3] = d;
		experimentalsetupvariables.L_c = computeLc(calibrationValues.c, normPD, calibrationValues.L_m, experimentalsetupvariables.L_s, experimentalsetupvariables.L_g, experimentalsetupvariables.L_t);
		Lengths[0] = experimentalsetupvariables.L_c;
		Lengths[2] = experimentalsetupvariables.L_t;

		auto tr5 = std::chrono::high_resolution_clock::now();
		std::vector<cv::Mat> nfieldandposition = CalculateN(GridX, GridY, calibrationValues.meanGridX, calibrationValues.meanGridY, DX, DY, experimentalsetupvariables.focal_length, Lengths, experimentalsetupvariables.Distance_From_Pixels_To_Meters, PlaneDefinition, experimentalsetupvariables.n_0, experimentalsetupvariables.n_1, experimentalsetupvariables.n_ref, inputvariables.directionsToInclude);
		store_matrix(inputvariables.path,"nfield", nfieldandposition[0]);
		store_matrix(inputvariables.path,"ksi", nfieldandposition[1]);
		store_matrix(inputvariables.path,"eta", nfieldandposition[2]);
		store_matrix(inputvariables.path,"zeta", nfieldandposition[3]);
		auto tr6= std::chrono::high_resolution_clock::now();
		std::cout << "n Calculation took " << std::chrono::duration_cast<std::chrono::milliseconds>(tr6-tr5).count()
		<< " milliseconds = " << std::chrono::duration_cast<std::chrono::seconds>(tr6-tr5).count() << " seconds = " << std::chrono::duration_cast<std::chrono::minutes>(tr6-tr5).count() << " minutes"<<std::endl;
		std::cout << std::endl << "\033[1;32mN Calculation Completed\033[0m\n" << std::endl;
	}
	auto trend = std::chrono::high_resolution_clock::now();
	std::cout << "Total took " << std::chrono::duration_cast<std::chrono::milliseconds>(trend-tr1).count()
	<< " milliseconds = " << std::chrono::duration_cast<std::chrono::seconds>(trend-tr1).count() << " seconds = " << std::chrono::duration_cast<std::chrono::minutes>(trend-tr1).count() << " minutes"<<std::endl;
    return 0;
}
/*--------------------------------------------------------------------------*/
void store_matrix(std::string path, std::string filename, cv::Mat Matrix_To_Be_Stored)
{
    std::ofstream myfile;
    myfile.open(path+"/"+filename+".csv");
    myfile << format(Matrix_To_Be_Stored,cv::Formatter::FMT_MATLAB);
    myfile.close();
}
/*--------------------------------------------------------------------------*/
static double calculateMean(const cv::Mat &Mat)
{
	Scalar meanMatrix = mean(Mat);
	return meanMatrix[0];
}
/*--------------------------------------------------------------------------*/
