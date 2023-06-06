#include "InputOut.hpp"

/*--------------------------------------------------------------------------*/
extern int readinput(unsigned int argc, char *argv[], InputVariables &inputvariables)
{
	unsigned int NumberOfInputArguments = 22;
	std::string commandline = "usage: DisplayImage.out <Image_Path> SplineDegree SubsetLength GridLength ShapeFunction PropagationFunction OrderingImages xStart xEnd yStart yEnd NumberOfThreads MaxPixelYVertical Tolerance MinCorrCoeffIG BlurSize directionsToInclude nref DICNeeded CalibrationNeeded CalculateRefractionIndex \n";
	// Checks number of Arguments
    if ( argc != NumberOfInputArguments)
    {
		printf("Too few arguments\n");
    printf("%s", commandline.c_str());
		return -1;
    }

	// Checks pathname
	const char *pathname;
	pathname = argv[1];
	struct stat info;

	if( stat( pathname, &info ) != 0 )
	{
		printf( "cannot access %s\n", pathname );
		return -1;
	}
	else if( info.st_mode & S_IFDIR )
	{
		//printf( "%s is a directory\n", pathname );
	}
	else
	{
		printf( "%s is no directory\n", pathname );
		return -1;
	}
	inputvariables.path = std::string(pathname);

	// Saves Command
	std::ofstream myfile;
	myfile.open (inputvariables.path+"/commandline.txt", std::ios_base::app);
	myfile << commandline;
	for (unsigned int i=0;i<NumberOfInputArguments;i++)
	{
		myfile << argv[i] << " ";
	}
	myfile << "\n";
	myfile.close();

	// Checks input
	unsigned int SplineDegree = atoi(argv[2]);
    inputvariables.SplineDegree = atoi(argv[2]);
    if (!(SplineDegree > 2 && SplineDegree < 9))
    {
        printf("Spline degree must be between 2 and 9 \n");
        return -1;
    }
	inputvariables.SplineDegree = SplineDegree;
    unsigned int SubsetLength = atoi(argv[3]);
    if (!(SubsetLength > 1))
    {
        printf("Subset Length must be larger than 1 \n");
        return -1;
    }
	inputvariables.SubsetLength = SubsetLength;
    unsigned int GridLength = atoi(argv[4]);
    if (!(GridLength > 0))
    {
        printf("GridLength must be larger than 0 \n");
        return -1;
    }
	inputvariables.GridLength = GridLength;
    unsigned int ShapeFunction = atoi(argv[5]);
    if (!(ShapeFunction==0 || ShapeFunction==1 || ShapeFunction==2 || ShapeFunction==3))
    {
        printf("Shape function must be 0, 1, 2 or 3 \n");
        return -1;
    }
	inputvariables.ShapeFunction = ShapeFunction;
    if (SubsetLength < GridLength)
    {
        printf("Warning: Subset Length smaller than GridLength: you are not using all available data. \n");
    }
    unsigned int propagationfunction = atoi(argv[6]);
    if (propagationfunction != 0 && propagationfunction != 1)
    {
        printf("Propagation Function must be on (1) or off (0)\n");
        return -1;
    }
	inputvariables.propagationfunction = propagationfunction;
	unsigned int ordering = atoi(argv[7]);
    if (ordering != 0 && ordering != 1)
    {
        printf("Ordering must be natural (0) or reverse (1)\n");
        return -1;
    }
	inputvariables.ordering = ordering;
	unsigned int xStart = atoi(argv[8]);
	unsigned int xEnd = atoi(argv[9]);
	unsigned int yStart = atoi(argv[10]);
	unsigned int yEnd = atoi(argv[11]);
	if (xStart > xEnd)
	{
		std::cout << "xStart larger than xEnd" << std::endl;
		return -1;
	}
	inputvariables.xStart = xStart;
	inputvariables.xEnd = xEnd;
	if (yStart > yEnd)
	{
		std::cout << "yStart larger than yEnd" << std::endl;
		return -1;
	}
	inputvariables.yStart = yStart;
	inputvariables.yEnd = yEnd;

    unsigned int offset = 2*((SplineDegree+1)/2) > 5 ? 2*((SplineDegree+1)/2) : 5;
	//std::cout << " offset = " << offset << std::endl;
    offset = GridLength > (SubsetLength/2+offset) ? GridLength-SubsetLength/2 : offset;
    //std::cout << "offset = " << offset << std::endl;
    offset++;
	inputvariables.offset = offset;
    //std::cout << "offset = " << offset << std::endl;
    //std::cout << " SubsetLength/2 = "  << SubsetLength/2 << std::endl;
	if ((yEnd-yStart) < SubsetLength/2+offset)
	{
		std::cout << "Vertical Range of Image is too small with this Subset" << std::endl;
		return -1;
	}
	if ((xEnd-xStart) < SubsetLength/2+offset)
	{
		std::cout << "Horizontal Range of Image is too small with this Subset" << std::endl;
		return -1;
	}
	inputvariables.Number_Of_Threads = atoi(argv[12]);
	if (inputvariables.Number_Of_Threads > std::thread::hardware_concurrency())
	{
		std::cout << "Specified Number of Threads larger than maximum available on this system. Using the system's maximum" << std::endl;
		inputvariables.Number_Of_Threads = std::thread::hardware_concurrency();
	}
	unsigned int MaxPixelYVertical = atoi(argv[13]);
	if (MaxPixelYVertical > (yEnd-yStart))
	{
		std::cout << "Maximum Vertical Pixel Allowed is larger than Vertical Image Size"<< std::endl;
		//std::cout << "TO DO: check if this is an error"<< std::endl;
		//return -1;
	}
	inputvariables.MaxPixelYVertical = MaxPixelYVertical;
	double tolerance = atof(argv[14]);
	if (tolerance < 0)
	{
		std::cout << "Tolerance is negative"<< std::endl;
		return -1;
	}
	if (tolerance >= 1)
	{
		std::cout << "Tolerance is too large" << std::endl;
		return -1;
	}
    // Stopping Criterion
    inputvariables.abs_tolerance_threshold = tolerance;
    inputvariables.rel_tolerance_threshold = tolerance;
	// Minimum Acceptable Correlation Coefficient for Initial Guess
	double minimum_corrcoeff_IG =  atof(argv[15]);
	if (minimum_corrcoeff_IG > 1)
	{
		std::cout << "Minumum Required Correlation Coefficient is too large" << std::endl;
		return -1;
	}
	if (minimum_corrcoeff_IG < -1)
	{
		std::cout << "Minumum Required Correlation Coefficient is too small" << std::endl;
		return -1;
	}
	inputvariables.minimum_corrcoeff_IG = minimum_corrcoeff_IG;
  inputvariables.xStart_ROI = xStart+SubsetLength/2+offset;
  unsigned int xEnd_ROI   = xEnd-SubsetLength/2-offset;
  inputvariables.yStart_ROI = yStart+SubsetLength/2+offset;
  unsigned int yEnd_ROI   = yEnd-SubsetLength/2-offset;
  inputvariables.horx_ROI = xEnd_ROI - inputvariables.xStart_ROI;
  inputvariables.very_ROI = yEnd_ROI - inputvariables.yStart_ROI;

	int BlurSize = atoi(argv[16]);
	if (BlurSize < 0)
	{
			printf("BlurSize must be no (0) or yes (size square filter)\n");
			return -1;
	}
	inputvariables.BlurSize = BlurSize;
	int directionsToInclude = atoi(argv[17]);
	if (directionsToInclude != 1 && directionsToInclude != 2)
	{
			printf("directionsToInclude must be X (1) or XYZ (2)\n");
			return -1;
	}
	inputvariables.directionsToInclude = directionsToInclude;
	double nref = atof(argv[18]);
	if (nref < 1 )
	{
			printf("Reference Index of Refraction is too small\n");
			return -1;
	}
	if (nref > 1.4 )
	{
			printf("Reference Index of Refraction is too large\n");
			return -1;
	}
	inputvariables.nref = nref;

	unsigned int DICneeded = atoi(argv[19]);
  if (DICneeded != 0 && DICneeded != 1)
  {
      printf("DICneeded must be no (0) or yes (1)\n");
      return -1;
  }
	inputvariables.DICNeeded = DICneeded;
	unsigned int Calibrationneeded = atoi(argv[20]);
    if (Calibrationneeded != 0 && Calibrationneeded != 1)
    {
        printf("CalibrationNeeded must be no (0) or yes (1)\n");
        return -1;
    }
	inputvariables.CalibrationNeeded = Calibrationneeded;
	unsigned int Refractionneeded = atoi(argv[21]);
    if (Refractionneeded != 0 && Refractionneeded != 1)
    {
        printf("CalculateRefractionIndex must be no (0) or yes (1)\n");
        return -1;
    }
	inputvariables.CalculateRefractionIndex = Refractionneeded;
	return 1;
}
/*--------------------------------------------------------------------------*/
extern int readImageDataFromFile(cv::Mat &img, cv::Mat &img1, InputVariables &inputvariables)
{
    // First see if a SubFolder 'Averaged' exists
    // If it does, load the csv files from 'Averaged'
    // Otherwise, load the tif files from Folder in path
    std::vector<cv::Mat> data;
    bool readimagesfromcurrentpath = 0;
    std::vector<cv::String> filenames;
    cv::String path1;
    path1 = inputvariables.path+"/Averaged/*.csv";
    cv::glob(path1,filenames,true); // recurse

    if (!filenames.empty())
    {
				// load_matrix is slowest function call here
        cv::Mat I0 = load_matrix(inputvariables.path+"/Averaged/", "A", 0);
        cv::Mat I1 = load_matrix(inputvariables.path+"/Averaged/", "B", 0);


			 //cv::Mat dst;               // dst must be a different Mat
			 //cv::flip(I1, dst, 0);

        data.push_back(I0);
        //data.push_back(dst);
        data.push_back(I1);
    }
    else
    {
				std::cout << "Loading Original Images" << std::endl;
        readimagesfromcurrentpath = 1;
    }
    if (readimagesfromcurrentpath)
    {
        std::vector<cv::String> filenames;
        cv::String path1;
        path1 = inputvariables.path+"/*.tif";
        cv::glob(path1,filenames,true); // recurse
        // Read First Image to Determine Size
        {
            cv::Mat im = cv::imread(filenames[0] , IMREAD_GRAYSCALE );
            if (im.empty())
            {
                std::cout << "No Images Found" << std::endl;
                return -1;
            }
            else
            {
                if (inputvariables.xEnd > static_cast<unsigned int>(im.cols) || inputvariables.yEnd > static_cast<unsigned int>(im.rows))
                {
                    std::cout << "Region of Interest (ROI) larger than Image" << std::endl;
                    std::cout << "Image Size = " << im.size() << std::endl;
                    return -1;
                }
            }
        }
        // Read all tif files in folder
        for (size_t k=0; k<filenames.size(); ++k)
        {
             cv::Mat im = cv::imread(filenames[k] , IMREAD_GRAYSCALE );
             if (im.empty()) continue; //only proceed if sucsessful
             // Change Here
             //im = im(Range(inputvariables.yStart,inputvariables.yEnd),Range(inputvariables.xStart,inputvariables.xEnd));
             data.push_back(im);
             std::cout << filenames[k] << std::endl;
        }
    }

    if (inputvariables.ordering==0)
    {
        (data.at(0)).copyTo(img);
        (data.at(1)).copyTo(img1);
    }
    else if (inputvariables.ordering==1)
    {
        (data.at(1)).copyTo(img);
        (data.at(0)).copyTo(img1);
    }

	// Restrict Template to specified Region
	img = img(Range(inputvariables.yStart,inputvariables.yEnd),Range(inputvariables.xStart,inputvariables.xEnd));

	// Set origin of coordinate system to middle of figure
	// Optical Axis = Middle of Image
	inputvariables.x_p0 = 0.0;//
	inputvariables.y_p0 = 0.0;//1089.0;

	inputvariables.xDiff = inputvariables.xStart;
	inputvariables.yDiff = inputvariables.yStart;

  imwrite(inputvariables.path+"/Template.png", img);
  imwrite(inputvariables.path+"/Deformed.png", img1);

  img.convertTo(img, CV_32F);
  img1.convertTo(img1, CV_32F);
	if (inputvariables.BlurSize != 0)
	{
		GaussianBlur(img,  img,  Size(inputvariables.BlurSize, inputvariables.BlurSize), 0);
		GaussianBlur(img1, img1, Size(inputvariables.BlurSize, inputvariables.BlurSize), 0);
	}
	return 1;
}
/*--------------------------------------------------------------------------*/
extern int checkImageDataFromFile(const cv::String &path, cv::Mat &img, cv::Mat &img1, const unsigned int &xStart, const unsigned int &xEnd, const unsigned int &yStart, const unsigned int &yEnd, const unsigned int &SubsetLength, const unsigned int &offset, unsigned int &xStart_ROI, unsigned int &yStart_ROI, unsigned int &horx_ROI, unsigned int &very_ROI, const unsigned int &ordering)
{
    std::vector<cv::String> filenames;
	cv::String path1;
	path1 = path+"/*.tif";
    cv::glob(path1,filenames,true); // recurse
	// Read First Image to Determine Size
	{
		cv::Mat im = cv::imread(filenames[0] , IMREAD_GRAYSCALE );
		if (im.empty())
		{
			std::cout << "No Images Found" << std::endl;
			return -1;
		}
		else
		{
			if (xEnd > static_cast<unsigned int>(im.cols) || yEnd > static_cast<unsigned int>(im.rows))
			{
				std::cout << "Region of Interest (ROI) larger than Image" << std::endl;
				std::cout << "Image Size = " << im.size() << std::endl;
				return -1;
			}
		}
	}
	std::vector<cv::Mat> data;
	// Read all tif files in folder
    for (size_t k=0; k<filenames.size(); ++k)
    {
         cv::Mat im = cv::imread(filenames[k] , IMREAD_GRAYSCALE );
         if (im.empty()) continue; //only proceed if sucsessful
		 // Change Here
         im = im(Range(yStart,yEnd),Range(xStart,xEnd));
         data.push_back(im);
         std::cout << filenames[k] << std::endl;
    }
	std::cout << std::endl;
    xStart_ROI = xStart+SubsetLength/2+offset;
    unsigned int xEnd_ROI   = xEnd-SubsetLength/2-offset;
    yStart_ROI = yStart+SubsetLength/2+offset;
    unsigned int yEnd_ROI   = yEnd-SubsetLength/2-offset;
    horx_ROI = xEnd_ROI - xStart_ROI;
    very_ROI = yEnd_ROI - yStart_ROI;

    if (ordering==0)
    {
        (data.at(0)).copyTo(img);
        (data.at(1)).copyTo(img1);
    }
    else if (ordering==1)
    {
        (data.at(1)).copyTo(img);
        (data.at(0)).copyTo(img1);
    }
	/*
	// Restrict Template to specified Region
	img = img(Range(yStart,yEnd),Range(xStart,xEnd));
	// We do not want to restrict the Deformed Image (we might miss data)
	// TO DO: Need to rewrite nonlineariteration.cpp to use img1.cols and img1.rows instead of img.cols and img.rows since interpol.c needs width and height of image.
	img1 = img1(Range(yStart,yEnd),Range(xStart,xEnd));
	 * */
    imwrite(path+"/Template.png", img);
    imwrite(path+"/Deformed.png", img1);


    img.convertTo(img, CV_32F);
    img1.convertTo(img1, CV_32F);
	return 1;
}
/*--------------------------------------------------------------------------*/
