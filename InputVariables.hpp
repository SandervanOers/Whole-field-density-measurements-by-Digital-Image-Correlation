#ifndef H_InputVariables
#define H_InputVariables
struct InputVariables {

	std::string path;
	unsigned int SplineDegree;
	unsigned int SubsetLength;
	unsigned int GridLength;
	unsigned int ShapeFunction;
	unsigned int ReliabilityGuidedDIC;
	unsigned int propagationfunction;
	unsigned int ordering;
	unsigned int xStart;
	unsigned int xEnd;
	unsigned int yStart;
	unsigned int yEnd;
	unsigned int offset;
	unsigned int Number_Of_Threads;
	unsigned int MaxPixelYVertical;
	double abs_tolerance_threshold; // not used
	double rel_tolerance_threshold; // not used
	double minimum_corrcoeff_IG;
	unsigned int xDiff;
	unsigned int yDiff;
	unsigned int x_p0;
	unsigned int y_p0;

	unsigned int xStart_ROI;
	unsigned int yStart_ROI;
	unsigned int horx_ROI;
	unsigned int very_ROI;

	unsigned int DICNeeded;
	unsigned int CalibrationNeeded;
	unsigned int CalculateRefractionIndex;
	unsigned int BlurSize;
	unsigned int directionsToInclude;

	double nref;
};
#endif
