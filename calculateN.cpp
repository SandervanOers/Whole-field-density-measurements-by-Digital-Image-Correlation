#include "calculateN.hpp"
/*--------------------------------------------------------------------------*/
static std::vector<cv::Mat> ForwardModel(const cv::Mat &GridX, const cv::Mat &GridY, const double &meanGridX, const double &meanGridY, const cv::Mat &Dx, const cv::Mat &Dy, const double &focal_length, const std::vector<double> &Lengths, const double &Distance_From_Pixels_To_Meters, const std::vector<double> &PlaneDefinition, const double &n_0, const double &n_1, const double &n);
static std::vector<std::vector<double>> getPlanes(const std::vector<double> &Lengths, const std::vector<double> &PlaneDefinition);
static PositionDirection getPositionDirectionPlane1(const cv::Mat &GridX, const cv::Mat &GridY, const double &meanGridX, const double &meanGridY, const cv::Mat &Dx, const cv::Mat &Dy, const double &focal_length, const double &Distance_From_Pixels_To_Meters, const std::vector<double> &PlaneDefinition);
/*--------------------------------------------------------------------------*/
static void calculate_HessianJacobian_ForwardModel(const cv::Mat &GridX, const cv::Mat &GridY, const double &meanGridX, const double &meanGridY, const cv::Mat &Dx, const cv::Mat &Dy, const cv::Mat &W, const double &focal_length, const std::vector<double> &Lengths, const double &Distance_From_Pixels_To_Meters, const std::vector<double> &PlaneDefinition, const double &n_0, const double &n_1, const double &n, const double &nref, cv::Mat &Jacobian, cv::Mat &Hessian, const double &lambda, const int &directionsToInclude);
static std::vector<cv::Mat> computeNumericalDerivativeForwardModel(const cv::Mat &GridX, const cv::Mat &GridY, const double &meanGridX, const double &meanGridY, const cv::Mat &Dx, const cv::Mat &Dy, const double &focal_length, const std::vector<double> &Lengths, const double &Distance_From_Pixels_To_Meters, const std::vector<double> &PlaneDefinition, const double &n_0, const double &n_1, const double &n);
/*--------------------------------------------------------------------------*/
static void calculate_HessianJacobian_ForwardModelwrtn(const cv::Mat &GridX, const cv::Mat &GridY, const double &meanGridX, const double &meanGridY, const cv::Mat &Dx, const cv::Mat &Dy, const double &focal_length, const std::vector<double> &Lengths, const double &Distance_From_Pixels_To_Meters, const std::vector<double> &PlaneDefinition, const double &n_0, const double &n_1, const double &n, const double &nref, cv::Mat &Jacobian, cv::Mat &Hessian, const double &lambda, const int &directionsToInclude);
static std::vector<cv::Mat> computeNumericalDerivativeForwardModelwrtn(const cv::Mat &GridX, const cv::Mat &GridY, const double &meanGridX, const double &meanGridY, const cv::Mat &DX, const cv::Mat &DY, const double &focal_length, const std::vector<double> &Lengths, const double &Distance_From_Pixels_To_Meters, const std::vector<double> &PlaneDefinition, const double &n_0, const double &n_1, const double &n, const int &directionsToInclude);
/*--------------------------------------------------------------------------*/
static double CalculateCostFunction(const cv::Mat &GridX, const cv::Mat &GridY, const double &meanGridX, const double &meanGridY, const cv::Mat &DX, const cv::Mat &DY, const cv::Mat &W, const double &focal_length, const std::vector<double> &Lengths, const double &Distance_From_Pixels_To_Meters, const std::vector<double> &PlaneDefinition, const double &n_0, const double &n_1, const double &n, const double &nref, const int &directionsToInclude);
/*--------------------------------------------------------------------------*/
static double computeXq0(const double &a, const double &c, const double &Lm, const double &f, const double &Sp);
static double computeYq0(const double &b, const double &c, const double &Lm, const double &f, const double &Sp);
static double calculateLf(const double &focal_length, const double &Lm);
/*--------------------------------------------------------------------------*/
static std::vector<cv::Mat> calculateDirectionCosines(const cv::Mat &X, const cv::Mat &Y, const double &L_f);
static std::vector<cv::Mat> calculateIntersectionPlaneLine(const std::vector<cv::Mat> &InitialPosition, const std::vector<cv::Mat> &DirectionCosines, const std::vector<double> &PlaneDefinition);
static PositionDirection calculateIntersectionConstantRefraction(const PositionDirection &InitialPositionDirection, const std::vector<double> &PlaneDefinition, const double &n_initial, const double &n_final);
/*--------------------------------------------------------------------------*/
static double calculateMean(const cv::Mat &Mat)
{
	Scalar meanMatrix = mean(Mat);
	return meanMatrix[0];
}
/*--------------------------------------------------------------------------*/
extern double calculateLf(const double &focal_length, const double &Lm)
{
	return 1.0/(1.0/focal_length-1.0/Lm);
}
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculateDirectionCosines(const cv::Mat &X, const cv::Mat &Y, const double &L_f)
{
  // GridX, GridY, Xq0, Yq0 are now normalized
  // We already subtract: X = (GridX - Xq0)*Sp;
	cv::Mat L_F(X.size(), CV_64FC1, L_f);

	cv::Mat normII = X.mul(X)+Y.mul(Y)+L_f*L_f;
	cv::sqrt(normII, normII);

	cv::Mat alpha = X.mul(-1.0/normII);
	cv::Mat beta = Y.mul(-1.0/normII);
	cv::Mat gamma = L_F.mul(1.0/normII);

	std::vector<cv::Mat> ReturnVector;
  ReturnVector.push_back(alpha);
  ReturnVector.push_back(beta);
  ReturnVector.push_back(gamma);
  return ReturnVector;
}
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> calculateIntersectionPlaneLine(const std::vector<cv::Mat> &InitialPosition, const std::vector<cv::Mat> &DirectionCosines, const std::vector<double> &PlaneDefinition)
{
	// Equation Plane
	// a x + b y + c z + d = 0
	double a = PlaneDefinition[0];
	double b = PlaneDefinition[1];
	double c = PlaneDefinition[2];
	double d = PlaneDefinition[3];
	// Equation Line
	// P = S + I L
	// P = Intersection Point (unknown)
	// S = Initial Position
	// I = Direction Cosines
	// L = Length (unknown)
	cv::Mat S_x = InitialPosition.at(0);
	cv::Mat S_y = InitialPosition.at(1);
	cv::Mat S_z = InitialPosition.at(2);
	cv::Mat alpha = DirectionCosines.at(0);
	cv::Mat beta = DirectionCosines.at(1);
	cv::Mat gamma = DirectionCosines.at(2);

	cv::Mat L(S_x.size(), CV_64FC1, Scalar(0));
	L = (d+a*S_x+b*S_y+c*S_z).mul(-1.0/(a*alpha+b*beta+c*gamma));

	// Calculate Intersection Point
	cv::Mat P_x = S_x+alpha.mul(L);
	cv::Mat P_y = S_y+beta.mul(L);
	cv::Mat P_z = S_z+gamma.mul(L);

	std::vector<cv::Mat> ReturnVector;
  ReturnVector.push_back(P_x);
  ReturnVector.push_back(P_y);
  ReturnVector.push_back(P_z);
  return ReturnVector;
}
/*--------------------------------------------------------------------------*/
extern PositionDirection calculateIntersectionConstantRefraction(const PositionDirection &InitialPositionDirection, const std::vector<double> &PlaneDefinition, const double &n_initial, const double &n_final)
{
	std::vector<cv::Mat> InitialPosition = InitialPositionDirection.Position;
	std::vector<cv::Mat> InitialDirection = InitialPositionDirection.Direction;
	std::vector<cv::Mat> Intersection = calculateIntersectionPlaneLine(InitialPosition, InitialDirection, PlaneDefinition);
	std::vector<cv::Mat> Refracted = SnellsLaw(InitialDirection, PlaneDefinition, n_initial, n_final);
	PositionDirection Return(Intersection, Refracted);
	return Return;
}
/*--------------------------------------------------------------------------*/
static std::vector<cv::Mat> ForwardModel(const cv::Mat &GridX, const cv::Mat &GridY, const double &meanGridX, const double &meanGridY, const cv::Mat &Dx, const cv::Mat &Dy, const double &focal_length, const std::vector<double> &Lengths, const double &Distance_From_Pixels_To_Meters, const std::vector<double> &PlaneDefinition, const double &n_0, const double &n_1, const double &n)
{
  std::vector<std::vector<double>> Planes = getPlanes(Lengths, PlaneDefinition);
  PositionDirection Plane1 = getPositionDirectionPlane1(GridX, GridY, meanGridX, meanGridY, Dx, Dy, focal_length, Distance_From_Pixels_To_Meters, PlaneDefinition);
  PositionDirection Plane2 = calculateIntersectionConstantRefraction(Plane1, Planes[2], n_0, n_1);
  PositionDirection Plane3 = calculateIntersectionConstantRefraction(Plane2, Planes[3], n_1, n);
  PositionDirection Plane4 = calculateIntersectionConstantRefraction(Plane3, Planes[4], n, n_1);
  PositionDirection Plane5 = calculateIntersectionConstantRefraction(Plane4, Planes[5], n_1, n_0);
  std::vector<cv::Mat> PositionPlane6 = calculateIntersectionPlaneLine(Plane5.Position, Plane5.Direction, PlaneDefinition);
	return PositionPlane6;
}
/*--------------------------------------------------------------------------*/
static std::vector<std::vector<double>> getPlanes(const std::vector<double> &Lengths, const std::vector<double> &PlaneDefinition)
{
  double L_c = Lengths[0];
  double L_g = Lengths[1];
  double L_t = Lengths[2];
  double L_s = Lengths[3];
  //double a = PlaneDefinition[0]; // is normalized
  //double b = PlaneDefinition[1];
  //double c = PlaneDefinition[2];
  double d = PlaneDefinition[3];

  double d5 = d-L_s;
  double d4 = d-(L_s+L_g);
  double d3 = d-(L_s+L_g+L_t);
  double d2 = d-(L_s+2.0*L_g+L_t);
  double d1 = d-(L_s+2.0*L_g+L_t+L_c);

  std::vector<double> PlaneDefinition0{0,0,0,0}; // Empty plane
  std::vector<double> PlaneDefinition1{PlaneDefinition[0],PlaneDefinition[1],PlaneDefinition[2],d1}; // Should pass through Origin (0,0,0)
  std::vector<double> PlaneDefinition2{PlaneDefinition[0],PlaneDefinition[1],PlaneDefinition[2],d2};
  std::vector<double> PlaneDefinition3{PlaneDefinition[0],PlaneDefinition[1],PlaneDefinition[2],d3};
  std::vector<double> PlaneDefinition4{PlaneDefinition[0],PlaneDefinition[1],PlaneDefinition[2],d4};
  std::vector<double> PlaneDefinition5{PlaneDefinition[0],PlaneDefinition[1],PlaneDefinition[2],d5};

  std::vector<std::vector<double>> Planes;
  Planes.push_back(PlaneDefinition0);
  Planes.push_back(PlaneDefinition1);
  Planes.push_back(PlaneDefinition2);
  Planes.push_back(PlaneDefinition3);
  Planes.push_back(PlaneDefinition4);
  Planes.push_back(PlaneDefinition5);
  return Planes;
}
/*--------------------------------------------------------------------------*/
static PositionDirection getPositionDirectionPlane1(const cv::Mat &GridX, const cv::Mat &GridY, const double &meanGridX, const double &meanGridY, const cv::Mat &Dx, const cv::Mat &Dy, const double &focal_length, const double &Distance_From_Pixels_To_Meters, const std::vector<double> &PlaneDefinition)
{
  double c = PlaneDefinition[2];
  double d = PlaneDefinition[3];

  double L_m = - d / c;
  double L_f = calculateLf(focal_length, L_m);

  cv::Mat X = GridX.clone();
  cv::Mat Y = GridY.clone();
  X = (X+Dx-meanGridX)*Distance_From_Pixels_To_Meters;
  Y = (Y+Dy-meanGridY)*Distance_From_Pixels_To_Meters;
	std::vector<cv::Mat> InitialDirection = calculateDirectionCosines(X, Y, L_f);

	std::vector<cv::Mat> InitialPosition;
	cv::Mat S_x(GridX.size(), CV_64FC1, Scalar(0));
	cv::Mat S_y(GridX.size(), CV_64FC1, Scalar(0));
	cv::Mat S_z(GridX.size(), CV_64FC1, Scalar(0));
	InitialPosition.push_back(S_x);
	InitialPosition.push_back(S_y);
	InitialPosition.push_back(S_z);

	PositionDirection Plane1(InitialPosition, InitialDirection);
  return Plane1;
}
/*--------------------------------------------------------------------------*/
static std::vector<cv::Mat> computeNumericalDerivativeForwardModelwrtn(const cv::Mat &GridX, const cv::Mat &GridY, const double &meanGridX, const double &meanGridY, const cv::Mat &DX, const cv::Mat &DY, const double &focal_length, const std::vector<double> &Lengths, const double &Distance_From_Pixels_To_Meters, const std::vector<double> &PlaneDefinition, const double &n_0, const double &n_1, const double &n, const int &directionsToInclude)
{
	std::vector<cv::Mat> ReturnMat;
	double h = 1e-8;

	double nplus = n+h;
	double nminus = n-h;
	std::vector<cv::Mat> nder1 = ForwardModel(GridX, GridY, meanGridX, meanGridY, DX, DY, focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, nplus);
	std::vector<cv::Mat> nder2 = ForwardModel(GridX, GridY, meanGridX, meanGridY, DX, DY, focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, nminus);

	cv::Mat nderX = (nder1[0] - nder2[0])/2.0/h;
	ReturnMat.push_back(nderX);
	if (directionsToInclude == 2)
	{
		cv::Mat nderY = (nder1[1] - nder2[1])/2.0/h;
		cv::Mat nderZ = (nder1[2] - nder2[2])/2.0/h;
		ReturnMat.push_back(nderY);
		ReturnMat.push_back(nderZ);
	}

	return ReturnMat;
}
/*--------------------------------------------------------------------------*/
static std::vector<cv::Mat> computeNumericalDerivativeForwardModel(const cv::Mat &GridX, const cv::Mat &GridY, const double &meanGridX, const double &meanGridY, const cv::Mat &DX, const cv::Mat &DY, const double &focal_length, const std::vector<double> &Lengths, const double &Distance_From_Pixels_To_Meters, const std::vector<double> &PlaneDefinition, const double &n_0, const double &n_1, const double &n)
{
  //For the Calibration Measurements the derivatives of n are zero: dn/dx=dx/dy=0
	std::vector<cv::Mat> ReturnMat;
	double h = 1e-8;

	//double L_c = Lengths[0];
	double L_g = Lengths[1];
	double L_t = Lengths[2];
	double L_s = Lengths[3];
	double a = PlaneDefinition[0];
	double b = PlaneDefinition[1];
	double c = PlaneDefinition[2];
	double d = PlaneDefinition[3];
	double L_m = - d / c;
	double normPD = calculateNorm(a, b, c);

	{
	// a derivative
	double aplus = a+h;
	double normPDplus = calculateNorm(aplus, b, c);
	double dplus = -c/normPDplus*L_m;
	double L_cplus = computeLc(c, normPDplus, L_m, L_s, L_g, L_t);
	std::vector<double> Lengthsnewplus(Lengths.begin(), Lengths.end());
	Lengthsnewplus[0] = L_cplus;
	std::vector<double> PlaneDefinitionplus{aplus/normPDplus, b/normPDplus, c/normPDplus, dplus};

	double aminus = a-h;
	double normPDminus = calculateNorm(aminus, b, c);
	double dminus = -c/normPDminus*L_m;
	double L_cminus = computeLc(c, normPDminus, L_m, L_s, L_g, L_t);
	std::vector<double> Lengthsnewminus(Lengths.begin(), Lengths.end());
	Lengthsnewminus[0] = L_cminus;
	std::vector<double> PlaneDefinitionminus{aminus/normPDminus, b/normPDminus, c/normPDminus, dminus};

	std::vector<cv::Mat> ader1 = ForwardModel(GridX, GridY, meanGridX, meanGridY, DX, DY, focal_length, Lengthsnewplus, Distance_From_Pixels_To_Meters, PlaneDefinitionplus, n_0, n_1, n);
	std::vector<cv::Mat> ader2 = ForwardModel(GridX, GridY, meanGridX, meanGridY, DX, DY, focal_length, Lengthsnewminus, Distance_From_Pixels_To_Meters, PlaneDefinitionminus, n_0, n_1, n);

	cv::Mat aderX = (ader1[0] - ader2[0])/2.0/h;
	cv::Mat aderY = (ader1[1] - ader2[1])/2.0/h;
	cv::Mat aderZ = (ader1[2] - ader2[2])/2.0/h;

	ReturnMat.push_back(aderX);
	ReturnMat.push_back(aderY);
	ReturnMat.push_back(aderZ);
	}

	{
	// b derivative
	double bplus = b+h;
	double normPDplus = calculateNorm(a, bplus, c);
	double dplus = -c/normPDplus*L_m;
	double L_cplus = computeLc(c, normPDplus, L_m, L_s, L_g, L_t);
	std::vector<double> Lengthsnewplus(Lengths.begin(), Lengths.end());
	Lengthsnewplus[0] = L_cplus;
	std::vector<double> PlaneDefinitionplus{a/normPDplus, bplus/normPDplus, c/normPDplus, dplus};

	double bminus = b-h;
	double normPDminus = calculateNorm(a, bminus, c);
	double dminus = -c/normPDminus*L_m;
	double L_cminus = computeLc(c, normPDminus, L_m, L_s, L_g, L_t);
	std::vector<double> Lengthsnewminus(Lengths.begin(), Lengths.end());
	Lengthsnewminus[0] = L_cminus;
	std::vector<double> PlaneDefinitionminus{a/normPDminus, bminus/normPDminus, c/normPDminus, dminus};

	std::vector<cv::Mat> bder1 = ForwardModel(GridX, GridY, meanGridX, meanGridY, DX, DY, focal_length, Lengthsnewplus, Distance_From_Pixels_To_Meters, PlaneDefinitionplus, n_0, n_1, n);
	std::vector<cv::Mat> bder2 = ForwardModel(GridX, GridY, meanGridX, meanGridY, DX, DY, focal_length, Lengthsnewminus, Distance_From_Pixels_To_Meters, PlaneDefinitionminus, n_0, n_1, n);

	cv::Mat bderX = (bder1[0] - bder2[0])/2.0/h;
	cv::Mat bderY = (bder1[1] - bder2[1])/2.0/h;
	cv::Mat bderZ = (bder1[2] - bder2[2])/2.0/h;

	ReturnMat.push_back(bderX);
	ReturnMat.push_back(bderY);
	ReturnMat.push_back(bderZ);
	}

	{
	// c derivative
	double cplus = c+h;
	double normPDplus = calculateNorm(a, b, cplus);
	double dplus = -cplus/normPDplus*L_m;
	double L_cplus = computeLc(cplus, normPDplus, L_m, L_s, L_g, L_t);
	std::vector<double> Lengthsnewplus(Lengths.begin(), Lengths.end());
	Lengthsnewplus[0] = L_cplus;
	std::vector<double> PlaneDefinitionplus{a/normPDplus, b/normPDplus, cplus/normPDplus, dplus};

	double cminus = c-h;
	double normPDminus = calculateNorm(a, b, cminus);
	double dminus = -cminus/normPDminus*L_m;
	double L_cminus = computeLc(cminus, normPDminus, L_m, L_s, L_g, L_t);
	std::vector<double> Lengthsnewminus(Lengths.begin(), Lengths.end());
	Lengthsnewminus[0] = L_cminus;
	std::vector<double> PlaneDefinitionminus{a/normPDminus, b/normPDminus, cminus/normPDminus, dminus};

	std::vector<cv::Mat> cder1 = ForwardModel(GridX, GridY, meanGridX, meanGridY, DX, DY, focal_length, Lengthsnewplus, Distance_From_Pixels_To_Meters, PlaneDefinitionplus, n_0, n_1, n);
	std::vector<cv::Mat> cder2 = ForwardModel(GridX, GridY, meanGridX, meanGridY, DX, DY, focal_length, Lengthsnewminus, Distance_From_Pixels_To_Meters, PlaneDefinitionminus, n_0, n_1, n);

	cv::Mat cderX = (cder1[0] - cder2[0])/2.0/h;
	cv::Mat cderY = (cder1[1] - cder2[1])/2.0/h;
	cv::Mat cderZ = (cder1[2] - cder2[2])/2.0/h;

	ReturnMat.push_back(cderX);
	ReturnMat.push_back(cderY);
	ReturnMat.push_back(cderZ);
	}

  {
	// Lm derivative
	double L_mplus = L_m + h;
	double dplus = -c/normPD*L_mplus;
	double L_cplus = computeLc(c, normPD, L_mplus, L_s, L_g, L_t);
	std::vector<double> Lengthsnewplus(Lengths.begin(), Lengths.end());
	Lengthsnewplus[0] = L_cplus;
	std::vector<double> PlaneDefinitionplus{a/normPD, b/normPD, c/normPD, dplus};

	double L_mminus = L_m-h;
	double dminus = -c/normPD*L_mminus;
	double L_cminus = computeLc(c, normPD, L_mminus, L_s, L_g, L_t);
	std::vector<double> Lengthsnewminus(Lengths.begin(), Lengths.end());
	Lengthsnewminus[0] = L_cminus;
	std::vector<double> PlaneDefinitionminus{a/normPD, b/normPD, c/normPD, dminus};

	std::vector<cv::Mat> Lmder1 = ForwardModel(GridX, GridY, meanGridX, meanGridY, DX, DY, focal_length, Lengthsnewplus, Distance_From_Pixels_To_Meters, PlaneDefinitionplus, n_0, n_1, n);
	std::vector<cv::Mat> Lmder2 = ForwardModel(GridX, GridY, meanGridX, meanGridY, DX, DY, focal_length, Lengthsnewminus, Distance_From_Pixels_To_Meters, PlaneDefinitionminus, n_0, n_1, n);

	cv::Mat LmderX = (Lmder1[0] - Lmder2[0])/2.0/h;
	cv::Mat LmderY = (Lmder1[1] - Lmder2[1])/2.0/h;
	cv::Mat LmderZ = (Lmder1[2] - Lmder2[2])/2.0/h;

	ReturnMat.push_back(LmderX);
	ReturnMat.push_back(LmderY);
	ReturnMat.push_back(LmderZ);
	}
	return ReturnMat;
}
/*--------------------------------------------------------------------------*/
static void calculate_HessianJacobian_ForwardModelwrtn(const cv::Mat &GridX, const cv::Mat &GridY, const double &meanGridX, const double &meanGridY, const cv::Mat &Dx, const cv::Mat &Dy, const double &focal_length, const std::vector<double> &Lengths, const double &Distance_From_Pixels_To_Meters, const std::vector<double> &PlaneDefinition, const double &n_0, const double &n_1, const double &n, const double &nref, cv::Mat &Jacobian, cv::Mat &Hessian, const double &lambda, const int &directionsToInclude)
{
  // Measurement with unknown n: dn/dx and dn/dy can be nonzero -> J
  // Turned the y and z terms off
	std::vector<cv::Mat> X1dn = computeNumericalDerivativeForwardModelwrtn(GridX, GridY, meanGridX, meanGridY, Dx, Dy, focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, n, directionsToInclude);

	cv::Mat DX0(GridX.size(), CV_64FC1, Scalar(0));
	cv::Mat DY0(GridX.size(), CV_64FC1, Scalar(0));
	std::vector<cv::Mat> X1 = ForwardModel(GridX, GridY, meanGridX, meanGridY, Dx, Dy, focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0,  n_1, n);
	std::vector<cv::Mat> X0 = ForwardModel(GridX, GridY, meanGridX, meanGridY, DX0, DY0, focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0,  n_1, nref);

	int numberofrows = GridX.rows*GridX.cols;
	cv::Mat Jx = X1dn[0].reshape(0, numberofrows);
	Hessian = Jx.t()*Jx;

	cv::Mat DX = X1[0] - X0[0];
	cv::Mat dy_x = -DX.reshape(0,numberofrows);
	Jacobian = Jx.t()*dy_x;

	if (directionsToInclude == 2)
	{
		cv::Mat Jy = X1dn[1].reshape(0, numberofrows);
		cv::Mat Jz = X1dn[2].reshape(0, numberofrows);
		cv::Mat DY = X1[1] - X0[1];
		cv::Mat DZ = X1[2] - X0[2];
		cv::Mat dy_y = -DY.reshape(0,numberofrows);
		cv::Mat dy_z = -DZ.reshape(0,numberofrows);
		Jacobian = Jacobian + Jy.t()*dy_y + Jz.t()*dy_z;
		Hessian = Hessian + Jy.t()*Jy + Jz.t()*Jz;
	}

	for (unsigned int i = 0; i < 1; i++)
	{
		Hessian.at<double>(i,i) *= (1.0+lambda);
	}
}
/*--------------------------------------------------------------------------*/
static void calculate_HessianJacobian_ForwardModel(const cv::Mat &GridX, const cv::Mat &GridY, const double &meanGridX, const double &meanGridY, const cv::Mat &Dx, const cv::Mat &Dy, const cv::Mat &W, const double &focal_length, const std::vector<double> &Lengths, const double &Distance_From_Pixels_To_Meters, const std::vector<double> &PlaneDefinition, const double &n_0, const double &n_1, const double &n, const double &nref, cv::Mat &Jacobian, cv::Mat &Hessian, const double &lambda, const int &directionsToInclude)
{
  // We use only x-direction
	std::vector<cv::Mat> X1da = computeNumericalDerivativeForwardModel(GridX, GridY, meanGridX, meanGridY, Dx, Dy, focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, n);
	cv::Mat DX0(GridX.size(), CV_64FC1, Scalar(0));
	cv::Mat DY0(GridX.size(), CV_64FC1, Scalar(0));
	std::vector<cv::Mat> X0da = computeNumericalDerivativeForwardModel(GridX, GridY, meanGridX, meanGridY, DX0, DY0, focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, nref);
  std::vector<cv::Mat> JacobianFull;
  // Current calculation is wrong: Weight is applied twice in the Hessian, once in the RHS (Jacobian). Not a problem now since the weights are 0 and 1.
	for (unsigned int i = 0; i < 12; i++)
	{
		cv::Mat JXa = X1da[i]-X0da[i];
    JXa = JXa.mul(W);
		JacobianFull.push_back(JXa);
	}
	int numberofrows = GridX.rows*GridX.cols;
	cv::Mat Jx;

	cv::Mat dfxda = JacobianFull[0].reshape(0, numberofrows);
	cv::Mat dfxdb = JacobianFull[3].reshape(0, numberofrows);
	cv::Mat dfxdc = JacobianFull[6].reshape(0, numberofrows);
	cv::Mat dfxdLm = JacobianFull[9].reshape(0, numberofrows);
	cv::Mat matArrayX[] = {dfxda, dfxdb, dfxdc, dfxdLm};
	cv::hconcat( matArrayX, 4, Jx);

  Hessian = Jx.t()*Jx;

	std::vector<cv::Mat> X1 = ForwardModel(GridX, GridY, meanGridX, meanGridY, Dx , Dy , focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0,  n_1, n);
	std::vector<cv::Mat> X0 = ForwardModel(GridX, GridY, meanGridX, meanGridY, DX0, DY0, focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0,  n_1, nref);

	cv::Mat DX = X1[0] - X0[0];
	cv::Mat dy_x = -DX.reshape(0,numberofrows);
	Jacobian = Jx.t()*dy_x;

	if (directionsToInclude == 2)
	{
		cv::Mat Jy, Jz;

		cv::Mat dfyda = JacobianFull[1].reshape(0, numberofrows);
		cv::Mat dfydb = JacobianFull[4].reshape(0, numberofrows);
		cv::Mat dfydc = JacobianFull[7].reshape(0, numberofrows);
		cv::Mat dfydLm = JacobianFull[10].reshape(0, numberofrows);
		cv::Mat matArrayY[] = {dfyda, dfydb, dfydc, dfydLm};
		cv::hconcat( matArrayY, 4, Jy);

		cv::Mat dfzda = JacobianFull[2].reshape(0, numberofrows);
		cv::Mat dfzdb = JacobianFull[5].reshape(0, numberofrows);
		cv::Mat dfzdc = JacobianFull[8].reshape(0, numberofrows);
		cv::Mat dfzdLm = JacobianFull[11].reshape(0, numberofrows);
		cv::Mat matArrayZ[] = {dfzda, dfzdb, dfzdc, dfzdLm};
		cv::hconcat( matArrayZ, 4, Jz);

		cv::Mat DY = X1[1] - X0[1];
		cv::Mat DZ = X1[2] - X0[2];
		cv::Mat dy_y = -DY.reshape(0,numberofrows);
		cv::Mat dy_z = -DZ.reshape(0,numberofrows);
		Jacobian = Jacobian + Jy.t()*dy_y + Jz.t()*dy_z;

		Hessian = Hessian + Jy.t()*Jy + Jz.t()*Jz;
	}
	for (unsigned int i = 0; i < 4; i++)
	{
		Hessian.at<double>(i,i) *= (1.0+lambda);
	}
}
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> CalculateN(const cv::Mat &GridX, const cv::Mat &GridY, const double &meanGridX, const double &meanGridY, const cv::Mat &Dx, const cv::Mat &Dy, const double &focal_length, const std::vector<double> &Lengths, const double &Distance_From_Pixels_To_Meters, const std::vector<double> &PlaneDefinition, const double &n_0, const double &n_1, const double &nref, const int &directionsToInclude)
{
	double n_ini = 1.0;
	double n_old = 1.0;
	double n_lower_bound = 1.3;
	double n_upper_bound = 1.4;
	double n_around = 0.0;
	if (nref == 1)
	{	// reference fluid is air
		// n range = 1.333 - 1.4
			n_ini = 1.337;
			n_old = 1.337;
			n_around = 0.0;
			n_lower_bound = 1.3;
			n_upper_bound = 1.4;
	}
	else
	{ // reference fluid is water
		// n = 1.333 + (0 - 0.666) (nref = 1.333)
		// n = n_around + dn
			n_ini = 0.0;
			n_old = 0.0;
			n_around = nref;
			n_lower_bound = -0.033;
			n_upper_bound = 0.666;
	}

	double abs_tolerance_threshold = 1e-12;
	double rel_tolerance_threshold = 1e-12;
	double max_val_nu = 1e7;
	unsigned int max_iterations = 1e4;

	std::vector<double> LengthsMiddleTank{experimentalsetupvariables.L_c, experimentalsetupvariables.L_g, experimentalsetupvariables.L_t/2.0, 0};

	cv::Mat n_field(GridX.size(), CV_64FC1, Scalar(0));
	cv::Mat ksi(GridX.size(), CV_64FC1, Scalar(0));
	cv::Mat eta(GridX.size(), CV_64FC1, Scalar(0));
	cv::Mat dzeta(GridX.size(), CV_64FC1, Scalar(0));
	double outputcount = 1;

	for( int i = 0; i < GridX.rows; ++i)
	{
		for( int j = 0; j < GridX.cols; ++j)
		{
			unsigned int iterations = 1;
			double abs_tolerance = 1;
			double rel_tolerance = 1;
			double nu = 2.0;
			double lambda = 1e-10;
			cv::Mat GX(1, 1, CV_64F, GridX.at<double>(i,j));
			cv::Mat GY(1, 1, CV_64F, GridY.at<double>(i,j));
			cv::Mat dx(1, 1, CV_64F, Dx.at<double>(i,j));
			cv::Mat dy(1, 1, CV_64F, Dy.at<double>(i,j));
			cv::Mat W(1, 1, CV_64FC1, Scalar(1));
			double Sold = CalculateCostFunction(GX, GY, meanGridX, meanGridY, dx, dy, W, focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, n_around+n_ini, nref, directionsToInclude);

			double Snew = Sold;
			double n = n_ini;
			double nnew = n_ini;
			cv::Mat Hessian;
			cv::Mat Jacobian;
			calculate_HessianJacobian_ForwardModelwrtn(GX, GY, meanGridX, meanGridY, dx, dy, focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, n_around+n_ini, nref, Jacobian, Hessian, lambda, directionsToInclude);
			while (iterations < max_iterations && nu < max_val_nu && (abs_tolerance > abs_tolerance_threshold || rel_tolerance > rel_tolerance_threshold))
			{
					iterations++;
					cv::Mat delta_alpha(1,1,CV_64F);
					cv::solve(Hessian, Jacobian, delta_alpha, cv::DECOMP_CHOLESKY);
					abs_tolerance = cv::sum(cv::abs(delta_alpha)).val[0];
					rel_tolerance = cv::sum(cv::abs(delta_alpha)).val[0]/(n_ini);
					if (isnan(abs_tolerance))
					{
						std::cout << "NaN Detected"<< std::endl;
						abs_tolerance = 1;
						rel_tolerance = 1;
						delta_alpha.at<double>(0) = 0;
					}
					nnew = n + delta_alpha.at<double>(0);
					Snew = CalculateCostFunction(GX, GY, meanGridX, meanGridY, dx, dy, W, focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, n_around+nnew, nref, directionsToInclude);
					if ((Snew < Sold ) && nnew >= n_lower_bound && nnew <= n_upper_bound)
					{
						n = nnew;
						Sold = Snew;
						nu = 2;
						lambda /= 3.0;
						calculate_HessianJacobian_ForwardModelwrtn(GX, GY, meanGridX, meanGridY, dx, dy, focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, n_around+n, nref, Jacobian, Hessian, lambda, directionsToInclude);
					}
					else
					{
						// No Improvement
						// Increase lambda => Less Gauss-Newton, more Gradient Search
						double lambda_new = lambda*nu;
						nu *= 2.0;
						// Scale Diagonal of Hessian
						// Jacobian stays the same
						Hessian.at<double>(0,0) *= (1.0+lambda_new)/(1.0+lambda);
						lambda = lambda_new;
					}
			}
			std::vector<cv::Mat> X60 = ForwardModel(GX, GY, meanGridX, meanGridY, dx, dy, focal_length, LengthsMiddleTank, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, nref);
			ksi.at<double>(i,j) = X60[0].at<double>(0,0);
			eta.at<double>(i,j) = X60[1].at<double>(0,0);
			dzeta.at<double>(i,j) = X60[2].at<double>(0,0);
			n_field.at<double>(i,j) = n_around+n;
			n_ini = n_old;
		}
		double computed = (double)i/GridX.rows*100.0;
		if (computed>outputcount)
		{
			outputcount++;
			outputcount++;
			outputcount++;
		  std::cout << "Points Computed: " << std::setprecision(2) << computed << "%" << std::endl;
      std::cout << std::setprecision(6);
		}
	}
  std::cout << "mean n = " << std::setprecision(6) << calculateMean(n_field) << std::endl;
  std::cout << "L2 error = " << norm(n_field-1.337, NORM_L2, noArray()) << std::endl;
  std::cout << "L2 error = " << std::setprecision(10) << std::sqrt(1.0/n_field.rows/n_field.cols) * norm(n_field-1.337, NORM_L2, noArray()) << std::endl;

	std::vector<cv::Mat> ReturnMat;
	ReturnMat.push_back(n_field);
	ReturnMat.push_back(ksi);
	ReturnMat.push_back(eta);
	ReturnMat.push_back(dzeta);
	return ReturnMat;
}
/*--------------------------------------------------------------------------*/
extern std::vector<double> Calibration(const cv::Mat &GridX, const cv::Mat &GridY, const cv::Mat &Dx, const cv::Mat &Dy, const cv::Mat &CorrelationCoefficient, const double &focal_length, const std::vector<double> &Lengths, const double &Distance_From_Pixels_To_Meters, const double &n_0, const double &n_1, const double &n, const double &nref, const std::string &path, const double &corr_cut_off, const int &directionsToInclude)
{
	double abs_tolerance = 1;
	double rel_tolerance = 1;
	double abs_tolerance_threshold = 1e-14;
	double rel_tolerance_threshold = 1e-14;
	double max_val_nu = 1e7;
	unsigned int max_iterations = 2e5;
	unsigned int iterations = 1;

	double L_g = Lengths[1];
	double L_t = Lengths[2];
	double L_s = Lengths[3];

	cv::Mat DX0(GridX.size(), CV_64FC1, Scalar(0));
	cv::Mat DY0(GridX.size(), CV_64FC1, Scalar(0));

	double nu  = 2.0;
	double lambda = 1e-10;
  // Initial Guess
	//////////////////////////////////////////////////////////////////////////////
	// Initialize vector q = (a, b, c, L_m) to your own values                  //
	//////////////////////////////////////////////////////////////////////////////
  double L_m = 1.188;          // in meters
  double a = -0.293;
  double b = -0.015;
  double c = -0.955;


  double normPD = calculateNorm(a, b, c);
  a = a/normPD;
  b = b/normPD;
  c = c/normPD;
  normPD = calculateNorm(a, b, c);

  double Sp = Distance_From_Pixels_To_Meters;
  double Xq0 = computeXq0(a, c, L_m, focal_length, Sp);  // in pixels
  double Yq0 = computeYq0(b, c, L_m, focal_length, Sp);  // in pixels

	double aold = a;
	double bold = b;
	double cold = c;
	double L_mold = L_m;
  double Xq0old = Xq0;
  double Yq0old = Yq0;

	double d = -c*L_m;
	double L_c = computeLc(c, normPD, L_m, L_s, L_g, L_t);
	std::vector<double> LengthsChanging(Lengths);
	LengthsChanging[0] = L_c;
	std::vector<double> PlaneDefinition{a/normPD, b/normPD, c/normPD, d/normPD};

	std::cout << "L_c = " << L_c << std::endl;
	std::cout << "L_m = " << L_m << " >  Ltot = " << L_c + 2*L_g + L_s +L_t << std::endl;
	double anew = a;
	double bnew = b;
	double cnew = c;
	double dnew = d;
	double L_mnew = L_m;
  double Xq0new = Xq0;
  double Yq0new = Yq0;
	double L_cnew = L_c;
	std::vector<double> Lengthsnew(LengthsChanging.begin(), LengthsChanging.end());
	std::cout << "Lengthsnew: " ;
	 for (auto i = Lengthsnew.begin(); i < Lengthsnew.end(); i++)
	{
		std::cout << *i << " " ;
	}
	std::cout << std::endl;
	std::cout << std::endl;
	std::vector<double> PlaneDefinitionnew{a/normPD, b/normPD, c/normPD, d/normPD};
	std::cout << "Original Plane: " ;
	 for (auto i = PlaneDefinition.begin(); i < PlaneDefinition.end(); i++)
	{
		std::cout << *i << " " ;
	}
	std::cout << std::endl;

	int numberofrows = GridX.rows*GridX.cols;
	cv::Mat W;
	CorrelationCoefficient.copyTo(W);
	for(int i = 0; i < W.rows; i++)
	{
		double* Wi = W.ptr<double>(i);
		for(int j = 0; j < W.cols; j++)
			if (Wi[j]<corr_cut_off)
			{
				Wi[j] = 0.0;
			}
	}
  double Sold = CalculateCostFunction(GridX, GridY, Xq0, Yq0, Dx, Dy, W, focal_length, LengthsChanging, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, n, nref, directionsToInclude);

	double Snew = Sold;
	double S_ini = Sold;
	std::cout << "S old = " << Sold << std::endl;

	cv::Mat Hessian;
	cv::Mat Jacobian;
	calculate_HessianJacobian_ForwardModel(GridX, GridY, Xq0, Yq0, Dx, Dy, W, focal_length, LengthsChanging, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, n, nref, Jacobian, Hessian, lambda, directionsToInclude);
	std::cout << "Size Jacobian = " << Jacobian.size() << std::endl;
	std::cout << "Size Hessian = " << Hessian.size() << std::endl;

	while (iterations < max_iterations && nu < max_val_nu && (abs_tolerance > abs_tolerance_threshold || rel_tolerance > rel_tolerance_threshold))
	{
            iterations++;
            cv::Mat delta_alpha(4,1,CV_64F);
            cv::solve(Hessian, Jacobian, delta_alpha, cv::DECOMP_CHOLESKY);
            abs_tolerance = cv::sum(cv::abs(delta_alpha)).val[0];
            rel_tolerance = cv::sum(cv::abs(delta_alpha)).val[0]/(abs(aold)+abs(bold)+abs(cold)+abs(L_mold));
				if (isnan(abs_tolerance))
				{
					std::cout << "NaN Detected"<< std::endl;
					abs_tolerance = 1;
					rel_tolerance = 1;
					delta_alpha.at<double>(0) = 0;
					delta_alpha.at<double>(1) = 0;
					delta_alpha.at<double>(2) = 0;
					delta_alpha.at<double>(3) = 0;
				}

		{
			anew = a+delta_alpha.at<double>(0);
			bnew = b+delta_alpha.at<double>(1);
			cnew = c+delta_alpha.at<double>(2);
			L_mnew = L_m+delta_alpha.at<double>(3);
      double normPDnew = calculateNorm(anew, bnew, cnew);
      anew = anew/normPDnew;
      bnew = bnew/normPDnew;
      cnew = cnew/normPDnew;
      normPDnew = calculateNorm(anew, bnew, cnew);
      Xq0new = computeXq0(anew, cnew, L_mnew, focal_length, Sp);  // in pixels
      Yq0new = computeYq0(bnew, cnew, L_mnew, focal_length, Sp);  // in pixels
			dnew = -cnew*L_mnew;
			double pd[] = {anew/normPDnew, bnew/normPDnew, cnew/normPDnew, dnew/normPDnew};
			PlaneDefinitionnew.assign(pd, pd+4);
    	L_cnew = computeLc(cnew, normPDnew, L_mnew, L_s, L_g, L_t);
			Lengthsnew[0] = L_cnew;
      Snew = CalculateCostFunction(GridX, GridY, Xq0new, Yq0new, Dx, Dy, W, focal_length, Lengthsnew, Distance_From_Pixels_To_Meters, PlaneDefinitionnew, n_0, n_1, n, nref, directionsToInclude);
		}
		if ((Snew < Sold ) && L_cnew > 0 && L_mnew > L_cnew + 2*L_g+L_s+L_t)
		{
			// Improvement
			Sold = Snew;
			a = anew;
			b = bnew;
			c = cnew;
			d = dnew;
			L_m = L_mnew;
			L_c = L_cnew;
			Xq0 = Xq0new;
			Yq0 = Yq0new;
			PlaneDefinition.assign(PlaneDefinitionnew.begin(), PlaneDefinitionnew.end());
			nu = 2;
			lambda /= 3.0;
			LengthsChanging[0] = L_c;
			std::cout << a << " " << b << " " << c << " " << L_m << " " << Xq0 << " " << Yq0 << " " << Sold << " " << lambda << " " << iterations << "\n";

			calculate_HessianJacobian_ForwardModel(GridX, GridY, Xq0, Yq0, Dx, Dy, W, focal_length, LengthsChanging, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, n, nref, Jacobian, Hessian, lambda, directionsToInclude);
		}
		else
		{
			// No Improvement
			// Increase lambda => Less Gauss-Newton, more Gradient Search
			double lambda_new = lambda*nu;
			nu *= 2.0;

			// Scale Diagonal of Hessian
			// Jacobian stays the same
			for (unsigned int i = 0; i < 4; i++)
			{
				Hessian.at<double>(i,i) *= (1.0+lambda_new)/(1.0+lambda);
			}
			lambda = lambda_new;
		}
	}
	std::cout << std::endl;
	std::cout << std::endl;
	if (iterations >= max_iterations)
	{
		std::cout << "max iterations reached" << std::endl;
	}
	if (nu >= max_val_nu)
	{
		std::cout << "max nu reached" << std::endl;
	}
	if (abs_tolerance <= abs_tolerance_threshold)
	{
		std::cout << "abs tolerance reached" << std::endl;
	}
	if (rel_tolerance <= rel_tolerance_threshold)
	{
		std::cout << "rel tolerance reached" << std::endl;
	}
	std::cout << std::endl;
	std::cout << iterations << " " << nu  << " " << abs_tolerance  << " " << rel_tolerance  << std::endl;
	std::cout << std::endl << "\033[1;32mResult Calibration\033[0m\n" << std::endl<< std::endl;
	std::cout << "aold = " <<   std::setprecision(4)  << aold << ", bold = " << bold << ", cold = " << cold << ", Lmold = " << L_mold << std::endl;
	std::cout << "S old = " << std::setprecision(2) << std::scientific << S_ini << std::endl;
	std::cout << "a = " << std::setprecision(8) << a << ", b = " << b << ", c = " << c << ", Lm = " << L_m << std::endl;
	std::cout << "S = " << std::setprecision(2) << std::scientific << Sold << std::endl << std::endl;
	std::cout << "L_c = " << std::setprecision(3)  << L_c << std::endl;
	std::cout << "L_tot = " << std::setprecision(3) << L_c+2.0*L_g+L_t+L_s << std::endl;

	std::cout << "Mean S old = " << S_ini/numberofrows << std::endl;
	std::cout << "Mean S new = " << Sold/numberofrows << std::endl;

  std::cout << "old <a,b,c> = <" << aold << ", " << bold << ", " << cold << ">" << std::endl;
  std::cout << "side should be " << -12.0/sqrt(12*12+34*34) << ", 0, " << -34.0/sqrt(12*12+34*34)<< std::endl;
	std::cout << std::endl;
	std::cout << "Xq0 old = " << Xq0old << std::endl;
	std::cout << "Yq0 old = " << Yq0old << std::endl;
	std::cout << "Xq0 new = " << Xq0 << std::endl;
	std::cout << "Yq0 new = " << Yq0 << std::endl;
	std::cout << std::endl;

	// Normalize plane normal vector
	normPD = calculateNorm(a, b, c);
	a = a/normPD;
	b = b/normPD;
	c = c/normPD;
	std::vector<double> Returnvector;
	Returnvector.push_back(a);
	Returnvector.push_back(b);
	Returnvector.push_back(c);
	Returnvector.push_back(L_m);
	Returnvector.push_back(Xq0);
	Returnvector.push_back(Yq0);
	std::ofstream myfilen;
  myfilen.open(path+"/CalibrationFigures.csv");
  myfilen << "a, b, c, Lm, Xq0, Yq0, n, meanS"<< std::endl;
  myfilen << std::setprecision(8) << a << ", " << b << ", " << c << ", " << L_m << ", " << Xq0 << ", " << Yq0 << ", " << n << ", " << Sold/numberofrows << std::endl;
	myfilen.close();
	return Returnvector;
}
/*--------------------------------------------------------------------------*/
extern void CalibrationFigures(const cv::Mat &GridX, const cv::Mat &GridY, const cv::Mat &Dx, const cv::Mat &Dy, const cv::Mat &CorrelationCoefficient, const double &focal_length, const std::vector<double> &Lengths, const double &Distance_From_Pixels_To_Meters, const double &n_0, const double &n_1, const double &n, const double &nref, const std::string &path, const double &corr_cut_off)
{
  // Check this function before using it

	double L_g = Lengths[1];
	double L_t = Lengths[2];
	double L_s = Lengths[3];
	// Grid Search
	unsigned int iterations = 300;

	cv::Mat W(GridX.rows, GridX.cols, CV_64F, cv::Scalar(1));
	for (int i = 0; i < GridX.rows; i++)
	{
		for (int j = 0; j < GridX.cols; j++)
		{
			if (CorrelationCoefficient.at<double>(i,j) < corr_cut_off)
			{
				W.at<double>(i,j) = 0.0;
			}
		}
	}

	// Loop over meanGridX, meanGridY, Lm
	double Xp0, Yp0, Lm, a, b, c, d;

  // Initial Guess
  //  Tueday Side
	double Lm_ini = 1;
	double Xp0_ini = 0.0;
	double Yp0_ini = 0.0;
  double a_ini = -10;
  double b_ini = 0;
  double c_ini = -34;

	Lm = Lm_ini;
	Xp0 = Xp0_ini;
	Yp0 = Yp0_ini;
	a = a_ini;
	b = b_ini;
	c = c_ini;
	d = -c*Lm;
  double normPD = calculateNorm(a, b, c);
	std::vector<double> PlaneDefinition{a/normPD, b/normPD, c/normPD, d/normPD};

	double Lc_ini = computeLc(c, normPD, Lm, L_s, L_g, L_t);
	double Lc = Lc_ini;

	std::vector<double> Lengthsnew(Lengths.begin(), Lengths.end());
	Lengthsnew[0] = Lc;
  double Xq0 = computeXq0(a, c, Lm, focal_length, Distance_From_Pixels_To_Meters);  // in pixels
  double Yq0 = computeYq0(b, c, Lm, focal_length, Distance_From_Pixels_To_Meters);  // in pixels

  std::ofstream myfileXp0Yp0;
  myfileXp0Yp0.open(path+"/SfileXp0Yp0.csv");
  myfileXp0Yp0 << "Xp0 Yp0 Lm S Lc Ltot a b c"<< std::endl;
  for (unsigned int i = 0; i <= iterations; i++)
  {
    Xp0 = -2500.0 + 1500.0*(double)i/iterations;
    for (unsigned int j = 0; j <= iterations; j++)
    {
      Yp0 = -1800.0 + 3600.0*(double)j/iterations;
      double S = CalculateCostFunction(GridX+Xp0, GridY+Yp0, Xq0, Yq0, Dx, Dy, W, focal_length, Lengthsnew, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, n, nref, 1);
			myfileXp0Yp0 << Xp0 << " " << Yp0 << " " << Lm << " " << S << " " << Lc << " " << (Lc+L_s+L_t+2.0*L_g) << " " << a << " " << b << " " << c << std::endl;
    }
  }
myfileXp0Yp0.close();
std::cout << std::endl << "\033[1;32m Xp0 Yp0 Calibration Complete\033[0m\n" << std::endl;
}
/*--------------------------------------------------------------------------*/
extern double calculateNorm(const double &a, const double &b, const double &c)
{
	return sqrt(a*a+b*b+c*c);
}
/*--------------------------------------------------------------------------*/
static double CalculateCostFunction(const cv::Mat &GridX, const cv::Mat &GridY, const double &meanGridX, const double &meanGridY, const cv::Mat &DX, const cv::Mat &DY, const cv::Mat &W, const double &focal_length, const std::vector<double> &Lengths, const double &Distance_From_Pixels_To_Meters, const std::vector<double> &PlaneDefinition, const double &n_0, const double &n_1, const double &n, const double &nref, const int &directionsToInclude)
{
  // Calculates Cost Function
  // directionsToInclude: 1 = only X, 2 = X, Y and Z
	cv::Mat DX0(GridX.size(), CV_64FC1, Scalar(0));
	cv::Mat DY0(GridX.size(), CV_64FC1, Scalar(0));
	std::vector<cv::Mat> X61 = ForwardModel(GridX, GridY, meanGridX, meanGridY, DX, DY, focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, n);
	std::vector<cv::Mat> X60 = ForwardModel(GridX, GridY, meanGridX, meanGridY, DX0, DY0, focal_length, Lengths, Distance_From_Pixels_To_Meters, PlaneDefinition, n_0, n_1, nref);
	cv::Mat DX6 = X61[0] - X60[0];
  cv::Mat WDX6 = DX6.mul(W);
  Scalar SS = DX6.dot(WDX6);
  if (directionsToInclude == 2)
	{
    cv::Mat DY6 = X61[1] - X60[1];
    cv::Mat WDY6 = DY6.mul(W);
	  cv::Mat DZ6 = X61[2] - X60[2];
    cv::Mat WDZ6 = DZ6.mul(W);
    Scalar S2 = DY6.dot(WDY6) + DZ6.dot(WDZ6);
    SS = SS +  S2;
  }
	return SS[0];
}
/*--------------------------------------------------------------------------*/
double computeXq0(const double &a, const double &c, const double &Lm, const double &f, const double &Sp)
{
  double Lf = calculateLf(f, Lm);
  return -a/c*Lf/Sp; // in pixels
}
double computeYq0(const double &b, const double &c, const double &Lm, const double &f, const double &Sp)
{
  double Lf = calculateLf(f, Lm);
  return -b/c*Lf/Sp; // in pixels
}
extern double computeLc(const double &c, const double &normPD, const double &L_m, const double &L_s, const double &L_g, const double &L_t)
{
  return -c/normPD*L_m-L_s-2.0*L_g-L_t;
}
