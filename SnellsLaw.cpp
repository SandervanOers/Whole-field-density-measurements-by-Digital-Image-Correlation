#include "SnellsLaw.hpp"

/*--------------------------------------------------------------------------*/
static cv::Mat signum(cv::Mat src)
{
    cv::Mat dst = (src >= 0) & 1;
    dst.convertTo(dst,CV_64F, 2.0, -1.0);
    return dst;
}
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> SnellsLaw(const std::vector<cv::Mat> &DirectionCosines, const std::vector<double> &PlaneDefinition, const double &n_0, const double &n_1)
{
	// Normal vector of Plane
	double a = PlaneDefinition[0];
	double b = PlaneDefinition[1];
	double c = PlaneDefinition[2];
	cv::Mat alpha = DirectionCosines.at(0);
	cv::Mat beta = DirectionCosines.at(1);
	cv::Mat gamma = DirectionCosines.at(2);

	// Angle of Incidence
	cv::Mat cosIncidence = -(alpha*a+beta*b+gamma*c);
	cv::Mat Sign = signum(cosIncidence);
	// Snell's law
	cv::Mat lhsSnell;
	cv::sqrt(1.0-cosIncidence.mul(cosIncidence),lhsSnell);
	cv::Mat sinRefracted =  n_0/n_1*lhsSnell;
	cv::Mat cosRefracted;
	cv::sqrt(1.0-sinRefracted.mul(sinRefracted),cosRefracted);
	// Refracted Ray
	cv::Mat T_x = n_0/n_1*alpha+(n_0/n_1*cosIncidence-Sign.mul(cosRefracted))*a;
	cv::Mat T_y = n_0/n_1*beta+(n_0/n_1*cosIncidence-Sign.mul(cosRefracted))*b;
	cv::Mat T_z = n_0/n_1*gamma+(n_0/n_1*cosIncidence-Sign.mul(cosRefracted))*c;

	std::vector<cv::Mat> ReturnVector;
	ReturnVector.push_back(T_x);
	ReturnVector.push_back(T_y);
	ReturnVector.push_back(T_z);
	return ReturnVector;
}
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> SnellsLaw(const std::vector<cv::Mat> &DirectionCosines, const std::vector<double> &PlaneDefinition, const double &n_0, const cv::Mat &n_1)
{
  // n_1 is a Matrix (entry to fluid domain)
  cv::Mat n_0_mat(n_1.size(), CV_64FC1, n_0);

	// Normal vector of Plane
	double a = PlaneDefinition[0];
	double b = PlaneDefinition[1];
	double c = PlaneDefinition[2];
	cv::Mat alpha = DirectionCosines.at(0);
	cv::Mat beta = DirectionCosines.at(1);
	cv::Mat gamma = DirectionCosines.at(2);

	// Angle of Incidence
	cv::Mat cosIncidence = -(alpha*a+beta*b+gamma*c);
	cv::Mat Sign = signum(cosIncidence);
	cv::Mat Ratio = n_0_mat.mul(1.0/n_1);
	// Snell's law
	cv::Mat lhsSnell;
	cv::sqrt(1.0-cosIncidence.mul(cosIncidence),lhsSnell);
  cv::Mat sinRefracted = Ratio.mul(lhsSnell);
	cv::Mat cosRefracted;
	cv::sqrt(1.0-sinRefracted.mul(sinRefracted),cosRefracted);
	// Refracted Ray
	cv::Mat T_x = Ratio.mul(alpha)+(Ratio.mul(cosIncidence)-Sign.mul(cosRefracted))*a;
	cv::Mat T_y = Ratio.mul(beta)+(Ratio.mul(cosIncidence)-Sign.mul(cosRefracted))*b;
	cv::Mat T_z = Ratio.mul(gamma)+(Ratio.mul(cosIncidence)-Sign.mul(cosRefracted))*c;

	std::vector<cv::Mat> ReturnVector;
	ReturnVector.push_back(T_x);
	ReturnVector.push_back(T_y);
	ReturnVector.push_back(T_z);
	return ReturnVector;
}
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> SnellsLaw(const std::vector<cv::Mat> &DirectionCosines, const std::vector<double> &PlaneDefinition, const cv::Mat &n_0, const double &n_1)
{
  // n_0 is a Matrix (exit from fluid domain)
  cv::Mat n_1_mat(n_0.size(), CV_64FC1, n_1);

	// Normal vector of Plane
	double a = PlaneDefinition[0];
	double b = PlaneDefinition[1];
	double c = PlaneDefinition[2];
	cv::Mat alpha = DirectionCosines.at(0);
	cv::Mat beta = DirectionCosines.at(1);
	cv::Mat gamma = DirectionCosines.at(2);

	// Angle of Incidence
	cv::Mat cosIncidence = -(alpha*a+beta*b+gamma*c);
	cv::Mat Sign = signum(cosIncidence);
	cv::Mat Ratio = n_0.mul(1.0/n_1_mat);
	// Snell's law
	cv::Mat lhsSnell;
	cv::sqrt(1.0-cosIncidence.mul(cosIncidence),lhsSnell);
  cv::Mat sinRefracted = Ratio.mul(lhsSnell);
	cv::Mat cosRefracted;
	cv::sqrt(1.0-sinRefracted.mul(sinRefracted),cosRefracted);
	// Refracted Ray
	cv::Mat T_x = Ratio.mul(alpha)+(Ratio.mul(cosIncidence)-Sign.mul(cosRefracted))*a;
	cv::Mat T_y = Ratio.mul(beta)+(Ratio.mul(cosIncidence)-Sign.mul(cosRefracted))*b;
	cv::Mat T_z = Ratio.mul(gamma)+(Ratio.mul(cosIncidence)-Sign.mul(cosRefracted))*c;

	std::vector<cv::Mat> ReturnVector;
	ReturnVector.push_back(T_x);
	ReturnVector.push_back(T_y);
	ReturnVector.push_back(T_z);
	return ReturnVector;
}
/*--------------------------------------------------------------------------*/
extern std::vector<cv::Mat> SnellsLaw(const std::vector<cv::Mat> &DirectionCosines, const std::vector<double> &PlaneDefinition, const cv::Mat &n_0, const cv::Mat &n_1)
{
	// Normal vector of Plane
	double a = PlaneDefinition[0];
	double b = PlaneDefinition[1];
	double c = PlaneDefinition[2];
	cv::Mat alpha = DirectionCosines.at(0);
	cv::Mat beta = DirectionCosines.at(1);
	cv::Mat gamma = DirectionCosines.at(2);

	// Angle of Incidence
	cv::Mat cosIncidence = -(alpha*a+beta*b+gamma*c);
	cv::Mat Sign = signum(cosIncidence);
	cv::Mat Ratio = n_0.mul(1.0/n_1);
	// Snell's law
	cv::Mat lhsSnell;
	cv::sqrt(1.0-cosIncidence.mul(cosIncidence),lhsSnell);
	cv::Mat sinRefracted = Ratio.mul(lhsSnell);
	cv::Mat cosRefracted;
	cv::sqrt(1.0-sinRefracted.mul(sinRefracted),cosRefracted);
	// Refracted Ray
	cv::Mat T_x = Ratio.mul(alpha)+(Ratio.mul(cosIncidence)-Sign.mul(cosRefracted))*a;
	cv::Mat T_y = Ratio.mul(beta)+(Ratio.mul(cosIncidence)-Sign.mul(cosRefracted))*b;
	cv::Mat T_z = Ratio.mul(gamma)+(Ratio.mul(cosIncidence)-Sign.mul(cosRefracted))*c;

	std::vector<cv::Mat> ReturnVector;
	ReturnVector.push_back(T_x);
	ReturnVector.push_back(T_y);
	ReturnVector.push_back(T_z);
	return ReturnVector;
}
/*--------------------------------------------------------------------------*/
