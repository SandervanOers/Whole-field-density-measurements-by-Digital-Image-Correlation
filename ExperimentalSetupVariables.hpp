#ifndef H_ExperimentalSetupVariables
#define H_ExperimentalSetupVariables
struct ExperimentalSetupVariables {

	double focal_length;
	double Distance_From_Pixels_To_Meters;
	// constants
	double n_0 = 1.0003; // water
	double n_1 = 1.52; // glass
	// variables
	double n_cal; // index of refraction of calibration fluid
	double n_ref; // index of refraction of reference fluid
	double n;
	// Lengths Small Tank
	double L_c;
	double L_t;
	double L_g;
	double L_s;
	//std::vector<double> Lengths{L_c, L_g, L_t, L_s};

};
#endif
