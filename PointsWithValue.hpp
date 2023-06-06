#ifndef PointsWithValue
#define PointsWithValue
struct Points_With_Value {
    double C_value;
    cv::Point  Loc;
    Points_With_Value(double k, cv::Point s) : C_value(k), Loc(s) {}
};
#endif