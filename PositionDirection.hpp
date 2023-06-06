#ifndef Position_Direction
#define Position_Direction
struct PositionDirection {
	std::vector<cv::Mat> Position;
	std::vector<cv::Mat> Direction;
    PositionDirection(std::vector<cv::Mat> k, std::vector<cv::Mat> s) : Position(k), Direction(s) {}
};
#endif