#pragma once
#include <utility>
#include <vector>


typedef std::pair<double, double> Point;


class RamerDouglasPeucker
{

public:

	static double PerpendicularDistancePointToSegment(const Point& pt, const Point& segmentStart, const Point& segmentEnd);

	static std::vector<Point> SimplifyCurve(const std::vector<Point>& pointList, double epsilon);

	RamerDouglasPeucker() = delete;

};