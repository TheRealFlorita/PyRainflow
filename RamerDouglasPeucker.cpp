#include <cmath>
#include <iostream>
#include <RamerDouglasPeucker.h>


double RamerDouglasPeucker::PerpendicularDistancePointToSegment(const Point& pt, const Point& segmentStart, const Point& segmentEnd)
{
	/* Get delta's */
	double dx = segmentEnd.first - segmentStart.first;
	double dy = segmentEnd.second - segmentStart.second;

	/* Normalise */
    double a = std::pow(std::pow(dx, 2.0) + std::pow(dy, 2.0), 0.5);
	if (a > 0.0)
	{
		dx /= a;
		dy /= a;
	}

	double pvx = pt.first - segmentStart.first;
	double pvy = pt.second - segmentStart.second;

	/* Get dot product (project pv on normalised direction) */
	double pvdot = dx * pvx + dy * pvy;

	/* Scale line direction vector */
	double dsx = pvdot * dx;
	double dsy = pvdot * dy;

	/* Subtract this from pv */
	double ax = pvx - dsx;
	double ay = pvy - dsy;

    return std::pow(std::pow(ax, 2.0) + std::pow(ay, 2.0), 0.5);
}


std::vector<Point> RamerDouglasPeucker::SimplifyCurve(const std::vector<Point>& pointList, double epsilon)
{
	std::vector<Point> simplified;

	/* Return original list if simplification is not possible */
	if (pointList.size() < 2)
	{
		simplified.emplace_back(pointList.front());
		simplified.emplace_back(pointList.back());
		return simplified;
	}

	/* Find the point with the maximum distance from line between start and end */
	double dmax = 0.0;
	size_t index = 0;
	size_t end = pointList.size() - 1;
	for (size_t i = 1; i < end; i++)
	{
		double d = PerpendicularDistancePointToSegment(pointList[i], pointList[0], pointList[end]);
		if (d > dmax)
		{
			index = i;
			dmax = d;
		}
	}

	/* If max distance is greater than epsilon, recursively simplify */
	if (dmax > epsilon)
	{
		/* Recursive call */
		std::vector<Point> firstLine(pointList.begin(), pointList.begin() + index + 1);
		std::vector<Point> lastLine(pointList.begin() + index, pointList.end());
		std::vector<Point> recResults1 = SimplifyCurve(firstLine, epsilon);
		std::vector<Point> recResults2 = SimplifyCurve(lastLine, epsilon);

		/* Build the result list */
		simplified.assign(recResults1.begin(), recResults1.end() - 1);
		simplified.insert(simplified.end(), recResults2.begin(), recResults2.end());

		if (simplified.size() < 2)
			std::cerr << "RDP: Invalid length of simplified curve" << std::endl;
	}
	else
	{
		/* Just return start and end points */
		simplified.clear();
		simplified.emplace_back(pointList.front());
		simplified.emplace_back(pointList.back());
	}
	return simplified;
}
