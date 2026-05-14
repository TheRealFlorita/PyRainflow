#pragma once
#include <map>
#include <memory>
#include <vector>


class Rainflow
{

public:

	double roundingTolerance, deltaCutOff;

	std::vector<double> peaks, rs;

	std::map<std::pair<double, double>, double> cycles;

	Rainflow(double tolerance = 1e-9, double cutoff = 1e-3);

	Rainflow(std::vector<double>& values, double tolerance = 1e-9, double cutoff = 1e-3);

	Rainflow(std::vector<double>& values, std::map<std::pair<double, double>, double>& cycles, double tolerance = 1e-9, double cutoff = 1e-3);

	~Rainflow();

	void setTolerance(double tolerance);

	void setCutOff(double cutoff);

	double round(double value);

	void addCycles(std::map<std::pair<double, double>, double>& counts);

	static bool isPeaksOnly(std::vector<double>& values);

	bool isPeaksOnly();

	static std::vector<double> getPeaks(std::vector<double>& values);

	void setPeaks(std::vector<double>& values);

	void filterPeaks();

	void rotatePeaks();

	static std::pair<std::vector<double>, std::vector<double>> getPeakLocations(std::vector<double> positions, std::vector<double> values);

	void rainflow3Points();

	void rainflow4Points(bool process_residue = true, size_t multiplier = 1);

	static void rainflow4PointConcurrent(double tolerance, double cutoff, unsigned it,
		std::vector<size_t>* ihistories,
		std::vector< std::unique_ptr<std::vector<double>> >* residues,
		std::vector< std::unique_ptr<std::vector<double>> >* tresidues,
		std::vector< std::unique_ptr<std::map<std::pair<double, double>, double>> >* tcycles);

};