#pragma once
#include <vector>
#include <map>
#include <memory>
#include <SDKDDKVer.h>

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#ifdef _DEBUG
	#undef _DEBUG
	#include <python.h>
	#define _DEBUG
#else
	#include <python.h>
#endif
#define Py_LIMITED_API


class PyRainflow
{

private:

	double roundingTolerance, deltaCutOff;

public:

	PyRainflow(double tolerance = 1e-9, double cutoff = 1e-3);

	~PyRainflow();

private:

	void setTolerance(double tolerance);

	void setCutOff(double cutoff);

	double round(double val);

	static size_t getAvailablePhysicalMemory();

	static bool isPeaksOnly(std::vector<double>& peaks);

	static std::vector<double> getPeaks(std::vector<double>& values);

	static std::vector<double> getPeaks(PyObject* values);

	static std::vector<double> rotatePeaks(std::vector<double>& peaks);

	void rainflow3Points(std::vector<double>& peaks, std::map<std::pair<double, double>, double>& cycles);

	std::vector<double> rainflow4Points(std::vector<double>& peaks, std::map<std::pair<double, double>, double>& cycles,
		bool residu = true, size_t multiplier = 1);

	static void randomiseOrder(PyObject* histories, std::vector<size_t>* ihistories, bool randomise = true);

	void rainflow4PointConcurrent(unsigned it, std::vector< std::unique_ptr<std::vector<double>> >* tresidues,
		std::vector<size_t>* ihistories, std::vector< std::unique_ptr<std::vector<double>> >* residues,
		std::vector< std::unique_ptr<std::map<std::pair<double, double>, double>> >* tcycles);

public:

	static PyObject* __stdcall PyGetPeaks(PyObject* values);

	static PyObject* __stdcall PyGetPeakLocations(PyObject* xpositions, PyObject* yvalues);

	static PyObject* __stdcall PyRainflowCounting(PyObject* values, int algorithm, double tolerance, double cutoff);

	static PyObject* __stdcall PyRainflowCountingRandomOrder(PyObject* histories, double tolerance, double cutoff, bool randomise, bool verify);

	static PyObject* __stdcall PyShift2ndLane(PyObject* histories_l1, PyObject* histories_l2, int step, int range, double tolerance);

};