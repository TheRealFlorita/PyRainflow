#pragma once
#include <vector>
#include <map>
#include <memory>
#include <SDKDDKVer.h>

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#define Py_LIMITED_API
#ifdef _DEBUG
	#undef _DEBUG
	#include <python.h>
	#define _DEBUG
#else
	#include <python.h>
#endif


class PyRainflow
{

private:

	double roundingTolerance, deltaCutOff;
	std::vector<double> peaks, rs;
	std::map<std::pair<double, double>, double> cycles;

public:

	PyRainflow(double tolerance = 1e-9, double cutoff = 1e-3);

	PyRainflow(PyObject* values, double tolerance = 1e-9, double cutoff = 1e-3);

	PyRainflow(std::vector<double>& values, double tolerance = 1e-9, double cutoff = 1e-3);

	PyRainflow(std::vector<double>& values, std::map<std::pair<double, double>, double>& cycles, double tolerance = 1e-9, double cutoff = 1e-3);

	~PyRainflow();

private:

	void setTolerance(double tolerance);

	void setCutOff(double cutoff);

	double round(double value);

	static size_t getAvailablePhysicalMemory();

	void addCycles(std::map<std::pair<double, double>, double>& counts);

	static bool isPeaksOnly(std::vector<double>& values);

	bool isPeaksOnly();

	static std::vector<double> getPeaks(std::vector<double>& values);

	static std::vector<double> getPeaks(PyObject* values);

	void setPeaks(std::vector<double>& values);

	void filterPeaks();

	void rotatePeaks();

	void rainflow3Points();

	void rainflow4Points(bool process_residue = true, size_t multiplier = 1);

	static void randomiseOrder(PyObject* histories, std::vector<size_t>* ihistories, bool randomise = true);

	static void rainflow4PointConcurrent(double tolerance, double cutoff, unsigned it,
		std::vector<size_t>* ihistories,
		std::vector< std::unique_ptr<std::vector<double>> >* residues,
		std::vector< std::unique_ptr<std::vector<double>> >* tresidues,
		std::vector< std::unique_ptr<std::map<std::pair<double, double>, double>> >* tcycles);


public:

	static PyObject* __stdcall PyGetPeaks(PyObject* values);

	static PyObject* __stdcall PyGetPeakLocations(PyObject* xpositions, PyObject* yvalues);

	static PyObject* __stdcall PyRainflowCounting(PyObject* values, int algorithm, double tolerance, double cutoff);

	static PyObject* __stdcall PyRainflowCountingRandomOrder(PyObject* histories, double tolerance, double cutoff, bool randomise, bool verify);

	static PyObject* __stdcall PyShift2ndLane(PyObject* histories_l1, PyObject* histories_l2, int step, int range, double tolerance);

};