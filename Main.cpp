#include <algorithm>
#include <iostream>
#include <memory>
#include <random>
#include <thread>
#include <Rainflow.h>
#include <RamerDouglasPeucker.h>


#define Py_LIMITED_API
#ifdef _DEBUG
	#undef _DEBUG
	#include <python.h>
	#define _DEBUG
#else
    #include <Python.h>
#endif


#if defined(_MSC_VER) || defined(WIN64) || defined(_WIN64) || defined(__WIN64__) || defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
	#define EXPORT __declspec(dllexport)
	#define IMPORT __declspec(dllimport)

	#define WIN32_LEAN_AND_MEAN
	#include <windows.h>

	static size_t getAvailablePhysicalMemory()
	{
		MEMORYSTATUSEX statex;
		statex.dwLength = sizeof(statex);
		GlobalMemoryStatusEx(&statex);
		size_t availableMBs = statex.ullAvailPhys / (1024 * 1024);
		return availableMBs;
	}

#else
	#define EXPORT __attribute__((visibility("default")))
	#define IMPORT __attribute__((visibility("default")))
	#define __stdcall

	#include <unistd.h>

	static size_t getAvailablePhysicalMemory()
	{
		long pages = sysconf(_SC_AVPHYS_PAGES);
		long page_size = sysconf(_SC_PAGE_SIZE);
		size_t availableMBs = pages * page_size / (1024 * 1024);
		return availableMBs;
	}

#endif


static std::vector<double> toVector(PyObject* values)
{
	/* Convert PyList to std::vector */
	size_t psize = PyList_Size(values);
	std::vector<double> vector;
	vector.reserve(psize);

	for (size_t i = 0; i < psize; ++i)
        vector.emplace_back(PyFloat_AsDouble(PyList_GetItem(values, i)));

    return vector;
}


static void randomiseOrder(PyObject* histories, std::vector<size_t>* ihistories, bool randomise)
{
	/* Get sizes */
	size_t nhistories = PyList_Size(histories);
	size_t hszmax = 0;
	size_t ihsize = 0;

	/* Iterate stress histories and counts */
	for (size_t i = 0; i < nhistories; ++i)
	{
		PyObject* item = PyList_GetItem(histories, i);
		size_t icount = size_t(PyLong_AsUnsignedLongLong(PyList_GetItem(item, 0)));

		/* Get size of longest stress history */
		hszmax = std::max<size_t>(hszmax, icount);

		/* Get total length of ihistories vector */
		ihsize += icount;
	}

	/* Indices of stress histories, reserve size */
	ihistories->reserve(ihsize);

	/* Define fixed order by iterating over stress histories first */
	for (size_t h = 0; h < hszmax; ++h)
		for (size_t i = 0; i < nhistories; ++i)
		{
			PyObject* item = PyList_GetItem(histories, i);
			size_t icount = size_t(PyLong_AsUnsignedLongLong(PyList_GetItem(item, 0)));

			if (h < icount)
				ihistories->emplace_back(i);
		}

	/* Randomise order/indices of stress histories */
	if (randomise)
	{
		std::random_device rd;
		std::default_random_engine engine(rd());
		std::shuffle(ihistories->begin(), ihistories->end(), engine);
	}
}

extern "C"
{
	EXPORT PyObject* __stdcall PyGetPeaks(PyObject* values)
	{
		std::vector<double> peaks = toVector(values);
		peaks = Rainflow::getPeaks(peaks);

		PyObject* output = PyList_New(peaks.size());
		for (int i = 0; i < peaks.size(); ++i)
			PyList_SetItem(output, i, PyFloat_FromDouble(peaks[i]));

		return output;
	}


	EXPORT PyObject* __stdcall PyGetPeakLocations(PyObject* xpositions, PyObject* yvalues)
	{
		/* Predefine output */
		PyObject* output = PyTuple_New(2);

		/* Get sizes */
		size_t psize = PyList_Size(xpositions);
		if (psize != size_t(PyList_Size(yvalues)) || psize <= 2)
		{
			if (psize != size_t(PyList_Size(yvalues)))
				std::cerr << "Positions and values have different lengths!" << std::endl;

			PyTuple_SetItem(output, 0, xpositions);
			PyTuple_SetItem(output, 1, yvalues);
			return output;
		}

		/* Convert PyList to std::vector */
		std::vector<double> positions = toVector(xpositions);
		std::vector<double> values = toVector(yvalues);

		/* Get peak locations */
		std::pair<std::vector<double>, std::vector<double>> points = Rainflow::getPeakLocations(positions, values);

		/* Output */
		psize = points.first.size();
		PyObject* pylocs = PyList_New(psize);
		PyObject* pypeaks = PyList_New(psize);

		for (size_t i = 0; i < psize; ++i)
		{
			PyList_SetItem(pylocs, i, PyFloat_FromDouble(points.first[i]));
			PyList_SetItem(pypeaks, i, PyFloat_FromDouble(points.second[i]));
		}

		PyTuple_SetItem(output, 0, pylocs);
		PyTuple_SetItem(output, 1, pypeaks);

		return output;
	}


	EXPORT PyObject* __stdcall PyRainflowCounting(PyObject* values, int algorithm, double tolerance, double cutoff)
	{
		/* Define tolerance for rounding and cutoff for output */
		std::vector<double> peaks = toVector(values);
		Rainflow Rf(peaks, tolerance, cutoff);

		/* Rainflow counting */
		switch (algorithm)
		{
			/* 3-points, non-periodic */
		case 1:
			Rf.rainflow3Points();
			break;

			/* 3-points, periodic */
		case 2:
			Rf.rotatePeaks();
			Rf.rainflow3Points();
			break;

			/* 4-points, periodic */
		default:
			Rf.rainflow4Points();
		}

		/* Output */
		PyObject* rdelta = PyList_New(Py_ssize_t(Rf.cycles.size()));
		PyObject* rmean = PyList_New(Py_ssize_t(Rf.cycles.size()));
		PyObject* rcount = PyList_New(Py_ssize_t(Rf.cycles.size()));

		size_t i = 0;
		for (std::map<std::pair<double, double>, double>::iterator itr = Rf.cycles.begin(); itr != Rf.cycles.end(); ++itr)
		{
			PyList_SetItem(rdelta, i, PyFloat_FromDouble(itr->first.first));
			PyList_SetItem(rmean, i, PyFloat_FromDouble(itr->first.second));
			PyList_SetItem(rcount, i, PyFloat_FromDouble(itr->second));
			++i;
		}
		PyObject* output = PyTuple_New(3);
		PyTuple_SetItem(output, 0, rdelta);
		PyTuple_SetItem(output, 1, rmean);
		PyTuple_SetItem(output, 2, rcount);

		return output;
	}


	EXPORT PyObject* __stdcall PyRainflowCountingRandomOrder(PyObject* histories, double tolerance, double cutoff, bool randomise, bool verify)
	{
		/* Randomise order of stress histories on a separate thread */
		size_t nhistories = PyList_Size(histories);
		std::vector<size_t> ihistories;
		std::thread t0(&randomiseOrder, histories, &ihistories, randomise);

		Rainflow RfAll(tolerance, cutoff);

		/* Perform rainflow counting on stress histories and save residue's for randomised order */
		std::vector< std::unique_ptr<std::vector<double>> > residues(nhistories);

		for (size_t i = 0; i < nhistories; ++i)
		{
			PyObject* item = PyList_GetItem(histories, i);
			size_t nmul = PyLong_AsUnsignedLongLong(PyList_GetItem(item, 0));
			PyObject* values = PyList_GetItem(item, 1);
			std::vector<double> peaks = toVector(values);

			Rainflow iRf(peaks, tolerance, cutoff);

			/* Iterate through turning points */
			iRf.rainflow4Points(false, nmul);

			/* Combine cycle counts */
			RfAll.addCycles(iRf.cycles);

			/* Store residue */
			residues[i] = std::make_unique<std::vector<double>>(iRf.rs);
		}

		/* Join thread before processing residues */
		t0.join();

		/* Define number of threads */
		unsigned nthreads = std::max<unsigned>(1, std::min<unsigned>(unsigned(ihistories.size()), std::thread::hardware_concurrency() / 2));
		if (verify)
			std::cout << "\t\tNumber of threads: " << nthreads << std::endl;

		/* Declare thread specific elements/structures */
		std::vector< std::unique_ptr<std::vector<double>> > tresidues(nthreads);
		std::vector< std::unique_ptr<std::map<std::pair<double, double>, double>> > tcycles(nthreads);

		/* Start threads */
		std::vector< std::unique_ptr<std::thread> > threads(nthreads);
		for (unsigned it = 0; it < nthreads; ++it)
			threads[it] = std::make_unique<std::thread>(std::thread(&Rainflow::rainflow4PointConcurrent, tolerance, cutoff, it, &ihistories, &residues, &tresidues, &tcycles));

		/* Wait for threads to finish */
		for (unsigned it = 0; it < nthreads; ++it)
			threads[it]->join();

		/* Free memory */
		residues.clear();
		threads.clear();
		if (!verify)
			ihistories.clear();

		/* Combine thread residues and cycle counts */
		std::vector<double> allresidues;
		allresidues.reserve(nthreads * 16);

		for (unsigned it = 0; it < nthreads; ++it)
		{
			allresidues.insert(allresidues.end(), tresidues[it]->begin(), tresidues[it]->end());
			RfAll.addCycles(*tcycles[it]);
		}

		/* Free memory */
		tresidues.clear();
		tcycles.clear();

		/* Process concatenated thread residues */
		RfAll.setPeaks(allresidues);
		RfAll.rainflow4Points();


		/* Verify results (SLOW) */
		if (verify)
		{
			Rainflow RfVerify(tolerance, cutoff);

			/* Get total size of concatenated stress histories */
			size_t sztotal = 0;
			for (size_t i = 0; i < nhistories; ++i)
			{
				PyObject* item = PyList_GetItem(histories, i);
				size_t nmul = PyLong_AsUnsignedLongLong(PyList_GetItem(item, 0));
				size_t length = PyList_Size(PyList_GetItem(item, 1));
				sztotal += (nmul * length);
			}

			/* Reserve memory */
			size_t availablememory = getAvailablePhysicalMemory();
			std::cout << "\t\tAvailable physical memory: " << availablememory << "MB" << std::endl;
			size_t nthreshold = std::min<size_t>(sztotal, std::max<size_t>(64, availablememory / 2 - 1024) * 131072);
			if (nthreshold / 131072 > 10)
				std::cout << "\t\tPhysical memory to be used: " << (2 * nthreshold / 131072) << "MB" << std::endl;
			else
				std::cout << "\t\tPhysical memory to be used: " << (2 * nthreshold / 128) << "KB" << std::endl;

			std::vector<double> v_allpeaks;
			v_allpeaks.reserve(nthreshold);
			std::vector<double> v_allresidues;
			v_allresidues.reserve(64);
			bool processed = false;

			/* Run on a single core, try to process entire signal at once if possible within memory */
			for (size_t i = 0; i < ihistories.size(); ++i)
			{
				/* Get stress history */
				PyObject* item = PyList_GetItem(histories, ihistories[i]);
				PyObject* values = PyList_GetItem(item, 1);
				std::vector<double> peaks = toVector(values);
				std::vector<double> vpeaks = Rainflow::getPeaks(peaks);

				/* Append stress history */
				v_allpeaks.insert(v_allpeaks.end(), vpeaks.begin(), vpeaks.end());

				/* Size of vector to be added in next iteration */
				size_t sznext = 0;
				if (i + 1 < ihistories.size())
					sznext = PyList_Size(PyList_GetItem(histories, ihistories[i + 1]));

				/* Process range of values and store residue */
				bool lastiteration = (i == ihistories.size() - 1);
				if ((lastiteration && processed) || (v_allpeaks.size() + sznext > nthreshold))
				{
					processed = true;
					RfVerify.setPeaks(v_allpeaks);
					RfVerify.rainflow4Points(false, 1);
					v_allresidues.insert(v_allresidues.end(), RfVerify.rs.begin(), RfVerify.rs.end());

					/* Process residue at last iteration */
					if (lastiteration)
					{
						RfVerify.setPeaks(v_allresidues);
						RfVerify.rainflow4Points();
					}
				}

				/* Process entire signal at once with 3-points, periodic algorithm */
				else if (lastiteration && !processed)
				{
					RfVerify.setPeaks(v_allpeaks);
					RfVerify.rotatePeaks();
					RfVerify.rainflow3Points();
					std::cout << "\t\tProcessed entire randomised stress history in one pass for verification" << std::endl;
				}
			}

			/* Compare cycles */
			if (RfVerify.cycles != RfAll.cycles)
				std::cerr << "\t\tCycle counts don't match!" << std::endl;
			else
				std::cout << "\t\tCycle counts successfully verified" << std::endl;

			/* Free memory */
			ihistories.clear();
		}

		/* Output */
		PyObject* rdelta = PyList_New(Py_ssize_t(RfAll.cycles.size()));
		PyObject* rmean = PyList_New(Py_ssize_t(RfAll.cycles.size()));
		PyObject* rcount = PyList_New(Py_ssize_t(RfAll.cycles.size()));

		size_t i = 0;
		for (std::map<std::pair<double, double>, double>::iterator itr = RfAll.cycles.begin(); itr != RfAll.cycles.end(); ++itr)
		{
			PyList_SetItem(rdelta, i, PyFloat_FromDouble(itr->first.first));
			PyList_SetItem(rmean, i, PyFloat_FromDouble(itr->first.second));
			PyList_SetItem(rcount, i, PyFloat_FromDouble(itr->second));
			++i;
		}
		PyObject* output = PyTuple_New(3);
		PyTuple_SetItem(output, 0, rdelta);
		PyTuple_SetItem(output, 1, rmean);
		PyTuple_SetItem(output, 2, rcount);

		return output;
	}


	EXPORT PyObject* __stdcall PyShift2ndLane(PyObject* histories_l1, PyObject* histories_l2, int step, int range, double tolerance)
	{
		/* Define tolerance for rounding  */
		Rainflow Rf(tolerance);

		/* Predefine output */
		PyObject* output = PyTuple_New(2);

		/* Get number of stress histories per lane/track */
		size_t n_his_lane1 = PyList_Size(histories_l1);
		size_t n_his_lane2 = PyList_Size(histories_l2);

		/* Check for equal lengths */
		size_t his_length = PyList_Size(PyList_GetItem(histories_l1, 0));

		for (size_t i = 0; i < n_his_lane1; ++i)
			if (his_length != PyList_Size(PyList_GetItem(histories_l1, i)))
			{
				std::cerr << "Different lengths of stress histories!" << std::endl;

				PyTuple_SetItem(output, 0, PyFloat_FromDouble(0.0));
				PyTuple_SetItem(output, 1, PyLong_FromLongLong(0));
				return output;
			}

		for (size_t i = 0; i < n_his_lane2; ++i)
			if (his_length != PyList_Size(PyList_GetItem(histories_l2, i)))
			{
				std::cerr << "Different lengths of stress histories!" << std::endl;
				return output;
			}

		/* Define range */
		long long i_step = std::max<long long>(1, step);
		long long i_range = (long long)std::floor(his_length / i_step) * i_step;
		if (range > 0)
			i_range = std::min<long long>(std::ceil(range / i_step) * i_step, i_range);
		size_t tot_length = i_range + his_length;

		/* Get averaged stress history for lane/track 1, rotated by one i_range */
		std::vector<double> shis_lane1(tot_length, 0.0);
		double sf1 = double(n_his_lane1);
		for (size_t i = 0; i < n_his_lane1; ++i)
		{
			PyObject* py_his = PyList_GetItem(histories_l1, i);

			for (size_t j = 0; j < his_length; ++j)
				shis_lane1[j + i_range] += PyFloat_AsDouble(PyList_GetItem(py_his, j)) / sf1;
		}

		/* Get averaged stress history for lane/track 2, not rotated */
		std::vector<double> shis_lane2(tot_length, 0.0);
		double sf2 = double(n_his_lane2);
		for (size_t i = 0; i < n_his_lane2; ++i)
		{
			PyObject* py_his = PyList_GetItem(histories_l2, i);

			for (size_t j = 0; j < his_length; ++j)
				shis_lane2[j] += PyFloat_AsDouble(PyList_GetItem(py_his, j)) / sf2;
		}

		/* Sum averaged stress histories */
		std::vector<double> shis_sum(tot_length, 0.0);
		for (size_t j = 0; j < shis_sum.size(); ++j)
			shis_sum[j] = shis_lane1[j] + shis_lane2[j];

		/* Start values for maximum stress delta and shift value */
		double delta_max = *std::max_element(shis_sum.begin(), shis_sum.end()) - *std::min_element(shis_sum.begin(), shis_sum.end());
		long long i_max = -1 * (long long)(i_range);

		/* Rotate each lane/track with fixed step size in opposite directions */
		for (long long i = 0; i < i_range; i += i_step)
		{
			/* Rotate/shift values of lane/track 1 by step size */
			std::rotate(shis_lane1.begin(), shis_lane1.begin() + i_step, shis_lane1.end());

			/* Recalculate sum of averaged stress histories */
			for (size_t j = 0; j < shis_sum.size(); ++j)
				shis_sum[j] = shis_lane1[j] + shis_lane2[j];

			/* Calculate maximum stress delta at current position */
			double delta = Rf.round(*std::max_element(shis_sum.begin(), shis_sum.end()) - *std::min_element(shis_sum.begin(), shis_sum.end()));
			long long i_cur = 2 * i + i_step - i_range;

			/* Save current position if stress delta is larger than current max */
			if (delta > delta_max)
			{
				delta_max = delta;
				i_max = i_cur;
			}
			/* Save current position if stress delta is equal to current max, but sift value is smaller */
			else if ((delta == delta_max) && (std::abs(i_cur) < std::abs(i_max)))
				i_max = i_cur;

			/* Rotate/shift values of lane/track 1 by step size, in opposite direction */
			std::rotate(shis_lane2.rbegin(), shis_lane2.rbegin() + i_step, shis_lane2.rend());

			/* Recalculate sum of averaged stress histories */
			for (size_t j = 0; j < shis_sum.size(); ++j)
				shis_sum[j] = shis_lane1[j] + shis_lane2[j];

			/* Calculate maximum stress delta at current position */
			delta = Rf.round(*std::max_element(shis_sum.begin(), shis_sum.end()) - *std::min_element(shis_sum.begin(), shis_sum.end()));
			i_cur += i_step;

			/* Save current position if stress delta is larger than current max */
			if (delta > delta_max)
			{
				delta_max = delta;
				i_max = i_cur;
			}
			/* Save current position if stress delta is equal to current max, but sift value is smaller */
			else if ((delta == delta_max) && (std::abs(i_cur) < std::abs(i_max)))
				i_max = i_cur;
		}

		/* Set return values */
		PyTuple_SetItem(output, 0, PyFloat_FromDouble(delta_max));
		PyTuple_SetItem(output, 1, PyLong_FromLongLong(i_max));

		return output;
	}


	EXPORT PyObject* __stdcall PySimplifyCurve(PyObject* xdata, PyObject* ydata, double epsilon)
	{
		/* Get curve data */
		size_t psize = std::min<size_t>(PyList_Size(xdata), PyList_Size(ydata));
		std::vector<Point> pointList;
		pointList.reserve(psize);

		for (size_t i = 0; i < psize; ++i)
			pointList.emplace_back(Point(PyFloat_AsDouble(PyList_GetItem(xdata, i)), PyFloat_AsDouble(PyList_GetItem(ydata, i))));

		/* Simplify curve data */
		std::vector<Point> simplified = RamerDouglasPeucker::SimplifyCurve(pointList, std::abs(epsilon));

		/* Define output */
		PyObject* output = PyTuple_New(2);

		PyObject* xout = PyList_New(simplified.size());
		PyObject* yout = PyList_New(simplified.size());

		for (int i = 0; i < simplified.size(); ++i)
		{
			PyList_SetItem(xout, i, PyFloat_FromDouble(simplified[i].first));
			PyList_SetItem(yout, i, PyFloat_FromDouble(simplified[i].second));
		}

		PyTuple_SetItem(output, 0, xout);
		PyTuple_SetItem(output, 1, yout);

		return output;
	}
}