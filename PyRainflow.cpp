#include <PyRainflow.h>
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>
#include <random>
#include <thread>
#include <windows.h>

#define Py_LIMITED_API

using namespace std;

double TOL = 1e-9;
double CUTOFF = 1e-9;


static size_t getAvailablePhysicalMemory()
{
	MEMORYSTATUSEX statex;
	statex.dwLength = sizeof(statex);
	GlobalMemoryStatusEx(&statex);
	size_t availableMBs = statex.ullAvailPhys / (1024 * 1024);
	return availableMBs;
}


static double mRound(double val)
{
	/* Round to multiples of tolerance */
	if (TOL > 0.0)
		return round(val / TOL) * TOL;
	return val;
}


static bool isPeaksOnly(vector<double>& peaks)
{
	/* Deal with vector size of two */
	if (peaks.size() == 2)
	{
		if (peaks[0] == peaks[1])
			return false;
	}
	/* Check for opposing signs of delta's between 3 points */
	else if (peaks.size() > 2)
	{
		for (size_t i = 1; i < peaks.size() - 1; ++i)
		{
			double delta1 = peaks[i] - peaks[i - 1];
			double delta2 = peaks[i + 1] - peaks[i];

			if ((delta1 * delta2) >= 0.0)
				return false;
		}
	}
	return true;
}


static vector<double> getPeaks(vector<double> &values)
{	
	/* Escape if vector contains only peaks */
	if (isPeaksOnly(values))
		return values;
		
	vector<double> peaks;

	/* Deal with small-sized vectors */
	if (values.size() <= 2)
	{
		peaks.insert(peaks.end(), values.begin(), values.end());
		if (values.size() == 2)
			if (values[0] == values[1])
				peaks.pop_back();
		return peaks;
	}

	/* Locate turning points and add to peaks vector */
	peaks.reserve(values.size());

	for (size_t i = 1; i < values.size() - 1; ++i)
	{
		double delta1 = values[i] - values[i - 1];
		double delta2 = values[i + 1] - values[i];

		if ((abs(delta1) > 0.0) && (delta1*delta2 <= 0.0))
		{
			if (peaks.size() == 0)
				peaks.emplace_back(values[i - 1]);
			else if (peaks.size() >= 2)
			{
				double point1 = peaks[peaks.size() - 2];
				double point2 = peaks.back();
				double point3 = values[i];

                if (point2 >= min(point1, point3) && point2 <= max(point1, point3))
					peaks.pop_back();
			}
			peaks.emplace_back(values[i]);
		}
	}

	/* Add last value to peaks vector if delta is non-zero */
	double delta = values.back() - peaks.back();
	if (abs(delta) > 0.0)
		peaks.emplace_back(values.back());

	return peaks;
}


static vector<double> getPeaks(PyObject* values)
{
	/* Convert PyList to std::vector */
	size_t psize = PyList_Size(values);
	vector<double> peaks;
	peaks.reserve(psize);

	for (size_t i = 0; i < psize; ++i)
		peaks.emplace_back(PyFloat_AsDouble(PyList_GetItem(values, i)));

	return getPeaks(peaks);
}


PyObject* __stdcall PyGetPeaks(PyObject* values)
{
	vector<double> peaks = getPeaks(values);

	PyObject* output = PyList_New(peaks.size());
	for (int i = 0; i < peaks.size(); ++i)
		PyList_SetItem(output, i, PyFloat_FromDouble(peaks[i]));

	return output;
}


static vector<double> rotatePeaks(vector<double> &peaks)
{
	/* Ensure the data only contains turning points */
	peaks = getPeaks(peaks);

	/* Rotate vector such that it starts with the absolute maximum value */
	if (peaks.size() > 2)
	{
		size_t idxmax = 0;
		double valmax = peaks[0];
		for (size_t i = 1; i < peaks.size() - 1; ++i)
		{
			if (abs(peaks[i]) > abs(valmax))
			{
				valmax = peaks[i];
				idxmax = i;
			}
		}

		rotate(peaks.begin(), peaks.begin() + idxmax, peaks.end());

		/* Append absolute maximum value to end to make signal/cycle count periodic */
		peaks.push_back(peaks[0]);
	}
	else if (peaks.size() == 2)
	{
		/* Rotate by one value if necessary */
		if (abs(peaks[0]) < abs(peaks[1]))
			rotate(peaks.begin(), peaks.begin() + 1, peaks.end());

		/* Append absolute maximum value */
		peaks.emplace_back(peaks[0]);
	}

    return getPeaks(peaks);
}


PyObject* __stdcall PyGetPeakLocations(PyObject* xpositions, PyObject* yvalues)
{
	/* Predefine output */
	PyObject* output = PyTuple_New(2);

	/* Get sizes */
	size_t psize = PyList_Size(xpositions);
	if (psize != size_t(PyList_Size(yvalues)) || psize <= 2)
	{
		if (psize != size_t(PyList_Size(yvalues)))
			cerr << "Positions and values have different lengths!" << endl;

		PyTuple_SetItem(output, 0, xpositions);
		PyTuple_SetItem(output, 1, yvalues);
		return output;
	}

	/* Convert PyList to std::vector */
	vector<double> positions;
	positions.reserve(psize);

	for (size_t i = 0; i < psize; ++i)
		positions.emplace_back(PyFloat_AsDouble(PyList_GetItem(xpositions, i)));

	/* Convert PyList to std::vector */
	vector<double> values;
	values.reserve(psize);

	for (size_t i = 0; i < psize; ++i)
		values.emplace_back(PyFloat_AsDouble(PyList_GetItem(yvalues, i)));

	/* Create vectors */
	vector<double> locations;
	locations.reserve(psize);

	vector<double> peaks;
	peaks.reserve(psize);

	/* Locate turning points and add to peaks vector */
	for (size_t i = 1; i < values.size() - 1; ++i)
	{
		double delta1 = values[i] - values[i - 1];
		double delta2 = values[i + 1] - values[i];

		if ((abs(delta1) > 0.0) && (delta1 * delta2 <= 0.0))
		{
			if (peaks.size() == 0)
			{
				locations.emplace_back(positions[i - 1]);
				peaks.emplace_back(values[i - 1]);
			}
			else if (peaks.size() >= 2)
			{
				double point1 = peaks[peaks.size() - 2];
				double point2 = peaks.back();
				double point3 = values[i];

				if (point2 >= min(point1, point3) && point2 <= max(point1, point3))
				{
					locations.pop_back();
					peaks.pop_back();
				}
			}

			locations.emplace_back(positions[i]);
			peaks.emplace_back(values[i]);
		}
	}

	/* Add last value to peaks vector if delta is non-zero */
	double delta = values.back() - peaks.back();
	if (abs(delta) > 0.0)
	{
		locations.emplace_back(positions.back());
		peaks.emplace_back(values.back());
	}

	/* Output */
	psize = locations.size();
	PyObject* locs = PyList_New(psize);
	PyObject* pks = PyList_New(psize);

	for (size_t i = 0; i < psize; ++i)
	{
		PyList_SetItem(locs, i, PyFloat_FromDouble(locations[i]));
		PyList_SetItem(pks, i, PyFloat_FromDouble(peaks[i]));
	}

	PyTuple_SetItem(output, 0, locs);
	PyTuple_SetItem(output, 1, pks);

	return output;
}


static void rainflow3Points(vector<double>& peaks, map<pair<double,double>, double>& cycles)
{
	/* Ensure the data only contains turning points */
	peaks = getPeaks(peaks);

	/* Create temporary vector */
	vector<double> rs;
	rs.reserve(64);

	/* Iterate through turning points */
	for (int i = 0; i < peaks.size(); ++i)
	{
		rs.emplace_back(peaks[i]);
		size_t sz = rs.size();

		while (sz >= 3 && abs(rs[sz - 3] - rs[sz - 2]) <= abs(rs[sz - 2] - rs[sz - 1]))
		{
			double delta = mRound(abs(rs[sz - 3] - rs[sz - 2]));
			double mean = mRound(0.5 * (rs[sz - 3] + rs[sz - 2]));

			/* Half cycle */
			if (sz == 3)
			{
				rs.erase(rs.begin());

				if (delta >= CUTOFF)
					cycles[pair<double,double>(delta, mean)] += 0.5;
			}

			/* Full cycle */
			else
			{
				rs.erase(rs.end() - 3, rs.end() - 1);

				if (delta >= CUTOFF)
					cycles[pair<double,double>(delta, mean)] += 1.0;
			}

			sz = rs.size();
		}
	}

	/* Residu */
	for (size_t i = 0; i < rs.size() - 1; ++i)
	{
		double delta = mRound(abs(rs[i] - rs[i + 1]));
		double mean = mRound(0.5 * (rs[i] + rs[i + 1]));

		if (delta >= CUTOFF)
			cycles[make_pair(delta, mean)] += 0.5;
	}

	peaks.clear();
}


static vector<double> rainflow4Points(vector<double>& peaks, map<pair<double,double>,double>& cycles, bool residu = true, size_t multiplier = 1)
{
	/* Ensure the data only contains turning points */
	peaks = getPeaks(peaks);

	/* Create temporary vector */
	vector<double> rs;
	rs.reserve(64);

	/* Iterate through turning points */
	for (size_t i = 0; i < peaks.size(); ++i)
	{
		rs.emplace_back(peaks[i]);
		size_t sz = rs.size();

        while ((sz >= 4) && (min(rs[sz - 3], rs[sz - 2]) >= min(rs[sz - 4], rs[sz - 1])) && (max(rs[sz - 3], rs[sz - 2]) <= max(rs[sz - 4], rs[sz - 1])))
		{
			double delta = mRound(abs(rs[sz - 3] - rs[sz - 2]));
			double mean = mRound(0.5 * (rs[sz - 3] + rs[sz - 2]));

			if (delta >= CUTOFF)
				cycles[pair<double,double>(delta, mean)] += 1.0 * multiplier;

			/* Full cycle */
			rs.erase(rs.end() - 3, rs.end() - 1);
			sz = rs.size();
		}
	}

	/* Process residu */
	if (residu)
	{
		rs.insert(rs.end(), rs.begin(), rs.end());
		rainflow4Points(rs, cycles, false);
		rs.clear();
	}

	peaks.clear();
	return rs;
}


PyObject* __stdcall PyRainflowCounting(PyObject* values, int algorithm, double tolerance, double cutoff)
{
	/* Define tolerance for rounding and cutoff for output */
    TOL = max<double>(tolerance, 1e-9);
	CUTOFF = max<double>(cutoff, 1e-9);

	/* Filter peaks and valleys */
	vector<double> peaks = getPeaks(values);

	/* Rainflow counting */
	map<pair<double,double>, double> cycles;

	switch (algorithm)
	{
	/* 3-points, non-periodic */
	case 1:
		rainflow3Points(peaks, cycles);
		break;

	/* 3-points, periodic */
	case 2:
		peaks = rotatePeaks(peaks);
		rainflow3Points(peaks, cycles);
		break;

	/* 4-points, periodic */
	default:
		rainflow4Points(peaks, cycles);
	}

	/* Output */
	PyObject* rdelta = PyList_New(Py_ssize_t(cycles.size()));
	PyObject* rmean = PyList_New(Py_ssize_t(cycles.size()));
	PyObject* rcount = PyList_New(Py_ssize_t(cycles.size()));

	size_t i = 0;
    for (map<pair<double,double>,double>::iterator itr = cycles.begin(); itr != cycles.end(); ++itr)
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


static void randomiseOrder(PyObject* histories, vector<size_t>* ihistories, bool randomise)
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
		hszmax = max(hszmax, icount);

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
		random_device rd;
		default_random_engine engine(rd());
		std::shuffle(ihistories->begin(), ihistories->end(), engine);
	}
}


void rainflow4PointConcurrent(unsigned it, vector< unique_ptr<vector<double>> >* tresidues, vector<size_t>* ihistories, vector< unique_ptr<vector<double>> >* residues, vector< unique_ptr<map<pair<double,double>,double>> >*  tcycles)
{
	/* Define range of loop */
	size_t nthreads = tresidues->size();
	size_t istep = max(1, size_t(round(double(ihistories->size()) / nthreads)));
	size_t istart = it * istep;
	size_t iend = min(istep * (it + 1), ihistories->size());

	/* Ensure last thread includes the tail of the vector */
	if (it == nthreads - 1)
		iend = ihistories->size();

	/* Only iterate when start value is smaller than end value */
	if (istart >= iend)
		return;

	/* Use unique pointers for dealing with memory allocation */
	unique_ptr<vector<double>> rs = make_unique<vector<double>>(vector<double>());
	rs->reserve(64);
	unique_ptr<map<pair<double, double>, double>> cycles = make_unique<map<pair<double, double>, double>>(map<pair<double, double>, double>());
	
	/* Determine size of peaks vector */
	size_t szsum = 0;
	size_t szmax = 0;
	for (size_t i = istart; i < iend; ++i)
	{
		szsum += residues->at(ihistories->at(i))->size();
		szmax = max(szmax, residues->at(ihistories->at(i))->size());
	}

	/* Try to limit size of peaks vector to 64 MB */
	size_t nthreshold = min(szsum, max(szmax, 64 * 131072));
	vector<double> peaks;
	peaks.reserve(nthreshold);

	/* Append turning points and process signal when size reaches nthreshold or end */
	for (size_t i = istart; i < iend; ++i)
	{
		/* Append next stress history/residue */
		peaks.insert(peaks.end(), residues->at(ihistories->at(i))->begin(), residues->at(ihistories->at(i))->end());

		/* Size of vector to be added in next iteration */
		size_t sznext = 0;
		if (i + 1 < ihistories->size())
			sznext = residues->at(ihistories->at(i + 1))->size();

		/* Process if size threshold will be exceeded in the next iteration or if this is the last iteration */
		if ((peaks.size() + sznext > nthreshold) || (i == iend - 1))
		{
			vector<double> irs = rainflow4Points(peaks, *cycles, false, 1);
			rs->insert(rs->end(), irs.begin(), irs.end());
		}
	}

	/* Store thread specific residue and cycle count */
	tresidues->at(it) = move(rs);
	tcycles->at(it) = move(cycles);
}


PyObject* __stdcall PyRainflowCountingRandomOrder(PyObject* histories, double tolerance, double cutoff, bool randomise, bool verify)
{
	/* Define tolerance for rounding and cutoff for output */
	TOL = max<double>(tolerance, 1e-9);
	CUTOFF = max<double>(cutoff, 1e-9);

	/* Randomise order of stress histories on a separate thread */
	size_t nhistories = PyList_Size(histories);
	vector<size_t> ihistories;
	thread t0(&randomiseOrder, histories, &ihistories, randomise);

	/* Perform rainflow counting on stress histories and save residu's for randomised order */
    map<pair<double,double>,double> cycles;
	vector< unique_ptr<vector<double>> > residues(nhistories);

	for (size_t i = 0; i < nhistories; ++i)
	{
		PyObject* item = PyList_GetItem(histories, i);
        size_t nmul = PyLong_AsUnsignedLongLong(PyList_GetItem(item, 0));
		vector<double> values = getPeaks(PyList_GetItem(item, 1));

		/* Iterate through turning points */
		unique_ptr<vector<double>> rs = make_unique<vector<double>>(rainflow4Points(values, cycles, false, nmul));

		/* Store residu */
		residues[i] = move(rs);
	}
	
	/* Join thread before processing residues */
	t0.join();

	/* Define number of threads */
    unsigned nthreads = max(1, min(unsigned(ihistories.size()), thread::hardware_concurrency()/2));
	if (verify)
		cout << "\t\tNumber of threads: " << nthreads << endl;

	/* Declare thread specific elements/structures */
	vector< unique_ptr<vector<double>> > tresidues(nthreads);
    vector< unique_ptr<map<pair<double,double>,double>> > tcycles(nthreads);

	/* Start threads */
	vector< unique_ptr<thread> > threads(nthreads);
	for (unsigned it = 0; it < nthreads; ++it)
		threads[it] = make_unique<thread>(thread(&rainflow4PointConcurrent, it, &tresidues, &ihistories, &residues, &tcycles));

	/* Wait for threads to finish */
	for (unsigned it = 0; it < nthreads; ++it)
		threads[it]->join();

	/* Free memory */
	residues.clear();
	threads.clear();
	if (!verify)
		ihistories.clear();

	/* Combine thread residues and cycle counts */
	vector<double> allresidues(nthreads*8);

	for (unsigned it = 0; it < nthreads; ++it)
	{
		allresidues.insert(allresidues.end(), tresidues[it]->begin(), tresidues[it]->end());

        for (map<pair<double,double>,double>::iterator itr = tcycles[it]->begin(); itr != tcycles[it]->end(); ++itr)
			cycles[itr->first] += itr->second;
	}

	/* Free memory */
	tresidues.clear();
	tcycles.clear();

	/* Process concatenated thread residues */
	rainflow4Points(allresidues, cycles);

	/* Verify results (SLOW) */
	if (verify)
	{
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
		cout << "\t\tAvailable physical memory: " << availablememory << "MB" << endl;
		size_t nthreshold = min(sztotal, max(64, availablememory/2 - 1024) * 131072);
		if (nthreshold / 131072 > 10)
			cout << "\t\tPhysical memory to be used: " << (2*nthreshold / 131072) << "MB" << endl;
		else
			cout << "\t\tPhysical memory to be used: " << (2*nthreshold / 128) << "KB" << endl;

		map<pair<double, double>, double> v_cycles;
		vector<double> v_allpeaks;
		v_allpeaks.reserve(nthreshold);
		vector<double> v_allresidues;
		v_allresidues.reserve(64);
		bool processed = false;

		/* Run on a single core, try to process entire signal at once if possible within memory */
		for (size_t i = 0; i < ihistories.size(); ++i)
		{
			/* Get stress history */
			PyObject* item = PyList_GetItem(histories, ihistories[i]);
			vector<double> peaks = getPeaks(PyList_GetItem(item, 1));
			
			/* Append stress history */
			v_allpeaks.insert(v_allpeaks.end(), peaks.begin(), peaks.end());

			/* Size of vector to be added in next iteration */
			size_t sznext = 0;
			if (i + 1 < ihistories.size())
				sznext = PyList_Size(PyList_GetItem(histories, ihistories[i+1]));

			/* Process range of values and store residue */
			bool lastiteration = (i == ihistories.size() - 1);
			if ((lastiteration && processed) || (v_allpeaks.size() + sznext > nthreshold))
			{
				processed = true;
				vector<double> rs = rainflow4Points(v_allpeaks, v_cycles, false, 1);
				v_allresidues.insert(v_allresidues.end(), rs.begin(), rs.end());

				/* Process residue at last iteration */
				if (lastiteration)
					rainflow4Points(v_allresidues, v_cycles);
			}

			/* Process entire signal at once with 3-points, periodic algorithm */
			else if (lastiteration && !processed)
			{
				v_allpeaks = rotatePeaks(v_allpeaks);
				rainflow3Points(v_allpeaks, v_cycles);
				cout << "\t\tProcessed entire randomised stress history in one pass for verification" << endl;
			}
		}

		/* Compare cycles */
		if (v_cycles != cycles)
			cerr << "\t\tCycle counts don't match!" << endl;
		else
			cout << "\t\tCycle counts successfully verified" << endl;
		
		/* Free memory */
		ihistories.clear();
	}

	/* Output */
	PyObject* rdelta = PyList_New(Py_ssize_t(cycles.size()));
	PyObject* rmean = PyList_New(Py_ssize_t(cycles.size()));
	PyObject* rcount = PyList_New(Py_ssize_t(cycles.size()));

	size_t i = 0;
	for (map<pair<double, double>, double>::iterator itr = cycles.begin(); itr != cycles.end(); ++itr)
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


PyObject* __stdcall PyShift2ndLane(PyObject* histories_l1, PyObject* histories_l2, int step, int range, double tolerance)
{
	/* Define tolerance for rounding  */
	TOL = max<double>(tolerance, 1e-9);

	/* Predefine output */
	PyObject* output = PyTuple_New(2);

	/* Get number of stress histories per lane */
	size_t n_his_lane1 = PyList_Size(histories_l1);
	size_t n_his_lane2 = PyList_Size(histories_l2);

	/* Check for equal lengths */
	size_t his_length = PyList_Size(PyList_GetItem(histories_l1, 0));

	for (size_t i = 0; i < n_his_lane1; ++i)
		if (his_length != PyList_Size(PyList_GetItem(histories_l1, i)))
		{
			cerr << "Different lengths of stress histories!" << endl;

			PyTuple_SetItem(output, 0, PyFloat_FromDouble(0.0));
			PyTuple_SetItem(output, 1, PyLong_FromLongLong(0));
			return output;
		}

	for (size_t i = 0; i < n_his_lane2; ++i)
		if (his_length != PyList_Size(PyList_GetItem(histories_l2, i)))
		{
			cerr << "Different lengths of stress histories!" << endl;
			return output;
		}

	/* Define range */
	size_t his_range = his_length;
	if (range > 0)
		size_t his_range = min(range, his_length);
	size_t his_total = his_range + his_length;

	/* Get averaged stress history for lane 1, rotated by one length */
	vector<double> shis_lane1(his_total, 0.0);
	for (size_t i = 0; i < n_his_lane1; ++i)
	{
		double sf1 = double(n_his_lane1);
		PyObject* py_his = PyList_GetItem(histories_l1, i);

		for (size_t j = 0; j < his_length; ++j)
			shis_lane1[j + his_range] += PyFloat_AsDouble(PyList_GetItem(py_his, j)) / sf1;
	}

	/* Get averaged stress history for lane 2 */
	vector<double> shis_lane2(his_total, 0.0);
	for (size_t i = 0; i < n_his_lane2; ++i)
	{
		double sf2 = double(n_his_lane2);
		PyObject* py_his = PyList_GetItem(histories_l2, i);

		for (size_t j = 0; j < his_length; ++j)
			shis_lane2[j] += PyFloat_AsDouble(PyList_GetItem(py_his, j)) / sf2;
	}

	/* Sum averaged stress histories */
	vector<double> shis_sum(his_total, 0.0);
	for (size_t j = 0; j < shis_sum.size(); ++j)
		shis_sum[j] = shis_lane1[j] + shis_lane2[j];

	/* Start values for maximum stress delta and shift value */
	double delta_max = *max_element(shis_sum.begin(), shis_sum.end()) - *min_element(shis_sum.begin(), shis_sum.end());
    long long i_max = -1 * (long long)(his_range);

	/* Rotate each lane with fixed stepsize in opposite directions */
    size_t i_step = size_t(max(1, step));

	for (size_t i = 0; i < his_range; i += i_step)
	{
		/* Rotate/shift values of lane 1 by stepsize */
		rotate(shis_lane1.begin(), shis_lane1.begin() + i_step, shis_lane1.end());

		/* Recalculate sum of averaged stress histories */
		for (size_t j = 0; j < shis_sum.size(); ++j)
			shis_sum[j] = shis_lane1[j] + shis_lane2[j];

		/* Calculate maximum stress delta at current position */
		double delta = mRound(*max_element(shis_sum.begin(), shis_sum.end()) - *min_element(shis_sum.begin(), shis_sum.end()));
        long long i_cur = long long(i * 2 + i_step) - (long long)(his_range);

		/* Save current position if stress delta is larger than current max */
		if (delta > delta_max)
		{
			delta_max = delta;
			i_max = i_cur;
		}
		/* Save current position if stress delta is equal to current max, but sift value is smaller */
		else if ((delta == delta_max) && (abs(i_cur) < abs(i_max)))
			i_max = i_cur;

		/* Rotate/shift values of lane 1 by stepsize, in opposite direction */
		rotate(shis_lane2.rbegin(), shis_lane2.rbegin() + i_step, shis_lane2.rend());

		/* Recalculate sum of averaged stress histories */
		for (size_t j = 0; j < shis_sum.size(); ++j)
			shis_sum[j] = shis_lane1[j] + shis_lane2[j];

		/* Calculate maximum stress delta at current position */
		delta = mRound(*max_element(shis_sum.begin(), shis_sum.end()) - *min_element(shis_sum.begin(), shis_sum.end()));
		i_cur = (long long)(i + i_step) * 2 - (long long)(his_range);

		/* Save current position if stress delta is larger than current max */
		if (delta > delta_max)
		{
			delta_max = delta;
			i_max = i_cur;
		}
		/* Save current position if stress delta is equal to current max, but sift value is smaller */
		else if ((delta == delta_max) && (abs(i_cur) < abs(i_max)))
			i_max = i_cur;
	}

	/* Set return values */
	PyTuple_SetItem(output, 0, PyFloat_FromDouble(delta_max));
	PyTuple_SetItem(output, 1, PyLong_FromLongLong(i_max));

	return output;
}
