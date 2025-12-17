#include <PyRainflow.h>
#include <iostream>
#include <algorithm>
#include <random>
#include <thread>


PyRainflow::PyRainflow(double tolerance, double cutoff)
{
	setTolerance(tolerance);
	setCutOff(cutoff);
};


PyRainflow::PyRainflow(PyObject* values, double tolerance, double cutoff)
{
	setTolerance(tolerance);
	setCutOff(cutoff);

	peaks = getPeaks(values);
};


PyRainflow::PyRainflow(std::vector<double>& values, double tolerance, double cutoff)
{
	setTolerance(tolerance);
	setCutOff(cutoff);

	peaks = getPeaks(values);
};


PyRainflow::PyRainflow(std::vector<double>& values, std::map<std::pair<double, double>, double>& counts, double tolerance, double cutoff)
{
	setTolerance(tolerance);
	setCutOff(cutoff);

	peaks = getPeaks(values);
	cycles = counts;
};


PyRainflow::~PyRainflow() {};


void PyRainflow::setTolerance(double tolerance)
{
	roundingTolerance = std::max<double>(1e-16, tolerance);
};


void PyRainflow::setCutOff(double cutoff)
{
	deltaCutOff = std::max<double>(1e-16, cutoff);
};


double PyRainflow::round(double value)
{
	return std::round(value / roundingTolerance) * roundingTolerance;
}


size_t PyRainflow::getAvailablePhysicalMemory()
{
	MEMORYSTATUSEX statex;
	statex.dwLength = sizeof(statex);
	GlobalMemoryStatusEx(&statex);
	size_t availableMBs = statex.ullAvailPhys / (1024 * 1024);
	return availableMBs;
}


void PyRainflow::addCycles(std::map<std::pair<double, double>, double>& counts)
{
	for (std::map<std::pair<double, double>, double>::iterator itr = counts.begin(); itr != counts.end(); ++itr)
		cycles[itr->first] += itr->second;

	return;
}


bool PyRainflow::isPeaksOnly(std::vector<double>& values)
{
	/* Deal with vector size of two */
	if (values.size() == 2)
	{
		if (values[0] == values[1])
			return false;
	}
	/* Check for opposing signs of delta's between 3 points */
	else if (values.size() > 2)
	{
		for (size_t i = 1; i < values.size() - 1; ++i)
		{
			double delta1 = values[i] - values[i - 1];
			double delta2 = values[i + 1] - values[i];

			if ((delta1 * delta2) >= 0.0)
				return false;
		}
	}
	return true;
}


bool PyRainflow::isPeaksOnly()
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


std::vector<double> PyRainflow::getPeaks(std::vector<double>& values)
{
	/* Escape if vector contains only peaks */
	if (isPeaksOnly(values))
		return values;

	std::vector<double> rpeaks;

	/* Deal with small-sized vectors */
	if (values.size() <= 2)
	{
		rpeaks.insert(rpeaks.end(), values.begin(), values.end());
		if (values.size() == 2)
			if (values[0] == values[1])
				rpeaks.pop_back();
		return rpeaks;
	}

	/* Locate turning points and add to peaks vector */
	rpeaks.reserve(values.size());

	for (size_t i = 1; i < values.size() - 1; ++i)
	{
		double delta1 = values[i] - values[i - 1];
		double delta2 = values[i + 1] - values[i];

		if ((std::abs(delta1) > 0.0) && (delta1 * delta2 <= 0.0))
		{
			if (rpeaks.size() == 0)
				rpeaks.emplace_back(values[i - 1]);
			else if (rpeaks.size() >= 2)
			{
				double point1 = rpeaks[rpeaks.size() - 2];
				double point2 = rpeaks.back();
				double point3 = values[i];

				if (point2 >= std::min<double>(point1, point3) && point2 <= std::max<double>(point1, point3))
					rpeaks.pop_back();
			}
			rpeaks.emplace_back(values[i]);
		}
	}

	/* Add last value to peaks vector if delta is non-zero */
	double delta = values.back() - rpeaks.back();
	if (std::abs(delta) > 0.0)
		rpeaks.emplace_back(values.back());

	return rpeaks;
}


std::vector<double> PyRainflow::getPeaks(PyObject* values)
{
	/* Convert PyList to std::vector */
	size_t psize = PyList_Size(values);
	std::vector<double> rpeaks;
	rpeaks.reserve(psize);

	for (size_t i = 0; i < psize; ++i)
		rpeaks.emplace_back(PyFloat_AsDouble(PyList_GetItem(values, i)));

	return getPeaks(rpeaks);
}


void PyRainflow::setPeaks(std::vector<double>& values)
{
	peaks = getPeaks(values);
	return;
}


void PyRainflow::filterPeaks()
{
	/* Escape if vector contains only peaks */
	if (isPeaksOnly())
		return;

	/* Reset vector */
	std::vector<double> values = peaks;
	peaks.clear();

	/* Deal with small-sized vectors */
	if (values.size() <= 2)
	{
		peaks.insert(peaks.end(), values.begin(), values.end());
		if (values.size() == 2)
			if (values[0] == values[1])
				peaks.pop_back();
		return;
	}

	/* Locate turning points and add to peaks vector */
	for (size_t i = 1; i < values.size() - 1; ++i)
	{
		double delta1 = values[i] - values[i - 1];
		double delta2 = values[i + 1] - values[i];

		if ((std::abs(delta1) > 0.0) && (delta1 * delta2 <= 0.0))
		{
			if (peaks.size() == 0)
				peaks.emplace_back(values[i - 1]);
			else if (peaks.size() >= 2)
			{
				double point1 = peaks[peaks.size() - 2];
				double point2 = peaks.back();
				double point3 = values[i];

				if (point2 >= std::min<double>(point1, point3) && point2 <= std::max<double>(point1, point3))
					peaks.pop_back();
			}
			peaks.emplace_back(values[i]);
		}
	}

	/* Add last value to peaks vector if delta is non-zero */
	double delta = values.back() - peaks.back();
	if (std::abs(delta) > 0.0)
		peaks.emplace_back(values.back());

	return;
}


void PyRainflow::rotatePeaks()
{
	/* Ensure the data only contains turning points */
	filterPeaks();

	/* Rotate vector such that it starts with the absolute maximum value */
	if (peaks.size() > 2)
	{
		size_t idxmax = 0;
		double valmax = peaks[0];
		for (size_t i = 1; i < peaks.size() - 1; ++i)
		{
			if (std::abs(peaks[i]) > std::abs(valmax))
			{
				valmax = peaks[i];
				idxmax = i;
			}
		}

		std::rotate(peaks.begin(), peaks.begin() + idxmax, peaks.end());

		/* Append absolute maximum value to end to make signal/cycle count periodic */
		peaks.push_back(peaks[0]);
	}
	else if (peaks.size() == 2)
	{
		/* Rotate by one value if necessary */
		if (std::abs(peaks[0]) < std::abs(peaks[1]))
			std::rotate(peaks.begin(), peaks.begin() + 1, peaks.end());

		/* Append absolute maximum value */
		peaks.emplace_back(peaks[0]);
	}

	/* Ensure the data only contains turning points */
	filterPeaks();

	return;
}


void PyRainflow::rainflow3Points()
{
	/* Ensure the data only contains turning points */
	filterPeaks();

	/* Reset residue vector */
	rs.clear();
	rs.reserve(64);

	/* Iterate through turning points */
	for (int i = 0; i < peaks.size(); ++i)
	{
		rs.emplace_back(peaks[i]);
		size_t sz = rs.size();

		while (sz >= 3 && std::abs(rs[sz - 3] - rs[sz - 2]) <= std::abs(rs[sz - 2] - rs[sz - 1]))
		{
			double delta = round(std::abs(rs[sz - 3] - rs[sz - 2]));
			double mean = round(0.5 * (rs[sz - 3] + rs[sz - 2]));

			/* Half cycle */
			if (sz == 3)
			{
				rs.erase(rs.begin());

				if (delta >= deltaCutOff)
					cycles[std::pair<double, double>(delta, mean)] += 0.5;
			}

			/* Full cycle */
			else
			{
				rs.erase(rs.end() - 3, rs.end() - 1);

				if (delta >= deltaCutOff)
					cycles[std::pair<double, double>(delta, mean)] += 1.0;
			}

			sz = rs.size();
		}
	}

	/* residue */
	for (size_t i = 0; i < rs.size() - 1; ++i)
	{
		double delta = round(std::abs(rs[i] - rs[i + 1]));
		double mean = round(0.5 * (rs[i] + rs[i + 1]));

		if (delta >= deltaCutOff)
			cycles[std::make_pair(delta, mean)] += 0.5;
	}

	rs.clear();
	peaks.clear();
}


void PyRainflow::rainflow4Points(bool process_residue, size_t multiplier)
{
	/* Ensure the data only contains turning points */
	filterPeaks();

	/* Reset residue vector */
	rs.clear();
	rs.reserve(64);

	/* Perform procedure twice, first peaks, then residue */
	for (int c = 0; c < 1 + int(process_residue); ++c)
	{
		/* Process residue */
		if (c > 0)
		{
			peaks = rs;
			peaks.insert(peaks.end(), peaks.begin(), peaks.end());
			rs.clear();

			filterPeaks();
		}

		/* Iterate through turning points */
		for (size_t i = 0; i < peaks.size(); ++i)
		{
			rs.emplace_back(peaks[i]);
			size_t sz = rs.size();

			while ((sz >= 4) && (std::min<double>(rs[sz - 3], rs[sz - 2]) >= std::min<double>(rs[sz - 4], rs[sz - 1])) && (std::max<double>(rs[sz - 3], rs[sz - 2]) <= std::max<double>(rs[sz - 4], rs[sz - 1])))
			{
				double delta = round(std::abs(rs[sz - 3] - rs[sz - 2]));
				double mean = round(0.5 * (rs[sz - 3] + rs[sz - 2]));

				if (delta >= deltaCutOff)
					cycles[std::pair<double, double>(delta, mean)] += 1.0 * multiplier;

				/* Full cycle */
				rs.erase(rs.end() - 3, rs.end() - 1);
				sz = rs.size();
			}
		}
	}

	/* Clear residue only if processed */
	if (process_residue)
		rs.clear();

	peaks.clear();
	return;
}


void PyRainflow::randomiseOrder(PyObject* histories, std::vector<size_t>* ihistories, bool randomise)
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


void PyRainflow::rainflow4PointConcurrent(double tolerance, double cutoff, unsigned it,
	std::vector<size_t>* ihistories,
	std::vector< std::unique_ptr<std::vector<double>> >* residues,
	std::vector< std::unique_ptr<std::vector<double>> >* tresidues,
	std::vector< std::unique_ptr<std::map<std::pair<double, double>, double>> >* tcycles)
{
	/* Define range of loop */
	size_t nthreads = tresidues->size();
	size_t istep = std::max<size_t>(1, size_t(std::round(double(ihistories->size()) / nthreads)));
	size_t istart = it * istep;
	size_t iend = std::min<size_t>(istep * (it + 1), ihistories->size());

	/* Ensure last thread includes the tail of the vector */
	if (it == nthreads - 1)
		iend = ihistories->size();

	/* Only iterate when start value is smaller than end value */
	if (istart >= iend)
		return;

	/* Determine size of peaks vector */
	size_t szsum = 0;
	size_t szmax = 0;
	for (size_t i = istart; i < iend; ++i)
	{
		szsum += residues->at(ihistories->at(i))->size();
		szmax = std::max<size_t>(szmax, residues->at(ihistories->at(i))->size());
	}

	/* Try to limit size of peaks vector to 64 MB */
	size_t nthreshold = std::min<size_t>(szsum, std::max<size_t>(szmax, 64 * 131072));
	std::vector<double> tpeaks;
	tpeaks.reserve(nthreshold);
	std::vector<double> trs;
	trs.reserve(64);

	PyRainflow iRf(tolerance, cutoff);

	/* Append turning points and process signal when size reaches nthreshold or end */
	for (size_t i = istart; i < iend; ++i)
	{
		/* Append next stress history/residue */
		tpeaks.insert(tpeaks.end(), residues->at(ihistories->at(i))->begin(), residues->at(ihistories->at(i))->end());

		/* Size of vector to be added in next iteration */
		size_t sznext = 0;
		if (i + 1 < ihistories->size())
			sznext = residues->at(ihistories->at(i + 1))->size();

		/* Process if size threshold will be exceeded in the next iteration or if this is the last iteration */
		if ((tpeaks.size() + sznext > nthreshold) || (i == iend - 1))
		{
			iRf.setPeaks(tpeaks);
			iRf.rainflow4Points(false, 1);
			trs.insert(trs.end(), iRf.rs.begin(), iRf.rs.end());
		}
	}

	/* Store thread specific residue and cycle count */
	tresidues->at(it) = std::make_unique<std::vector<double>>(trs);
	tcycles->at(it) = std::make_unique<std::map<std::pair<double, double>, double>>(iRf.cycles);
}


PyObject* __stdcall PyRainflow::PyGetPeaks(PyObject* values)
{
	std::vector<double> pypeaks = getPeaks(values);

	PyObject* output = PyList_New(pypeaks.size());
	for (int i = 0; i < pypeaks.size(); ++i)
		PyList_SetItem(output, i, PyFloat_FromDouble(pypeaks[i]));

	return output;
}


PyObject* __stdcall PyRainflow::PyGetPeakLocations(PyObject* xpositions, PyObject* yvalues)
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
	std::vector<double> positions;
	positions.reserve(psize);

	for (size_t i = 0; i < psize; ++i)
		positions.emplace_back(PyFloat_AsDouble(PyList_GetItem(xpositions, i)));

	/* Convert PyList to std::vector */
	std::vector<double> values;
	values.reserve(psize);

	for (size_t i = 0; i < psize; ++i)
		values.emplace_back(PyFloat_AsDouble(PyList_GetItem(yvalues, i)));

	/* Create vectors */
	std::vector<double> xlocs;
	xlocs.reserve(psize);

	std::vector<double> ypeaks;
	ypeaks.reserve(psize);

	/* Locate turning points and add to peaks vector */
	for (size_t i = 1; i < values.size() - 1; ++i)
	{
		double delta1 = values[i] - values[i - 1];
		double delta2 = values[i + 1] - values[i];

		if ((std::abs(delta1) > 0.0) && (delta1 * delta2 <= 0.0))
		{
			if (ypeaks.size() == 0)
			{
				xlocs.emplace_back(positions[i - 1]);
				ypeaks.emplace_back(values[i - 1]);
			}
			else if (ypeaks.size() >= 2)
			{
				double point1 = ypeaks[ypeaks.size() - 2];
				double point2 = ypeaks.back();
				double point3 = values[i];

				if (point2 >= std::min<double>(point1, point3) && point2 <= std::max<double>(point1, point3))
				{
					xlocs.pop_back();
					ypeaks.pop_back();
				}
			}

			xlocs.emplace_back(positions[i]);
			ypeaks.emplace_back(values[i]);
		}
	}

	/* Add last value to peaks vector if delta is non-zero */
	double delta = values.back() - ypeaks.back();
	if (std::abs(delta) > 0.0)
	{
		xlocs.emplace_back(positions.back());
		ypeaks.emplace_back(values.back());
	}

	/* Output */
	psize = xlocs.size();
	PyObject* pylocs = PyList_New(psize);
	PyObject* pypeaks = PyList_New(psize);

	for (size_t i = 0; i < psize; ++i)
	{
		PyList_SetItem(pylocs, i, PyFloat_FromDouble(xlocs[i]));
		PyList_SetItem(pypeaks, i, PyFloat_FromDouble(ypeaks[i]));
	}

	PyTuple_SetItem(output, 0, pylocs);
	PyTuple_SetItem(output, 1, pypeaks);

	return output;
}


PyObject* __stdcall PyRainflow::PyRainflowCounting(PyObject* values, int algorithm, double tolerance, double cutoff)
{
	/* Define tolerance for rounding and cutoff for output */
	PyRainflow Rf(values, tolerance, cutoff);

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


PyObject* __stdcall PyRainflow::PyRainflowCountingRandomOrder(PyObject* histories, double tolerance, double cutoff, bool randomise, bool verify)
{
	/* Randomise order of stress histories on a separate thread */
	size_t nhistories = PyList_Size(histories);
	std::vector<size_t> ihistories;
	std::thread t0(&randomiseOrder, histories, &ihistories, randomise);

	PyRainflow RfAll(tolerance, cutoff);

	/* Perform rainflow counting on stress histories and save residue's for randomised order */
	std::vector< std::unique_ptr<std::vector<double>> > residues(nhistories);

	for (size_t i = 0; i < nhistories; ++i)
	{
		PyObject* item = PyList_GetItem(histories, i);
		size_t nmul = PyLong_AsUnsignedLongLong(PyList_GetItem(item, 0));
		PyObject* values = PyList_GetItem(item, 1);

		PyRainflow iRf(values, tolerance, cutoff);

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
		threads[it] = std::make_unique<std::thread>(std::thread(&PyRainflow::rainflow4PointConcurrent, tolerance, cutoff, it, &ihistories, &residues, &tresidues, &tcycles));

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
		PyRainflow RfVerify(tolerance, cutoff);

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
			std::vector<double> vpeaks = getPeaks(values);

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


PyObject* __stdcall PyRainflow::PyShift2ndLane(PyObject* histories_l1, PyObject* histories_l2, int step, int range, double tolerance)
	{
		/* Define tolerance for rounding  */
		PyRainflow Rf(tolerance);

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
		long long i_range = (long long) std::floor(his_length / i_step) * i_step;
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