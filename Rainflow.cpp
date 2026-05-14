#include <algorithm>
#include <cmath>
#include <iostream>
#include <Rainflow.h>


Rainflow::Rainflow(double tolerance, double cutoff)
{
	setTolerance(tolerance);
	setCutOff(cutoff);
};


Rainflow::Rainflow(std::vector<double>& values, double tolerance, double cutoff)
{
	setTolerance(tolerance);
	setCutOff(cutoff);

	peaks = getPeaks(values);
};


Rainflow::Rainflow(std::vector<double>& values, std::map<std::pair<double, double>, double>& counts, double tolerance, double cutoff)
{
	setTolerance(tolerance);
	setCutOff(cutoff);

	peaks = getPeaks(values);
	cycles = counts;
};


Rainflow::~Rainflow() {};


void Rainflow::setTolerance(double tolerance)
{
	roundingTolerance = std::max<double>(1e-16, tolerance);
};


void Rainflow::setCutOff(double cutoff)
{
	deltaCutOff = std::max<double>(1e-16, cutoff);
};


double Rainflow::round(double value)
{
	return std::round(value / roundingTolerance) * roundingTolerance;
}


void Rainflow::addCycles(std::map<std::pair<double, double>, double>& counts)
{
	for (std::map<std::pair<double, double>, double>::iterator itr = counts.begin(); itr != counts.end(); ++itr)
		cycles[itr->first] += itr->second;

	return;
}


bool Rainflow::isPeaksOnly(std::vector<double>& values)
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


bool Rainflow::isPeaksOnly()
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


std::vector<double> Rainflow::getPeaks(std::vector<double>& values)
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


void Rainflow::setPeaks(std::vector<double>& values)
{
	peaks = getPeaks(values);
	return;
}


void Rainflow::filterPeaks()
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


void Rainflow::rotatePeaks()
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


std::pair<std::vector<double>, std::vector<double>> Rainflow::getPeakLocations(std::vector<double> positions, std::vector<double> values)
{
	/* Get sizes */
	size_t psize = positions.size();

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

	return std::pair<std::vector<double>, std::vector<double>>(xlocs, ypeaks);
}


void Rainflow::rainflow3Points()
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


void Rainflow::rainflow4Points(bool process_residue, size_t multiplier)
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


void Rainflow::rainflow4PointConcurrent(double tolerance, double cutoff, unsigned it,
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

	Rainflow iRf(tolerance, cutoff);

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