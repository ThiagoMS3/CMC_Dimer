#include "CmplxOrderParam.h"

// ---> Constructors
CmplxOrderParam::CmplxOrderParam(input_params& input,int histArgSize, int histMagSize, double histMaxMag)
{
	m_histArgSize = histArgSize;
	m_histMagSize = histMagSize;
	m_histMaxMag  = histMaxMag;

	m_hist.resize(histMagSize*histArgSize,0);
	m_histDeltaArg = 2*M_PI/histArgSize;
	m_histDeltaMag = histMaxMag/histMagSize;

	m_mean = 0;
	m_nbOfMeas = 0;
}

void CmplxOrderParam::SetFolder(const char folder[])
{
	m_folder = folder;
}

void CmplxOrderParam::SetType(const char dataType[])
{
	m_Type = dataType;
}

void CmplxOrderParam::CreateFiles(const char folder[],const char dataType[])
{
	SetFolder(folder);
	SetType(dataType);
}

void CmplxOrderParam::SetVt(double Vt)
{
	m_Vt = Vt;
}


// ---> Procedures
// >>>> Add data to the mean/binning
void CmplxOrderParam::AddData(complex<double> measurement)
{
        // Add mean value
        m_mean += measurement;

        // Add point to histogram
        int posArg = floor((M_PI+arg(measurement))/m_histDeltaArg);
        int posMag = min(floor(abs(measurement)/m_histDeltaMag),m_histMagSize-1);
		uint idx = idxConv(m_histArgSize,posMag,posArg);
		if(idx > m_hist.size()-1)
		{
			cout << "Warning !!! " << idx << endl;
		}
        ++m_hist[idxConv(m_histArgSize,posMag,posArg)];
}

// >>>> Reset the mean value for a new series of measurements
void CmplxOrderParam::RestartMean(int extraMeas)
{
	m_mean = m_mean*(double)m_nbOfMeas;
	m_nbOfMeas += extraMeas;
	m_mean = m_mean/(double)m_nbOfMeas;
}

void CmplxOrderParam::RestartMean()
{
	m_mean = m_mean*(double)m_nbOfMeas;
}

void CmplxOrderParam::RaiseMeasures(int extraMeas)
{
	m_nbOfMeas += extraMeas;
};

// >>>> Get mean values (normal and RMS)
complex<double> CmplxOrderParam::GetMean()
{
	return m_mean/(double)m_nbOfMeas;
};

void			CmplxOrderParam::GetHist(vector<int> & histOut)
{
	histOut.resize(m_hist.size());
	histOut = m_hist;
};

complex<double> CmplxOrderParam::ConvertMean()
{
	return m_mean*(double)m_nbOfMeas;
};

int CmplxOrderParam::GetNbOfMeas()
{
	return m_nbOfMeas;
}
