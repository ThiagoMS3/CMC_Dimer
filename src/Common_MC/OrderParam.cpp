#include "OrderParam.h"

// ---> Constructors
OrderParam::OrderParam(input_params& input)
{
	m_data.resize(2*input.minBinNb,0);
	m_bin.resize(2*input.minBinNb,0);
	m_error.resize(0);
	m_errorBinSize.resize(0);

	m_AC.resize(0);
	m_IntegAC = 0;
	m_ACInterval = input.AutocorrLength;

	m_minBinNb = input.minBinNb;
	m_totalPoints =0;
	m_binCount = 0;
	m_binSize = 1;
	m_fillBin = 0;
	m_mean = 0;
	m_nbOfMeas = 0;

	m_corrSize = 0;
}

void OrderParam::SetFolder(const char folder[])
{
	m_folder = folder;
}

void OrderParam::SetType(const char dataType[])
{
	m_Type = dataType;
}

void OrderParam::SetCorrLength(int size)
{
	m_corrSize = size;
	m_corr.resize(m_corrSize,0);
}

int OrderParam::GetCorrLength()
{
	return m_corrSize;
}

void OrderParam::CreateFiles()
{
	m_outErr  = m_folder + "/err_" + m_Type;
	m_outAC   = m_folder + "/AC_" + m_Type + ".dat";

	output.open(m_outAC.c_str(),ios::trunc);
	output.close();
}

void OrderParam::CreateFiles(const char folder[],const char dataType[])
{
	SetFolder(folder);
	SetType(dataType);

	m_outErr  = m_folder + "/err_" + m_Type;
	m_outAC   = m_folder + "/AC_" + m_Type + ".dat";

	output.open(m_outAC.c_str(),ios::trunc);
	output.close();
}

void OrderParam::SetVt(double Vt)
{
	m_Vt = Vt;
}


// ---> Procedures
// >>>> Add data to the mean/binning
void OrderParam::AddData(double measurement)
{
	m_data[m_binCount] 	+= measurement/m_binSize;
	m_mean 			+= measurement;
	++m_fillBin;

	// > If fillBin == binSize, then we must go to the next bin
	if(m_fillBin==m_binSize)
	{
		++m_binCount;
		m_fillBin = 0;

		// > If binCounter == 2*minBin, then we have to compact the binning data
		if(m_binCount==2*m_minBinNb)
		{
			// Compress
			for(int iii = 0; iii < m_minBinNb; ++iii)
			{
				m_data[iii] = (m_data[2*iii] + m_data[2*iii+1])/2;
			}

			// Resize bins
			m_binSize = 2*m_binSize;

			// Clear bins
			m_binCount = m_minBinNb;
			for(int iii = m_minBinNb; iii < 2*m_minBinNb; ++iii)
			{
				m_data[iii] = 0;
			}
		}
	}
}

void OrderParam::AddCorr(vector<double>& inp_corr)
{
	for(int iii = 0; iii < m_corrSize; ++iii)
	{
		m_corr[iii] += inp_corr[iii];
	}
}

// >>>> Do binning
void OrderParam::CalculateBinning(int& totalPoints)
{
	int l_binNumber = m_binCount/2;
	m_totalPoints = log(m_binCount)/log(2);
	totalPoints = m_totalPoints;

	m_bin.resize(2*l_binNumber,0);
	m_error.resize(m_totalPoints,0);
	m_errorBinSize.resize(m_totalPoints,1);

	for(int iii = 0; iii < 2*l_binNumber; ++iii)
	{
		m_bin[iii] = m_data[iii];
	}

	double sum = 0;
	double sumSqr = 0;

	for(int nnn = 0; nnn < m_totalPoints; ++nnn)
	{
		sum = 0;
		sumSqr = 0;

		for(int iii = 0; iii < l_binNumber; ++iii)
		{
			sum += m_bin[2*iii] + m_bin[2*iii + 1];
			sumSqr += m_bin[2*iii]*m_bin[2*iii] + m_bin[2*iii + 1]*m_bin[2*iii + 1];
			m_bin[iii] = (m_bin[2*iii] + m_bin[2*iii + 1])/2;
		}

		m_error[nnn] = sqrt( abs( sumSqr/(2.*l_binNumber) - pow(sum/(2.*l_binNumber),2) )/(2.*l_binNumber) );
		m_errorBinSize[nnn] = m_binSize*pow(2,nnn);

		l_binNumber = l_binNumber/2;
	}
}

void OrderParam::CalculateAC()
{
	m_AC.resize(m_ACInterval+1);
	double Qti = 0;
	double Qi = 0;
	double Qt = 0;

	int Npoints;
	double var = 0;
	m_IntegAC = 0.5;

	// Autocorrelation for the point 0 = variance
	for(int kkk = 0; kkk < m_binCount; ++kkk)
	{
		Qti += m_data[kkk]*m_data[kkk];
		Qi += m_data[kkk];
	}
	Qti = Qti/m_binCount;
	Qi = Qi/m_binCount;

	var = (Qti - Qi*Qi);
	m_AC[0] = 1;

	// Now, autocorrelation for the other points
	for(int offset = 1; offset < m_ACInterval + 1; ++offset)
	{
		Npoints = m_binCount - offset;
		m_AC[offset] = 0;
		Qti = 0;
		Qi = 0;
		Qt = 0;

		for(int kkk = 0; kkk < Npoints; ++kkk)
		{
			Qti += m_data[kkk]*m_data[kkk+offset];
			Qi += m_data[kkk];
			Qt += m_data[kkk+offset];
		}
		Qti = Qti/Npoints;
		Qi = Qi/Npoints;
		Qt = Qt/Npoints;

		// Autocorrelation
		if(abs(var)>1E-10)
			m_AC[offset] = (Qti - Qt*Qi)/var;
		else
			m_AC[offset] = (Qti - Qt*Qi);

		// Integrated autocorrelation
		m_IntegAC += m_AC[offset];
	}
}

void OrderParam::CalculateError(int& totalPoints)
{
	CalculateBinning(totalPoints);
	CalculateAC();
}

// >>>> Do binning
void OrderParam::CalculateRMSBinning(int& totalPoints)
{
	int l_binNumber = m_binCount/2;
	m_totalPoints = log(m_binCount)/log(2);
	totalPoints = m_totalPoints;

	m_bin.resize(2*l_binNumber,0);
	m_error.resize(m_totalPoints,0);
	m_errorBinSize.resize(m_totalPoints,1);

	for(int iii = 0; iii < 2*l_binNumber; ++iii)
	{
		m_bin[iii] = sqrt(m_data[iii]);
	}

	double sum = 0;
	double sumSqr = 0;

	for(int nnn = 0; nnn < m_totalPoints; ++nnn)
	{
		sum = 0;
		sumSqr = 0;

		for(int iii = 0; iii < l_binNumber; ++iii)
		{
			sum += m_bin[2*iii] + m_bin[2*iii + 1];
			sumSqr += m_bin[2*iii]*m_bin[2*iii] + m_bin[2*iii + 1]*m_bin[2*iii + 1];
			m_bin[iii] = (m_bin[2*iii] + m_bin[2*iii + 1])/2;
		}

		m_error[nnn] = sqrt( abs( sumSqr/(2.*l_binNumber) - pow(sum/(2.*l_binNumber),2) )/(2.*l_binNumber) );
		m_errorBinSize[nnn] = m_binSize*pow(2,nnn);

		l_binNumber = l_binNumber/2;
	}
}

void OrderParam::CalculateRMSAC()
{
	m_AC.resize(m_ACInterval+1);
	double Qti = 0;
	double Qi = 0;
	double Qt = 0;

	int Npoints;
	double var = 0;
	m_IntegAC = 0.5;

	// Autocorrelation for the point 0 = variance
	for(int kkk = 0; kkk < m_binCount; ++kkk)
	{
		Qti += m_data[kkk];
		Qi += sqrt(m_data[kkk]);
	}
	Qti = Qti/m_binCount;
	Qi = Qi/m_binCount;

	var = (Qti - Qi*Qi);
	m_AC[0] = 1;

	// Now, autocorrelation for the other points
	for(int offset = 1; offset < m_ACInterval + 1; ++offset)
	{
		Npoints = m_binCount - offset;
		m_AC[offset] = 0;
		Qti = 0;
		Qi = 0;
		Qt = 0;

		for(int kkk = 0; kkk < Npoints; ++kkk)
		{
			Qti += sqrt(m_data[kkk])*sqrt(m_data[kkk+offset]);
			Qi += sqrt(m_data[kkk]);
			Qt += sqrt(m_data[kkk+offset]);
		}
		Qti = Qti/Npoints;
		Qi = Qi/Npoints;
		Qt = Qt/Npoints;

		// Autocorrelation
		if(abs(var)>1E-10)
			m_AC[offset] = (Qti - Qt*Qi)/var;
		else
			m_AC[offset] = (Qti - Qt*Qi);

		// Integrated autocorrelation
		m_IntegAC += m_AC[offset];
	}
}

void OrderParam::CalculateRMSError(int& totalPoints)
{
	CalculateRMSBinning(totalPoints);
	CalculateRMSAC();
}

// >>>> Reset the mean value for a new series of measurements
void OrderParam::RestartMean(int extraMeas)
{
	m_mean = m_mean*m_nbOfMeas;
	m_nbOfMeas += extraMeas;
	m_mean = m_mean/m_nbOfMeas;
}

void OrderParam::RestartMean()
{
	m_mean = m_mean*m_nbOfMeas;
}

void OrderParam::RestartCorr(int extraMeas)
{
	int oldMeas = m_nbOfMeas;
	m_nbOfMeas += extraMeas;
	for(int iii = 0; iii < m_corrSize; ++iii)
	{
		m_corr[iii] = (m_corr[iii]*oldMeas)/m_nbOfMeas;
	}
}

void OrderParam::RestartCorr()
{
	for(int iii = 0; iii < m_corrSize; ++iii)
	{
		m_corr[iii] = m_corr[iii]*m_nbOfMeas;
	}
}

void OrderParam::RaiseMeasures(int extraMeas)
{
	m_nbOfMeas += extraMeas;
};

// ---> Getters
// // >>>> Get parameters
// void OrderParam::GetParam(vector<double>& data, int& minBinNb, int& binCount, int& binSize, int& fillBin, double& mean, int& nbOfMeas)
// {
// 	data.resize(2*m_minBinNb);
// 	data = m_data;
//
// 	minBinNb = m_minBinNb;
// 	binCount = m_binCount;
// 	binSize = m_binSize;
// 	fillBin = m_fillBin;
// 	mean = m_mean;
// 	nbOfMeas = m_nbOfMeas;
// }

// >>>> Get errors
void OrderParam::GetError(vector<double>& error, vector<double>& sizes)
{
	error.resize(m_error.size());
	sizes.resize(m_errorBinSize.size());

	error = m_error;
	sizes = m_errorBinSize;
}

// >>>> Get mean values (normal and RMS)
double OrderParam::GetMean()
{
	return m_mean/m_nbOfMeas;
}

double OrderParam::ConvertMean()
{
	return m_mean*m_nbOfMeas;
}

int OrderParam::GetNbOfMeas()
{
	return m_nbOfMeas;
}

void OrderParam::GetCorr(vector<double>& outCorr)
{
	for(unsigned int iii = 0; iii < m_corr.size(); ++iii)
	{
		outCorr[iii] = m_corr[iii]/m_nbOfMeas;
	}
}

void OrderParam::GetCorr(vector<double>& outCorr,int length)
{
	for(int iii = 0; iii < length; ++iii)
	{
		outCorr[iii] = m_corr[iii];
	}
}

void OrderParam::CalculateRMSMean()
{
	m_mean = sqrt(m_mean);
}

void OrderParam::ResetRMSMean()
{
	m_mean = m_mean*m_mean;
}
// ---> I/O
// void OrderParam::PrintMean()
// {
// 	output.open(m_outMean.c_str(), ios::out | ios::app);
// 	output << m_Vt << " " << m_mean << endl;
// 	output.close();
// }

// void OrderParam::PrintRMSMean()
// {
// 	output.open(m_outMean.c_str(), ios::out | ios::app);
// 	output << m_Vt << " " << sqrt(m_mean) << endl;
// 	output.close();
// }

void OrderParam::PrintError()
{
	PrintBin();
	PrintAC();
}

void OrderParam::PrintAC()
{
	output.open(m_outAC.c_str(), ios::out | ios::app);
	for(int iii = 0; iii < m_ACInterval + 1; ++iii)
	{
		output << " " << m_AC[iii];
	}
	output << endl;
	output.flush();
	output.close();
}

// void OrderParam::PrintIntegAC()
// {
// 	output.open(m_outACInteg.c_str(), ios::out | ios::app);
// 	output << m_Vt << " " << m_IntegAC << endl;
// 	output.flush();
// 	output.close();
// }

void OrderParam::PrintBin()
{
	string fileErr;

	// Output file name
	fileErr = m_outErr + ".dat";

	// Print
	output.open(fileErr.c_str(),ios::trunc);
	for(int iii = 0; iii < m_totalPoints; ++iii)
	{
		output << m_errorBinSize[iii] << " " << m_error[iii] << endl;
	}
	output.flush();
	output.close();
}

void OrderParam::PrintDebug(ostream& out)
{
	out << " ### DEBUG Print OrderParam ###" << endl;
	out << " Type                   = " << m_Type << endl;
	out << " minBinNb               = " << m_minBinNb << endl;
	out << " TotalPoints            = " << m_totalPoints << endl;
	out << " binCount               = " << m_binCount << endl;
	out << " binSize                = " << m_binSize << endl;
	out << " fillBin                = " << m_fillBin << endl;
	out << " ACInterval             = " << m_ACInterval << endl << endl;

	out << " Mean                   = " << m_mean << endl;
	out << " nbOfMeas               = " << m_nbOfMeas << endl;

	out << " V/t                    = " << m_Vt << endl;
}
