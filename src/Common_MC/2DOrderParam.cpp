#include "2DOrderParam.h"

// ---> Constructors
void TwoDOrderParam::Initialize(int row_number, int col_number)
{
	m_nbOfMeas = 0;

	m_row_number = row_number;
	m_col_number = col_number;

	int neuroticBug = col_number*row_number;
	m_mean.resize(neuroticBug,0);

}

void TwoDOrderParam::SetFolder(const char folder[])
{
	m_folder = folder;
	m_parent = m_folder + "/Obj_conf.odat";
}

void TwoDOrderParam::SetType(const char dataType[])
{
	m_Type = dataType;
}

void TwoDOrderParam::CreateFiles(const char folder[],const char dataType[])
{
	SetFolder(folder);
	SetType(dataType);
}

void TwoDOrderParam::SetVt(double Vt)
{
	m_Vt = Vt;
}

// ---> Procedures
// >>>> Add data to the mean/binning
void TwoDOrderParam::AddData(vector<double> & measurement)
{
	for(int iii = 0; iii < m_row_number*m_col_number; ++iii)
	{
		m_mean[iii] += measurement[iii];
	}
}

// >>>> Reset the mean value for a new series of measurements
void TwoDOrderParam::RestartMean(int extraMeas)
{
	for(int iii = 0; iii < m_row_number*m_col_number; ++iii)
	{
		m_mean[iii] = m_mean[iii]*(double)m_nbOfMeas/(m_nbOfMeas + extraMeas);
	}
	m_nbOfMeas += extraMeas;
}

void TwoDOrderParam::RestartMean()
{
	for(int iii = 0; iii < m_row_number*m_col_number; ++iii)
	{
		m_mean[iii] = m_mean[iii]*(double)m_nbOfMeas;
	}
}

void TwoDOrderParam::RaiseMeasures(int extraMeas)
{
	m_nbOfMeas += extraMeas;
};
// >>>> Get mean values (normal and RMS)
void TwoDOrderParam::GetMean(vector<double> & output)
{
	for(int iii = 0; iii < m_row_number*m_col_number; ++iii)
	{
		output[iii] = m_mean[iii]/m_nbOfMeas;
	}
}

void TwoDOrderParam::ConvertMean()
{
	for(int iii = 0; iii < m_row_number*m_col_number; ++iii)
	{
		m_mean[iii] = m_mean[iii]*m_nbOfMeas;
	}
}

int TwoDOrderParam::GetNbOfMeas()
{
	return m_nbOfMeas;
}

void TwoDOrderParam::GetParent(string& filename)
{
	filename = m_parent;
}

void TwoDOrderParam::Export2D_Data(string& outputString)
{
	ofstream output(outputString.c_str());
	int counter = 0;

	for(int iii = 0; iii < m_row_number; ++iii)
	{
		output << m_x[iii] << " " << m_y[iii];
		for(int jjj = 0; jjj < m_col_number; ++jjj)
		{
			output << " " << m_mean[counter]/m_nbOfMeas;
			++counter;
		}
		output << endl;
	}

	output.close();
}

void TwoDOrderParam::Export2D_Matrix(string& outputString, int nxColMax, int nyRowMax)
{
	ofstream output(outputString.c_str());
	int counter = 0;

	for(int iii = 0; iii < nyRowMax; ++iii)
	{
		for(int jjj = 0; jjj < nxColMax; ++jjj)
		{
			output << m_mean[counter] << " ";
			++counter;
		}
		output << endl;
	}

	output.close();
};
void TwoDOrderParam::ReadPositions(string& posTable)
{
	ifstream input(posTable.c_str());

	m_x.resize(m_row_number,0);
	m_y.resize(m_row_number,0);

	for(int iii = 0; iii < m_row_number; ++iii)
	{
		input >> m_x[iii];
		input >> m_y[iii];
	}

	input.close();
}

void TwoDOrderParam::SetPositions(input_params& inputPars)
{
	m_x.resize(m_row_number,0);
	m_y.resize(m_row_number,0);

	if(inputPars.runType.compare("Part")==0)
	{
		int counter = 0;
		vector<double> init_x(inputPars.p + inputPars.nx - 1,-1);
		vector<double> init_y(inputPars.p + inputPars.nx - 1,-1);

		vector<int> lineStart(inputPars.p + inputPars.nx - 1,-1);
		vector<int> lineSize(inputPars.p + inputPars.nx - 1,-1);

		// > Set limits
		for(int iii = 0; iii < inputPars.p; ++iii)
		{
			lineStart[iii] = 0;
			init_x[iii] = 0;
			init_y[iii] = (inputPars.p - 1 - iii);
		}
		for(int iii = inputPars.p; iii < inputPars.p + inputPars.nx -1; ++iii)
		{
			lineStart[iii] = (iii - inputPars.p + 1);
			init_x[iii] = (iii - inputPars.p + 1)*sqrt(3)/2;
			init_y[iii] = -(iii - inputPars.p + 1)*0.5;
		}

		for(int iii = 0; iii < inputPars.nx; ++iii)
		{
			lineSize[iii] = (inputPars.ny + iii) - lineStart[iii];
		}
		for(int iii = inputPars.nx; iii < inputPars.p + inputPars.nx - 1; ++iii)
		{
			lineSize[iii] = (inputPars.nx - 1) + inputPars.ny - lineStart[iii];
		}

		for(int iii = 0; iii < inputPars.p + inputPars.nx - 1; ++iii)
		{
			for(int jjj = 0; jjj < lineSize[iii]; ++jjj)
			{
				m_x[counter] = init_x[iii] + jjj*sqrt(3)/2;
				m_y[counter] = init_y[iii] + 0.5*jjj;
				++counter;
			}
		}
	}
//	else if(inputPars.runType.compare("Moessner")==0)
//	{
//		int counter = 0;
//		double init_x = 0;
//		double init_y = 0;
//
//		for(int iii = 0; iii < inputPars.ny; ++iii)
//		{
//			init_x = 0;
//			init_y = inputPars.ny - iii;
//
//			for(int jjj = 0; jjj < inputPars.nx; ++jjj)
//			{
//				m_x[counter] = init_x + jjj*sqrt(3)/2;
//				m_y[counter] = init_y - 0.5*jjj;
//
//				++counter;
//			}
//		}
//	}
	else if(inputPars.runType.compare("Moessner")==0)
	{
		int counter = 0;
		double init_x = 0;
		double init_y = 0;

		for(int jjj = 0; jjj < inputPars.ny; ++jjj)
		{
			init_x = jjj*0.5;
			init_y = jjj*sqrt(3)/2;

			for(int iii = 0; iii < inputPars.nx; ++iii)
			{
				m_x[counter] = init_x + iii;
				m_y[counter] = init_y;
				++counter;
			}
		}
	}
//	else if(inputPars.SimType.compare("Manual")==0)
//	{
//		int counter = 0;
//
//		for(int iii = 0; iii < inputPars.ny/3; ++iii)
//		{
//			// First line
//
//			x = 4*iii*sqrt(3)/2;
//			y = 0;
//
//			for(int jjj = 0; jjj < inputPars.nx/2; ++jjj)
//			{
//				output << x << " " << y;
//				for(int kkk = 0; kkk < m_col_number; ++kkk)
//				{
//					dataPoint = m_mean[counter];
//
//					output << " " << dataPoint;
//					++counter;
//				}
//				output  << endl;
//
//				// Move to next
//				x += sqrt(3)/2;
//				y +=0.5;
//
//				output << x << " " << y;
//				for(int kkk = 0; kkk < m_col_number; ++kkk)
//				{
//					dataPoint = m_mean[counter];
//
//					output << " " << dataPoint;
//					++counter;
//				}
//				output  << endl;
//
//				// Move to the next
//
//				x += 0;
//				y += 1;
//			}
//
//			// Second line
//
//			x = (1 + 4*iii)*sqrt(3)/2;
//			y = -0.5;
//
//			for(int jjj = 0; jjj < inputPars.nx/2; ++jjj)
//			{
//				output << x << " " << y;
//				for(int kkk = 0; kkk < m_col_number; ++kkk)
//				{
//					dataPoint = m_mean[counter];
//
//					output << " " << dataPoint;
//					++counter;
//				}
//				output  << endl;
//
//				// Move to next
//				x += sqrt(3)/2;
//				y +=0.5;
//
//				output << x << " " << y;
//				for(int kkk = 0; kkk < m_col_number; ++kkk)
//				{
//					dataPoint = m_mean[counter];
//
//					output << " " << dataPoint;
//					++counter;
//				}
//				output  << endl;
//
//				// Move to the next
//
//				x += 0;
//				y += 1;
//			}
//
//			// First line
//
//			x = (3 + 4*iii)*sqrt(3)/2;
//			y = -0.5;
//
//			for(int jjj = 0; jjj < inputPars.nx/2; ++jjj)
//			{
//				output << x << " " << y;
//				for(int kkk = 0; kkk < m_col_number; ++kkk)
//				{
//					dataPoint = m_mean[counter];
//
//					output << " " << dataPoint;
//					++counter;
//				}
//				output  << endl;
//
//				// Move to the next
//				x += 0;
//				y += 1;
//
//				output << x << " " << y;
//				for(int kkk = 0; kkk < m_col_number; ++kkk)
//				{
//					dataPoint = m_mean[counter];
//
//					output << " " << dataPoint;
//					++counter;
//				}
//				output  << endl;
//
//				// Move to next
//				x += sqrt(3)/2;
//				y +=0.5;
//			}
//		}
//	}
}

void TwoDOrderParam::RectManualGetSuperlatticeNf(vector<double>& output,int nx, int ny, vector<int>& neighTable)
{
	int idx = 2*ny+6;
	int originalIdx = idx;

	for(uint iii = 0; iii < output.size();++iii)
	{
		output[iii] = 0;
	}

	for(int jjj = 1; jjj < nx/2-1; ++jjj)
	{
		// First line
		for(int iii = 1; iii < ny/3-1; ++iii)
		{
			// A
			output[0] += m_mean[idx];
			idx = neighTable[idxConv(6,idx,2)];

			// B
			output[1] += m_mean[idx];
			idx = neighTable[idxConv(6,idx,2)];

			// C
			output[2] += m_mean[idx];
			idx = neighTable[idxConv(6,idx,2)];
		}

		// Move to the side
		idx = neighTable[idxConv(6,originalIdx,1)];
		originalIdx = idx;

		// Second line
		for(int iii = 1; iii < ny/3-1; ++iii)
		{
			// C
			output[2] += m_mean[idx];
			idx = neighTable[idxConv(6,idx,2)];

			// A
			output[0] += m_mean[idx];
			idx = neighTable[idxConv(6,idx,2)];

			// B
			output[1] += m_mean[idx];
			idx = neighTable[idxConv(6,idx,2)];
		}

		// Move to the side
		idx = neighTable[idxConv(6,originalIdx,0)];
		originalIdx = idx;
	}

	for(uint iii = 0; iii < output.size();++iii)
	{
		output[iii] = output[iii]/(2.*(nx/2 - 2)*(ny/3 - 2));
	}

// 	sort(output.begin(), output.end());
}

void TwoDOrderParam::RhombusManualGetSuperlatticeNf(vector<double>& output, vector<int>& SubSize,int nx, int ny, vector<int>& neighTable)
{

	int idx;

	bool countA = true;
	bool countB = true;
	bool countC = true;

	for(uint iii = 0; iii < output.size();++iii)
	{
		output[iii] = 0;
		SubSize[iii] = 0;
	}

	int L = nx*ny;
	for(int lll = 0; lll < L/3; ++lll)
	{
		countA = true;
		countB = true;
		countC = true;

		idx = 3*lll;
		for(int kkk = 0; kkk < 6; ++kkk)
		{
			if(neighTable[idxConv(6,idx,kkk)]>L-1)
			{
				countA = false;
				break;
			}
		}

		if(countA)
		{
			++SubSize[0];
			output[0] += m_mean[idx];
		}

		idx = 3*lll+1;
		for(int kkk = 0; kkk < 6; ++kkk)
		{
			if(neighTable[idxConv(6,idx,kkk)]>L-1)
			{
				countB = false;
				break;
			}
		}

		if(countB)
		{
			++SubSize[1];
			output[1] += m_mean[idx];
		}

		idx = 3*lll+2;
		for(int kkk = 0; kkk < 6; ++kkk)
		{
			if(neighTable[idxConv(6,idx,kkk)]>L-1)
			{
				countC = false;
				break;
			}
		}

		if(countC)
		{
			++SubSize[2];
			output[2] += m_mean[idx];
		}
	}

	for(uint iii = 0; iii < output.size();++iii)
	{
		output[iii] = output[iii]/SubSize[iii];
	}

// 	sort(output.begin(), output.end());
}

void TwoDOrderParam::GetSuperlatticeNf(vector<double>& output,int nx, int ny)
{
	int idx = -1;
	for(uint iii = 0; iii < output.size();++iii)
	{
		output[iii] = 0;
	}

	for(int iii = 0; iii < ny/3; ++iii)
	{
		for(int jjj = 0; jjj < nx/3; ++jjj)
		{
			// First line
			// A
			idx = idxConv(nx,3*iii,3*jjj);
			output[0] += m_mean[idx];

			// B
			idx = idxConv(nx,3*iii,3*jjj+1);
			output[1] += m_mean[idx];

			// C
			idx = idxConv(nx,3*iii,3*jjj+2);
			output[2] += m_mean[idx];

			// Second line
			// C
			idx = idxConv(nx,3*iii+1,3*jjj);
			output[2] += m_mean[idx];

			// A
			idx = idxConv(nx,3*iii+1,3*jjj+1);
			output[0] += m_mean[idx];

			// B
			idx = idxConv(nx,3*iii+1,3*jjj+2);
			output[1] += m_mean[idx];

			// Third line
			// B
			idx = idxConv(nx,3*iii+2,3*jjj);
			output[1] += m_mean[idx];

			// C
			idx = idxConv(nx,3*iii+2,3*jjj+1);
			output[2] += m_mean[idx];

			// A
			idx = idxConv(nx,3*iii+2,3*jjj+2);
			output[0] += m_mean[idx];
		}
	}

	for(uint iii = 0; iii < output.size();++iii)
	{
		output[iii] = 3.*output[iii]/m_row_number;
	}

// 	sort(output.begin(), output.end());
}


