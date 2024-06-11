#include"CModulate.h"
//High order modulation table 
//Gray map
static float table_qpsk[2] = { -0.707107f, 0.707107f };
static float	table_16qam[4] = { -0.316228f, -0.948683f,  0.316228f,  0.948683f };
static float	table_64qam[8] = { -0.462910f, -0.154303f, -0.771517f, -1.08012f,  0.462910f,  0.154303f,  0.771517f,  1.08012f };
static float    table_256qam[16] = { -0.383482f, -0.536875f, -0.230089f, -0.076696f, -0.843661f, -0.690268f, -0.997054f, -1.150447f, 0.383482f, 0.536875f, 0.230089f, 0.076696f, 0.843661f, 0.690268f, 0.997054f, 1.150447f };
CModulate::CModulate()
	: ModSeq(nullptr)
	, DemodSeq(nullptr)
	,InterLeaveSeq(nullptr)
	,ILSeq(nullptr)
	,DILSeq(nullptr)
	,DeInterLeaveSeq(nullptr)
	, SourceLen(0)
	, SymbolLen(0)
	, BPSKModSeq(nullptr)

{

}

CModulate::~CModulate()
{
	if (ModSeq)
	{
		free(ModSeq);
	}
	if (DemodSeq)
	{
		free(DemodSeq);
	}
	if (BPSKModSeq)
	{
		free(BPSKModSeq);
	}
	if (InterLeaveSeq)
	{
		free(InterLeaveSeq);
	}
	if (DeInterLeaveSeq)
	{
		free(DeInterLeaveSeq);
	}
	if (ILSeq)
	{
		free(ILSeq);
	}
	if (DILSeq)
	{
		free(DILSeq);
	}

}

void CModulate::Initial(unsigned long codelen)
{
	SourceLen = codelen;
	DemodSeq = (float*)vec_malloc(sizeof(float)*SourceLen);//new float[SourceLen];
	InterLeaveSeq = (int8_t*)vec_malloc(sizeof(int8_t) * 32 * BitsOverChannel);
	DeInterLeaveSeq = (float*)vec_malloc(sizeof(float) * 32 * BitsOverChannel);
	ILSeq = (int8_t*)vec_malloc(sizeof(int8_t) * 32 * BitsOverChannel);
	DILSeq = (float*)vec_malloc(sizeof(float) * 32 * BitsOverChannel);
	switch (ModulationType)
	{
	case 1:
		SymbolLen = SourceLen;
		BPSKModSeq = (float*)vec_malloc(sizeof(float)*SymbolLen);
		break;
	case 2:
		SymbolLen = SourceLen / 2;
		ModSeq = (MKL_Complex8*)vec_malloc(sizeof(MKL_Complex8)*SymbolLen);//new MKL_Complex8[SymbolLen];
		break;
	case 4:
		SymbolLen = SourceLen / 4;
		ModSeq = (MKL_Complex8*)vec_malloc(sizeof(MKL_Complex8)*SymbolLen);// new MKL_Complex8[SymbolLen];
		break;
	case 6:
		SymbolLen = SourceLen / 6;
		ModSeq = (MKL_Complex8*)vec_malloc(sizeof(MKL_Complex8)*SymbolLen);//new MKL_Complex8[SymbolLen];
		break;
	case 8:
		SymbolLen = SourceLen / 8;
		ModSeq = (MKL_Complex8*)vec_malloc(sizeof(MKL_Complex8)*SymbolLen);//new MKL_Complex8[SymbolLen];
		break;

	default:
		std::cout << "No such CQI\n" << std::endl;
		std::cin.get();
		break;
	}

}

void CModulate::BeforeModulationInterleaver(int8_t *codeseq)
{
	int m, i, j;
	int k = 0;
	
	/*恢复为32*（信息+校验）的格式*/
	//the information bits
	for (i = 0; i < 32; i++)
	{
		for(j=0;j< (NmoinsK - _PunctureBits - _ShortenBits);j++)
	    {
			ILSeq[i*BitsOverChannel + j] = codeseq[i*(NmoinsK - _PunctureBits - _ShortenBits) + j];
	    }
	 }	
	//for (i = 0; i < (NmoinsK - _PunctureBits - _ShortenBits); i += 1)
		//{
		//	for (j = 0; j < 32; j += 1)
		//	{
		//		ILSeq[j *(BitsOverChannel) + i] = codeseq[32 * i + j];//ILSeq(0)=outputbits(0) ILSeq(Bitsoverchannel)=outputbits(1) i=1 IL(1)=outputbits(32) j=3i=1 The first bit of the third frame outputbits(35)=ILseq(3*bitsoverchannel+1)
		//																						
		//	}
		//}
	
	// The check Bits
	for (i = 0; i < 32; i++)
	{
		for (j = 0; j < (_NoCheck); j++)
		{
			ILSeq[i*BitsOverChannel + (NmoinsK - _PunctureBits - _ShortenBits) + j] = codeseq[32* (NmoinsK - _PunctureBits - _ShortenBits)+i*_NoCheck + j];
		}
	}

		//for (i = 0; i < _NoCheck; ++i)
		//{
		//	for (j = 0; j < 32; ++j)
		//	{

		//		ILSeq[j*BitsOverChannel+ (NmoinsK - _PunctureBits - _ShortenBits)+i] = codeseq[(NmoinsK - _PunctureBits - _ShortenBits) * 32 + 32 * i + j];//j=3 i=5 The fifth check of the third frame ILSeq(3*bitsoverchannel+(NmoinsK - _PunctureBits - _ShortenBits)+5)=ptr(163)
		//	}
		//}
	
	/*Interleave each frame of data separately*/
	
	for (m = 0; m < 32; m++)
		{

			for (j = 0; j < (BitsOverChannel / InterleaveModType); j++)
				{
					for (i = 0; i < InterleaveModType; i++)
					{
						InterLeaveSeq[k] = ILSeq[BitsOverChannel / InterleaveModType * i + j + BitsOverChannel * m];
						k++;
					}
				}
		}


}



void CModulate::AfterDeModulationDeInterleaver()
{
	/*The format of DIL[k] after deinterleaving is 32*(information+check)*/
	int m, i, j;
	int k = 0;
	for (m = 0; m < 32; m++)
		{
				for (j = 0; j < InterleaveModType; j++)
				{
					for (i = 0; i < (BitsOverChannel / InterleaveModType); i++)
					{
						DILSeq[k] = DemodSeq[InterleaveModType*i + j + BitsOverChannel * m];
						k++;
					}
				 }
		 }

	/*Restore [32*(information+check)] to [(32*information)+(32*check)]*/
	//the information bits
	for (i = 0; i < 32; ++i)

	{
		for (j = 0; j < (NmoinsK - _PunctureBits - _ShortenBits); ++j)
		{

			DeInterLeaveSeq[i*(NmoinsK - _PunctureBits - _ShortenBits)+j] = DILSeq[BitsOverChannel * i+ j];//
		}
	}
		//for (i = 0; i < (NmoinsK - _PunctureBits - _ShortenBits); i += 1)
		//{
		//	for (j = 0; j < 32; j += 1)
		//	{
		//		DeInterLeaveSeq[32 * i + j] = DILSeq[j *BitsOverChannel + i];//j=1 i=1 The first bit of the first frame should be stored in DIL (Bitsoverchannel+1). Stored in DeInterLeave (33) j=0 i=5 The fifth bit of the 0th frame is stored in DIL (5) DeInterLeave (160)
		//																						
		//	}
		//}
	
	// The check Bits
	for (i = 0; i < 32; ++i)

	{
		for (j = 0; j < _NoCheck ; ++j)
		{

			DeInterLeaveSeq[32* (NmoinsK - _PunctureBits - _ShortenBits)+i*_NoCheck + j] = DILSeq[ i * BitsOverChannel  + (NmoinsK - _PunctureBits - _ShortenBits)+ j];//
		}
	}
			
	//for (i = 0; i < _NoCheck; ++i)
			//{
			//	for (j = 0; j < 32; ++j)
			//	{

			//		DeInterLeaveSeq[(NmoinsK - _PunctureBits - _ShortenBits) * 32 + i * 32 + j] = DemodSeq[BitsOverChannel * j+(NmoinsK - _PunctureBits - _ShortenBits) + i];//
			//	}
			//}
}



void CModulate::Modulation(int8_t* ReceivedSeq)
{
	unsigned long i, j;
	unsigned long tmp_i, tmp_q;
	int half_sym;
	float* map_table = NULL;
	half_sym = ModulationType / 2;
	switch (ModulationType)
	{

	case 2:
		map_table = table_qpsk;
		break;
	case 4:
		map_table = table_16qam;
		break;
	case 6:
		map_table = table_64qam;
		break;
	case 8:
		map_table = table_256qam;
		break;
	default:
		std::cout << "Invalid modulation type: 1:BPSK 2:QPSK, 4:16-QAM, 6:64-QAM! \n";
		break;

	}

	for (i = 0; i < SymbolLen; i++)//16QAM SymbolLen = SourceLen / 4 
								   //SourceLen = codelen=32*Bitsoverchannel
								   //R11/12 9216 bits are modulated into 2304 complex numbers
	{
		tmp_i = 0;
		tmp_q = 0;

		for (j = 0; j < half_sym; j++)//16QAM half_sym=2 For example, 11001100

		{

			tmp_i = tmp_i + (ReceivedSeq[i*ModulationType + 2 * j] << (half_sym - j - 1));// tmp_i=0+0<<1(shift left one bit)tmp_i=0+1<<1(shift left one bit)=2
			tmp_q = tmp_q + (ReceivedSeq[i*ModulationType + 2 * j + 1] << (half_sym - j - 1));
		}

		ModSeq[i].real = map_table[tmp_i];//The real part of the first complex number is the number in the table
		ModSeq[i].imag = map_table[tmp_q];

	}

}
/*
	Demodulation:
	Max-Log-Map Demodulation
*/

void CModulate::Demodulation(MKL_Complex8* received_symbols)
{
	unsigned long i;
	switch (ModulationType)
	{
	case 2:
		for (i = 0; i < SymbolLen; i++) {

			DemodSeq[i*ModulationType] = received_symbols[i].real;
			DemodSeq[i*ModulationType + 1] = received_symbols[i].imag;

		}
		break;
	case 4:
		for (i = 0; i < SymbolLen; i++) {


			DemodSeq[i*ModulationType] = received_symbols[i].real;
			DemodSeq[i*ModulationType + 1] = received_symbols[i].imag;

			DemodSeq[i*ModulationType + 2] = fabs(received_symbols[i].real) - 0.6324555;
			DemodSeq[i*ModulationType + 3] = fabs(received_symbols[i].imag) - 0.6324555;

		}
		//for (i = 0; i < SymbolLen; i++) 
		//{
		//received_symbols[i].real *= sqrt(10.0);
		//received_symbols[i].imag *= sqrt(10.0);
		////b0
		//DemodSeq[i*ModulationType] = received_symbols[i].real;
		////b1
		//DemodSeq[i*ModulationType + 1] = received_symbols[i].imag;
		////b2
		//// Calc b2				
		//if (DemodSeq[i*ModulationType] >= 2)
		//	DemodSeq[i*ModulationType + 2] = DemodSeq[i*ModulationType]-2;
		//else if (DemodSeq[i*ModulationType] < -2)
		//	DemodSeq[i*ModulationType + 2] = 0-(DemodSeq[i*ModulationType] + 2);
		//else if (DemodSeq[i*ModulationType] < 2 && DemodSeq[i*ModulationType] >= 0)
		//	DemodSeq[i*ModulationType + 2] = DemodSeq[i*ModulationType] - 2;
		//else
		//	DemodSeq[i*ModulationType + 2] = 0-(DemodSeq[i*ModulationType] + 2);

		//// Calc b3				
		//if (DemodSeq[i*ModulationType + 1] >= 2)
		//	DemodSeq[i*ModulationType + 3] = DemodSeq[i*ModulationType + 1] - 2;
		//else if (DemodSeq[i*ModulationType + 1] < -2)
		//	DemodSeq[i*ModulationType + 3] = 0-(DemodSeq[i*ModulationType + 1] + 2);
		//else if (DemodSeq[i*ModulationType + 1] < 2 && DemodSeq[i*ModulationType + 1] >= 0)
		//	DemodSeq[i*ModulationType + 3] = DemodSeq[i*ModulationType + 1] - 2;
		//else
		//	DemodSeq[i*ModulationType + 3] = 0-(DemodSeq[i*ModulationType + 1] + 2);
		//}

		break;
	case 6:
		for (i = 0; i < SymbolLen; i++) {

			DemodSeq[i*ModulationType] = received_symbols[i].real;
			DemodSeq[i*ModulationType + 1] = received_symbols[i].imag;

			DemodSeq[i*ModulationType + 2] = fabs(DemodSeq[i*ModulationType]) - 0.6172134;
			DemodSeq[i*ModulationType + 3] = fabs(DemodSeq[i*ModulationType + 1]) - 0.6172134;

			DemodSeq[i*ModulationType + 4] = fabs(DemodSeq[i*ModulationType + 2]) - 0.3086067;
			DemodSeq[i*ModulationType + 5] = fabs(DemodSeq[i*ModulationType + 3]) - 0.3086067;

		}
		break;

	case 8:
		for (i = 0; i < SymbolLen; i++) {

			DemodSeq[i*ModulationType] = received_symbols[i].real;
			DemodSeq[i*ModulationType + 1] = received_symbols[i].imag;

			DemodSeq[i*ModulationType + 2] = fabs(DemodSeq[i*ModulationType]) - 0.613568;
			DemodSeq[i*ModulationType + 3] = fabs(DemodSeq[i*ModulationType + 1]) - 0.613568;

			DemodSeq[i*ModulationType + 4] = fabs(DemodSeq[i*ModulationType + 2]) - 0.306784;
			DemodSeq[i*ModulationType + 5] = fabs(DemodSeq[i*ModulationType + 3]) - 0.306784;

			DemodSeq[i*ModulationType + 6] = fabs(DemodSeq[i*ModulationType + 4]) - 0.153392;
			DemodSeq[i*ModulationType + 7] = fabs(DemodSeq[i*ModulationType + 5]) - 0.153392;

		}
		break;

	default:
		std::cout << "Invalid modulation type: 2:QPSK, 4:16-QAM, 6:64-QAM! \n";
		break;
	}
}
void CModulate::BPSKModulation(int8_t* SourceSeq)
{
	unsigned long i;
	for (i = 0; i < SymbolLen; ++i)
	{
		BPSKModSeq[i] = 2 * SourceSeq[i] - 1;
	}
}


void CModulate::BPSKDemodulation(float *Received_Seq)
{
	unsigned long i;
	for (i = 0; i < SymbolLen; ++i)
	{
		DemodSeq[i] = Received_Seq[i];
	}
}

ModStatistic CModulate::ModCalErr(float *demodseq, int8_t *encodeseq)
{
	ModStatistic ModTest;
	CLDPC *ldpc;
	CChannel *channel;
	Parameter_Simulation p_simulation;
	ofstream eeout;
	ReadProfile(&p_simulation);
	ModTest.ErrorBits = 0;
	ModTest.ErrorSymbol = 0;
	ModTest.ErrorFrame = 0;
	unsigned long errorBits = 0;
	unsigned long errorSymbol = 0;

	int errorsymbolflag = 0;
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	int m = 0;
	float* demodout;
	demodout = (float*)malloc(sizeof(float)*BitsOverChannel * 32);

	for (i = 0; i < 32; ++i)//32 frame serial comparison
	{
		errorBits = 0;
		errorSymbol = 0;
		//Hard decision before decoding, only look at the Information bit part
		for (j = 0; j < (NmoinsK - _ShortenBits) / ModulationType; ++j)
		{
			for (m = 0; m < ModulationType; m++)//Determine modulationtype bits
			{
				errorsymbolflag = 0;
				demodout[i*BitsOverChannel + ModulationType * j + m] = (demodseq[i*BitsOverChannel + ModulationType * j + m] > 0) ? 1 : 0;
				
				if (demodout[i*BitsOverChannel + ModulationType * j + m] != encodeseq[i*BitsOverChannel + ModulationType * j + m])//encodeseq places outputbits size as bitsoverchannel
				//if (demodout[i*BitsOverChannel + ModulationType * j + m] != encodeseq[i*(NmoinsK - _ShortenBits) + ModulationType * j + m])//encodeseq puts input bits size as NmoinsK - _ShortenBits
				{
					errorBits++;
					errorsymbolflag++;
				}
			}
			if (errorsymbolflag > 0)
			{
				errorSymbol++;//If there is an error in the modulationtype bit, the error sign is +1
			}
		}
		//One frame statistics
		if (errorBits > 0)
		{
			ModTest.ErrorBits += errorBits;
			ModTest.ErrorSymbol += errorSymbol;
			ModTest.ErrorFrame++;

		}
	}


	//int* testout;
	//int* testout1;
	//testout = (int*)malloc(sizeof(int)*BitsOverChannel * 32);
	//testout1 = (int*)malloc(sizeof(int)*BitsOverChannel * 32);
	//eeout.open("input.txt", std::ios::app);
	//for (i = 0; i < 32; ++i)
	//{
	//	for (j = 0; j < NmoinsK; j++)
	//	{
	//		testout[l] = (int)demodseq[i*NmoinsK + j];

	//		l++;
	//	}
	//	for (j = 0; j < BitsOverChannel; j++)
	//	{
	//		testout1[k] = (float)encodeseq[i*BitsOverChannel + j];
	//		k++;
	//	}
	//}


	//	eeout << "input=[ ";//The value after passing through the awgn channel, information+check
	//	for (j = 0; j < l; j++)

	//	{
	//		eeout << testout[j] << "\t";

	//	}
	//	eeout << "];";
	//	eeout << endl;

	//	eeout << "output=[ ";//The value after passing through the awgn channel, information+check
	//	for (j = 0; j < k; j++)

	//	{
	//		eeout << testout1[j] << "\t";

	//	}
	//	eeout << "];";
	//	eeout << endl;
	//	eeout << endl;
	//	eeout.close();
	//	getchar();
	//	
	//
	//free(testout);
	//free(testout1);


    free(demodout);
	return ModTest;
}
