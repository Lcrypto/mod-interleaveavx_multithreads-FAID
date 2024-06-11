#ifndef CTOOL_H
#define CTOOL_H

extern int MAX_THREADS; // The Maximum Threads
#define CONTINUE_SEED 0 // Whether to resume simulation from a breakpoint, i.e. use the seed exported from Temp.txt
#define REGULAR_COL_WEIGHT 3 // Uniformly weighted column portion of the matrix, used for DTBF

#include<cstdint>
#include<fstream>
#include<iostream>
#include<cstdlib>
#include<string>
#include<xmmintrin.h>
#include<smmintrin.h>
#include<immintrin.h>
#include"poisx_memalign.h"
#include<cerrno>
using namespace std;
/* 
	Struct: Parameter_Simulation

*/
struct Parameter_Simulation
{
	float snr_start;//Start SNR
	float snr_pass;// SNR Step
	float snr_end;// End SNR
	float scale;// Quantize Scale
	int decode_method; // 0: NMS, 1: OMS
	int Max_Iteration;//Maximum Iteration
	int mod_type;//Modulation Type
	int interleavemod_type;
	int Factor_1;//Factor_1
	int Factor_2;//Factor_2
	int nb_frames;//Parallel Decoding 
	char *fileName;
	int Z;// Circulant size
	int ce;//Whether to collect error information
};
//extern const unsigned short PosNoeudsVariable[9240];
//int errorFrameAnalyzer(int8_t * encodedBits, int8_t * detectedBits, int inForlength, int codeLength, int frames);
void uchar_transpose_avx(__m256i *src, __m256i *dst, int n);// Interlever
void uchar_itranspose_avx(__m256i *src, __m256i *dst, int n);//DeinterweaveDeinterleaving


void *vec_malloc(uint32_t size);
void ReadProfile(Parameter_Simulation * p);

extern int collectflag;

#endif // !CTOOL_H
