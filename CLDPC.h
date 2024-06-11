#ifndef CLDPC_H
#define CLDPC_H

#include "./Constants/Constants_SSE.h"
#include "CChannel.h"
#include "CTool.h"
#include "mkl.h"
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <emmintrin.h>
#include <immintrin.h>
#include <iomanip>
#include <iostream>
#include <smmintrin.h>
#include <xmmintrin.h>

// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html

#define TYPE __m256i

#define VECTOR_LOAD(ptr) (_mm256_load_si256(ptr))
#define VECTOR_UNCACHED_LOAD(ptr) (_mm256_stream_load_si256(ptr))
#define VECTOR_STORE(ptr, v) (_mm256_store_si256(ptr, v))
#define VECTOR_ADD(a, b) (_mm256_adds_epi8(a, b))
#define VECTOR_SUB(a, b) (_mm256_subs_epi8(a, b))
#define VECTOR_ABS(a) (_mm256_abs_epi8(a))
#define VECTOR_MAX(a, b) (_mm256_max_epi8(a, b))
#define VECTOR_MIN(a, b) (_mm256_min_epi8(a, b))
#define VECTOR_XOR(a, b) (_mm256_xor_si256(a, b))
#define VECTOR_OR(a, b) (_mm256_or_si256(a, b))
#define VECTOR_AND(a, b) (_mm256_and_si256(a, b))
#define VECTOR_ANDNOT(a, b) (_mm256_andnot_si256(a, b)) // (!a) & b
#define VECTOR_MIN_1(a, min1) (_mm256_min_epi8(a, min1))
#define VECTOR_SIGN(a, b) (_mm256_sign_epi8(a, b))

/*
_mm256_sign_epi8(a,b)
operation:
FOR j := 0 to 31
i := j*8
IF b[i+7:i] < 0
dst[i+7:i] := NEG(a[i+7:i])
ELSE IF b[i+7:i] = 0
dst[i+7:i] := 0
ELSE
dst[i+7:i] := a[i+7:i]
FI
ENDFOR
dst[MAX:256] := 0
*/
#define VECTOR_invSIGN2(val, sig) (_mm256_sign_epi8(val, sig)) // Inverse sign
#define VECTOR_EQUAL(a, b) (_mm256_cmpeq_epi8(a, b)) // Compare equal
#define VECTOR_ZERO (_mm256_setzero_si256())
#define VECTOR_SET1(a) (_mm256_set1_epi8(a))
#define VECTOR_SET2(b) (_mm256_set1_epi16(b))
#define VECTOR_PACK(hi, lo) (_mm256_packs_epi16(lo, hi)) // Note the differenece
#define VECTOR_UNPACK_HIGH(a) (_mm256_unpackhi_epi8(a, VECTOR_ZERO))

#define VECTOR_UNPACK_LOW(a) (_mm256_unpacklo_epi8(a, VECTOR_ZERO))
#define VECTOR_MUL(a, b) (_mm256_mullo_epi16(a, b))
#define VECTOR_DIV(a, b) (_mm256_div_epi16(a, b))
#define VECTOR_DIV32(a) (_mm256_srli_epi16(a, 5))
#define VECTOR_SUB_AND_SATURATE_VAR_8bits(a, b, mina) (VECTOR_MAX(VECTOR_SUB(a, b), mina))

// search the second min value
#define VECTOR_MIN_2(val, old_min1, min2) (VECTOR_MIN(min2, VECTOR_MAX(old_min1, val)))

#define VECTOR_GET_SIGN_BIT(a, b) (_mm256_and_si256(a, b))

#define VECTOR_SBU(a, b) (_mm256_subs_epu8(a, b))
#define VECTOR_ADD_AND_SATURATE_VAR_8bits(a, b, mina) (_mm256_max_epi8(_mm256_adds_epi8(a, b), mina))

#define VECTOR_SIGN(a, b) (_mm256_sign_epi8(a, b)) // if b < 0, return -a; if b == 0, return 0; if b > 0, return a
// MASK: returns __mmask32
#define VECTOR_GT_MASK(a, b) (_mm256_cmpgt_epi8_mask(a, b))
#define VECTOR_GE_MASK(a, b) (_mm256_cmpge_epi8_mask(a, b))
#define VECTOR_LT_MASK(a, b) (_mm256_cmplt_epi8_mask(a, b))
#define VECTOR_LE_MASK(a, b) (_mm256_cmple_epi8_mask(a, b))
#define VECTOR_EQ_MASK(a, b) (_mm256_cmpeq_epi8_mask(a, b))
#define VECTOR_ADD_MASK(m, a, b) (_mm256_mask_adds_epi8(a, m, a, b))
#define VECTOR_SUB_MASK(m, a, b) (_mm256_mask_subs_epi8(a, m, a, b))
#define VECTOR_NEQ_MASK(a, b) (_mm256_cmpneq_epi8_mask(a, b)) // == : 0; != : 1
#define VECTOR_MIN_MASK(src, m, a, b) _mm256_mask_min_epi8(src, m, a, b)
#define VECTOR_MAX_MASK(src, m, a, b) _mm256_mask_max_epi8(src, m, a, b)

#define VECTOR_MIN_ADD_MASK(m, a, b, min_var) VECTOR_ADD_MASK(m, VECTOR_MIN_MASK(a, m, a, min_var), b)
#define VECTOR_MIN_SUB_MASK(m, a, b, min_var) VECTOR_SUB_MASK(m, VECTOR_MIN_MASK(a, m, a, min_var), b)

#define VECTOR_GTU_MASK(a, b) (_mm256_cmpgt_epu8_mask(a, b))
#define VECTOR_LTU_MASK(a, b) (_mm256_cmplt_epu8_mask(a, b))
#define VECTOR_ADDU_MASK(m, a, b) (_mm256_mask_adds_epu8(a, m, a, b))

// if m[i] == 1: ret[i] = a[i]; else: ret[i] = src[i]
#define VECTOR_MOV_MASK(src, m, a) (_mm256_mask_mov_epi8(src, m, a))

#define PETIT 1
#define MANUAL_PREFETCH 1

using namespace std;
// For data statistic
struct Statistic {
    unsigned long ErrorFrame;
    unsigned long ErrorBits;
    unsigned long LT3ErrBitFrame;
    char padding[48]; // padding data no use
};
// The LDPC class
class CLDPC {
public:
    double m_Rate; // code ratio
    int8_t* inputBits; // Information bits
    int* errorbitblock; //The block where the error bit is located
    int* errorbitindex; //The location of the block where the error bit is located
    int* errorcheckblock; //The block where the error check is located
    int* errorcheckindex; //The location of the block where the error check is located
    float* errorfloat; //The value of the error bit after passing through the channel
    int* errorchar; //The quantized value of the error bit
    int8_t* outputBits; // Encoder ouput removes the encoding information after puncturing and shortening
    __m256i* var_nodes; // AVX decoder for store LLR plus likelihood ratios for puncture and shorten
    __m256i** p_vn_adr; // The pointer for LLR
    __m256i* var_msgs; // Lmn
    __m256i* encoder_varnodes; // for parallel encoding
    int8_t* decodedBits; // hard decision bits, sequential storage of each frame
    int8_t* fixInput; // fixed LLR for the length is N*m_frames; remove the original log-likelihood ratio of puncture and shorten
    int8_t* VN_weight_; // length: _NoVar

    int nb_iteration; // Maximum iterations
    int m_M; // check
    int m_N; // code length
    int m_K; // Information bits
    int m_PunLen; // Puncture
    int m_ShortenLen; // Shorten
    int m_frame; // no of frames for parallel
    char padding[28]; // For Padding no use
public:
    CLDPC(); //
    ~CLDPC();
    void GenMsgSeq(); // random information bits generator

    void FakeEncoder(); // the all zero encoder

    void Encode(); // Random Encoder

    void Decode(); // AVX Decoder
    void Decode1(); // AVX Decoder
    void Decode_OMS(); // OMS Decoder
    void Decode_FAID(); // FAID Decoder
    int Decode_OMSBF(); // OMS+BF Decoder
    int Decode_OMS_DTBF(); // OMS+DTBF Decoder
    void Decode_FAID_2B1C(); // FAID+2B1C Decoder

    void VN_weight_count(); // count the weight of each VN

    void float2LimitChar_6bit(int8_t* output, const float* input, float scale, int length); // Quantize 6 bit

    void float2LimitChar_5bit(int8_t* output, const float* input, float scale, int length); // Quantize 5 bit

    void float2LimitChar_4bit(int8_t* output, const float* input, float scale, int length); // Quantize 4 bit

    void float2LimitChar_3bit(int8_t* output, const float* input, float scale, int length); // Quantize 3 bit

    void float2LimitChar_2bit(int8_t* output, const float* input, float scale, int length); // Quantize 2 bit

    void float2LimitChar_1bit(int8_t* output, const float* input, float scale, int length); // Quantize 1 bit

    void Initial(int nb_frame, int MaxItertion); // Initial the variable
    Statistic CalculateErrors(float* bpskinput, int8_t* charinput, int collectflag);
    // Calculate the errors
};

extern int vSAT_NEG_MSG;
extern int vSAT_POS_MSG;

extern int vSAT_NEG_VAR;
extern int vSAT_POS_VAR;

extern int vSAT_NEG_LLR;
extern int vSAT_POS_LLR;

#endif // !CTRAME_H
