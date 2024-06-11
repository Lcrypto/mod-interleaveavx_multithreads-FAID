
#include "CLDPC.h"

#define OMS_MODE 0 // 0: simple; 1: selective
#define STOP_EARLY 1
#define EF_ELIMINATION 1 
#define FAID2_SIGN_BACKTRACK 1

const int offset = 0; // for simple OMS

// V2C_map[weight = 3, 6, 11, others][0~6, other]
const int8_t V2C_map_it1_[4][SAT_POS_MSG + 1] = {
    { 0, 0, 1, 2, 3, 3, 3, 3 }, // weight = 3
    { 0, 0, 1, 2, 3, 3, 3, 3 }, // weight = 6
    { 0, 0, 1, 2, 3, 3, 3, 3 }, // weight = 11
    { 0, 0, 1, 2, 3, 3, 3, 3 }, // weight = others
};
const int8_t V2C_map_it2_[4][SAT_POS_MSG + 1] = {
    { 0, 1, 1, 2, 3, 3, 3, 3 }, // weight = 3
    { 0, 1, 1, 2, 3, 3, 3, 3 }, // weight = 6
    { 0, 1, 1, 2, 3, 3, 3, 3 }, // weight = 11
    { 0, 1, 1, 2, 3, 3, 3, 3 }, // weight = others
};
const int8_t V2C_map_it3_[4][SAT_POS_MSG + 1] = {
    { 0, 1, 1, 2, 3, 3, 3, 3 }, // weight = 3
    { 0, 1, 1, 2, 3, 3, 3, 3 }, // weight = 6
    { 0, 1, 1, 2, 3, 3, 3, 3 }, // weight = 11
    { 0, 1, 1, 2, 3, 3, 3, 3 }, // weight = others
};
const int8_t V2C_map_it4_[4][SAT_POS_MSG + 1] = {
    { 0, 1, 1, 3, 3, 4, 4, 4 }, // weight = 3
    { 0, 1, 1, 3, 3, 4, 4, 4 }, // weight = 6
    { 0, 1, 1, 3, 3, 4, 4, 4 }, // weight = 11
    { 0, 1, 1, 3, 3, 4, 4, 4 }, // weight = others
};
const int8_t V2C_map_it5_[4][SAT_POS_MSG + 1] = {
    { 0, 1, 1, 3, 3, 3, 6, 6 }, // weight = 3
    { 0, 1, 1, 3, 3, 3, 6, 6 }, // weight = 6
    { 0, 1, 1, 3, 3, 3, 6, 6 }, // weight = 11
    { 0, 1, 1, 3, 3, 3, 6, 6 }, // weight = others
};
const int8_t V2C_map_it6_[4][SAT_POS_MSG + 1] = {
    { 0, 1, 1, 3, 3, 3, 7, 7 }, // weight = 3
    { 0, 1, 1, 3, 3, 3, 7, 7 }, // weight = 6
    { 0, 1, 1, 3, 3, 3, 7, 7 }, // weight = 11
    { 0, 1, 1, 3, 3, 3, 7, 7 }, // weight = others
};

// Tables for eliminating error floor
const int8_t V2C_map_it1_ef[4][SAT_POS_MSG + 1] = {
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = 3
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = 6
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = 11
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = others
};
const int8_t V2C_map_it2_ef[4][SAT_POS_MSG + 1] = {
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = 3
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = 6
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = 11
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = others
};
const int8_t V2C_map_it3_ef[4][SAT_POS_MSG + 1] = {
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = 3
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = 6
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = 11
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = others
};
const int8_t V2C_map_it4_ef[4][SAT_POS_MSG + 1] = {
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = 3
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = 6
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = 11
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = others
};
const int8_t V2C_map_it5_ef[4][SAT_POS_MSG + 1] = {
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = 3
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = 6
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = 11
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = others
};
const int8_t V2C_map_it6_ef[4][SAT_POS_MSG + 1] = {
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = 3
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = 6
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = 11
    { 2, 3, 3, 4, 5, 6, 6, 7 }, // weight = others
};

const int _delta = 1; // Th is the descending step size. The smaller the step size, the better the performance and the higher the complexity.
const int _L0 = 100; // The maximum threshold has the maximum number of iterations. The larger the threshold, the better the performance and the higher the complexity.
const int _L1 = 0; // The maximum threshold has the maximum number of iterations. The larger the threshold, the better the performance and the higher the complexity.
const int _alpha = 1;

/**
 * @brief OMS decoder
 *
 */
void CLDPC::Decode_FAID_2B1C()
{
    Parameter_Simulation p_simulation;
    ReadProfile(&p_simulation);

    int8_t* ptr;
    int8_t* pp;
    const TYPE factor_1 = VECTOR_SET1(p_simulation.Factor_1);
    const TYPE factor_2 = VECTOR_SET1(p_simulation.Factor_2);
    const TYPE zero = VECTOR_ZERO;
    const TYPE ones = VECTOR_SET1(1);
    const TYPE twos = VECTOR_SET1(2);
    // #define SAT_NEG_MSG(-(0x0001 << (NB_BITS_MESSAGES - 1)) + 1) -31(6bit quantization)
    const TYPE minMesg = VECTOR_SET1(SAT_NEG_MSG);

    // When the number of remaining iterations is only floor_iter_thresh, if the number of decoding checksum errors is less than floor_err_count, change the offset setting value
#if EF_ELIMINATION == 0
    const int8_t floor_err_count = 0; // 0~127, the conditions for enabling the leveling algorithm
    const int floor_iter_thresh = -1; // 0~MaxIteration-1, the conditions for enabling the leveling algorithm
#endif
#if EF_ELIMINATION == 1
    const int8_t floor_err_count = 50; // 0~127, the conditions for enabling the leveling algorithm
    const int floor_iter_thresh = 6; // 0~MaxIteration-1, the conditions for enabling the leveling algorithm?The number of iterations to eliminate flat layers = this value + 1
#endif
#if EF_ELIMINATION == 2
    const int8_t floor_err_count = 20; // Experience value <100, wrong number of checksums, conditions for enabling the leveling algorithm
    const int floor_iter_thresh = 6; // 0~MaxIteration-1, the conditions for enabling the leveling algorithm?The number of iterations to eliminate flat layers = this value + 1
#endif
    __mmask32 l_mask_eq = 0; // == : 0; != : 1, save the last verification result, initially all wrong
    __mmask32 l_m_error_sum = 0; // Is the last error count less than floor_err_count? If yes, it is 1
    __mmask32 l_checksum_[_NoCheck] = { 0 }; //== : 0; != : 1, initially both are correct
    __mmask32 era_[_NoVar] = { 0 }; // Record whether the information output by the current VN has been erased (only valid for EF_ELIMINATION=2)
    const int _maxBFiter = 10; // BF Iterations
    size_t i, j, z;

    /* Lmn Initialization */
    for (i = 0; i < MESSAGE; ++i) {
        var_msgs[i] = zero;
    }
    // The information part is interleaved. The demodulated initial likelihood ratio is arranged in a format that is only stored in parallel, not really interleaved.
    // The first codeword has the first log-likelihood ratio, the second codeword has the first, ... the 32nd codeword has the first, the first codeword has the second, ...
    if ((NmoinsK - _PunctureBits - _ShortenBits) % 32 == 0) {
        uchar_transpose_avx(
            (TYPE*)fixInput, (TYPE*)(var_nodes + _PunctureBits), (NmoinsK - _PunctureBits - _ShortenBits));
    } else {
        ptr = (int8_t*)(var_nodes + _PunctureBits); // The pointer starts after the puncture, which means the puncture is at the front.
        for (i = 0; i < (NmoinsK - _PunctureBits - _ShortenBits); ++i) {
            for (j = 0; j < 32; ++j) {
                // Original likelihood ratio without puncture and shortening
                ptr[i * 32 + j] = fixInput[j * (NmoinsK - _PunctureBits - _ShortenBits) + i];
            }
        }
    }
    // The check part
    pp = fixInput + (NmoinsK - _PunctureBits - _ShortenBits) * 32;
    // pp = fixInput + (NmoinsK - _PunctureBits - _ShortenBits);
    if (_NoCheck % 32 == 0) {
        uchar_transpose_avx((TYPE*)(pp), (TYPE*)(var_nodes + NmoinsK), _NoCheck);
    } else {
        ptr = (int8_t*)(var_nodes + NmoinsK);
        for (i = 0; i < _NoCheck; ++i) {
            for (j = 0; j < 32; ++j) {
                ptr[i * 32 + j] = pp[j * _NoCheck + i];
            }
        }
    }

    // // puncture  llr=0
    // for (i = 0; i < _PunctureBits; ++i) {
    //     var_nodes[i] = zero;
    // }

    // // shorten last llr=minMesg
    // for (i = NmoinsK - _ShortenBits; i < NmoinsK; ++i) {
    //     var_nodes[i] = minMesg;
    // }

    for (i = 0; i < 384; ++i) {
        var_nodes[_NoVar-1-i] = zero;
    }

    // Memory Channel LLR
    TYPE llr_ch[_NoVar];
    for (i = 0; i < _NoVar; i++) {
        llr_ch[i] = var_nodes[i];
    }

    int nombre_iterations = nb_iteration;
    // Fixed iterations
    while (nombre_iterations--) { // An iteration starts, without dynamic termination
        TYPE* p_msg1r = var_msgs; // read Lmn, C2V
        TYPE* p_msg1w = var_msgs; // write Lmn
#if PETIT == 1
        TYPE** p_indice_nod1 = p_vn_adr; // p_vn_adr: Stores the En,read, LLR
        TYPE** p_indice_nod2 = p_vn_adr; // write

#else
        const unsigned short* p_indice_nod1 = PosNoeudsVariable;
        const unsigned short* p_indice_nod2 = PosNoeudsVariable;
#endif

        const TYPE min_var = VECTOR_SET1(vSAT_NEG_VAR); // semimimmum value 8bit -127
        const TYPE max_msg = VECTOR_SET1(vSAT_POS_MSG); // maximum value 6bit +31
        const TYPE max_var = VECTOR_SET1(vSAT_POS_VAR); // maximum value 8 bit 127

        int cn_count = 0;

#if STOP_EARLY
        const unsigned short* pCN = PosNoeudsVariable;
        const unsigned short* pCN2 = PosNoeudsVariable;
        __mmask32 mask_sum; // == : 0; != : 1, initially both are correct
        __m256i error_sum = zero; // Add 1 to every wrong number, up to 127
        for (i = 0; i < DEG_1_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_1; j++) {
                // > 0 ? (== 0 ? ch > 0)
                mask_sum ^= VECTOR_GT_MASK(var_nodes[*pCN], zero);
                // | (VECTOR_EQ_MASK(var_nodes[*pCN], zero) & VECTOR_GT_MASK(llr_ch[*pCN], zero));
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADD_MASK(mask_sum, error_sum, ones);
        }
#if NB_DEGRES >= 2
        for (i = 0; i < DEG_2_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_2; j++) {
                mask_sum ^= VECTOR_GT_MASK(var_nodes[*pCN], zero);
                // | (VECTOR_EQ_MASK(var_nodes[*pCN], zero) & VECTOR_GT_MASK(llr_ch[*pCN], zero));
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADD_MASK(mask_sum, error_sum, ones);
        }
#endif
#if NB_DEGRES >= 3
        for (i = 0; i < DEG_3_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_3; j++) {
                mask_sum ^= VECTOR_GT_MASK(var_nodes[*pCN], zero);
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADD_MASK(mask_sum, error_sum, ones);
        }
#endif
#if NB_DEGRES >= 4
        for (i = 0; i < DEG_4_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_4; j++) {
                mask_sum ^= VECTOR_GT_MASK(var_nodes[*pCN], zero);
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADD_MASK(mask_sum, error_sum, ones);
        }
#endif
#if NB_DEGRES >= 5
        for (i = 0; i < DEG_5_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_5; j++) {
                mask_sum ^= VECTOR_GT_MASK(var_nodes[*pCN], zero);
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADD_MASK(mask_sum, error_sum, ones);
        }
#endif
#if NB_DEGRES >= 6
        for (i = 0; i < DEG_6_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_6; j++) {
                mask_sum ^= VECTOR_GT_MASK(var_nodes[*pCN], zero);
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADD_MASK(mask_sum, error_sum, ones);
        }
#endif
#if NB_DEGRES >= 7
        for (i = 0; i < DEG_7_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_7; j++) {
                mask_sum ^= VECTOR_GT_MASK(var_nodes[*pCN], zero);
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADD_MASK(mask_sum, error_sum, ones);
        }
#endif
#if NB_DEGRES >= 8
        for (i = 0; i < DEG_8_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_8; j++) {
                mask_sum ^= VECTOR_GT_MASK(var_nodes[*pCN], zero);
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADD_MASK(mask_sum, error_sum, ones);
        }
#endif
#if NB_DEGRES >= 9
        for (i = 0; i < DEG_9_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_9; j++) {
                mask_sum ^= VECTOR_GT_MASK(var_nodes[*pCN], zero);
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADD_MASK(mask_sum, error_sum, ones);
        }
#endif
#if NB_DEGRES >= 10
        for (i = 0; i < DEG_10_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_10; j++) {
                mask_sum ^= VECTOR_GT_MASK(var_nodes[*pCN], zero);
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADD_MASK(mask_sum, error_sum, ones);
        }
#endif
#if NB_DEGRES >= 11
        for (i = 0; i < DEG_11_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_11; j++) {
                mask_sum ^= VECTOR_GT_MASK(var_nodes[*pCN], zero);
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADD_MASK(mask_sum, error_sum, ones);
        }
#endif
#if NB_DEGRES >= 12
        for (i = 0; i < DEG_12_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_12; j++) {
                mask_sum ^= VECTOR_GT_MASK(var_nodes[*pCN], zero);
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADD_MASK(mask_sum, error_sum, ones);
        }
#endif
#if NB_DEGRES >= 13
        for (i = 0; i < DEG_13_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_13; j++) {
                mask_sum ^= VECTOR_GT_MASK(var_nodes[*pCN], zero);
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADD_MASK(mask_sum, error_sum, ones);
        }
#endif
#if NB_DEGRES >= 14
        for (i = 0; i < DEG_14_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_14; j++) {
                mask_sum ^= VECTOR_GT_MASK(var_nodes[*pCN], zero);
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADD_MASK(mask_sum, error_sum, ones);
        }
#endif
#if NB_DEGRES >= 15
        for (i = 0; i < DEG_15_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_15; j++) {
                mask_sum ^= VECTOR_GT_MASK(var_nodes[*pCN], zero);
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADD_MASK(mask_sum, error_sum, ones);
        }
#endif
#if NB_DEGRES >= 16
        for (i = 0; i < DEG_16_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_16; j++) {
                mask_sum ^= VECTOR_GT_MASK(var_nodes[*pCN], zero);
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADD_MASK(mask_sum, error_sum, ones);
        }
#endif
#if NB_DEGRES >= 17
        for (i = 0; i < DEG_17_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_17; j++) {
                mask_sum ^= VECTOR_GT_MASK(var_nodes[*pCN], zero);
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADD_MASK(mask_sum, error_sum, ones);
        }
#endif
#if NB_DEGRES >= 18
        for (i = 0; i < DEG_18_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_18; j++) {
                mask_sum ^= VECTOR_GT_MASK(var_nodes[*pCN], zero);
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADD_MASK(mask_sum, error_sum, ones);
        }
#endif
#if NB_DEGRES >= 19
        for (i = 0; i < DEG_19_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_19; j++) {
                mask_sum ^= VECTOR_GT_MASK(var_nodes[*pCN], zero);
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADD_MASK(mask_sum, error_sum, ones);
        }
#endif
#if NB_DEGRES >= 20
        for (i = 0; i < DEG_20_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_20; j++) {
                mask_sum ^= VECTOR_GT_MASK(var_nodes[*pCN], zero);
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADD_MASK(mask_sum, error_sum, ones);
        }
#endif
        if (VECTOR_GT_MASK(error_sum, zero) == 0) { // all the frames are right
            break;
        }
        l_m_error_sum = VECTOR_LT_MASK(error_sum, VECTOR_SET1(floor_err_count));
#endif

        cn_count = 0;
#if EF_ELIMINATION == 2
        // At the beginning of each iteration, the erase flag is set to 0.
        for (i = 0; i < _NoVar; i++) {
            era_[i] = 0;
        }
#endif
        // DEG_1_COMPUTATIONS: The number of the rows of degree DEG_1
        /**************************************DEG_1****************************************************************/
        for (i = 0; i < DEG_1_COMPUTATIONS; i++) {
            TYPE tab_vContr[DEG_1]; // DEG1
            TYPE temp_vContr[DEG_1]; // V2C before FAID mapping
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR); // 8bit maximum value  +127
            TYPE min2 = min1; // min2=min1=127
            TYPE _sign[DEG_1];

#if (DEG_1 & 0x01) == 1
            const unsigned char sign8 = 0x80; // sign8 =128=1000000B
            // isign =12*16=11000000B, the reason why X100,0000 is used is that if its sign is 0, then
            // _mm_sign_epi8 operation meets the requirements.
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80;
            const unsigned char isign8b = 0x40; // isign8b=01000000B=64
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

#if PETIT == 1
#if MANUAL_PREFETCH == 1
            _mm_prefetch((const char*)(p_indice_nod1[DEG_1]), _MM_HINT_T0);
            _mm_prefetch((const char*)(&p_msg1r[DEG_1]), _MM_HINT_T0);
#endif
#endif

#pragma unroll(DEG_1)
            // Kernal 2 in algorithm
            for (j = 0; j < DEG_1; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1); // vNoeud:in algorithm:En
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                unsigned short col = PosNoeudsVariable[p_indice_nod1 - p_vn_adr];
                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                vContr = VECTOR_MIN(vContr, max_var);

#if FAID2_SIGN_BACKTRACK == 1
                TYPE cSign = VECTOR_GET_SIGN_BIT(VECTOR_ADD_MASK(VECTOR_EQ_MASK(vContr, zero), vContr, vNoeud), msign8);
#else
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
#endif
                _sign[j] = cSign;
                sign = VECTOR_XOR(sign, cSign); // 0: pos, 1: neg
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                temp_vContr[j] = vContr;

                int idx1 = 0;
                switch (VN_weight_[col]) {
                case 3:
                    idx1 = 0;
                    break;
                case 6:
                    idx1 = 1;
                    break;
                case 11:
                    idx1 = 2;
                    break;
                default:
                    idx1 = 3;
                    break;
                }
                TYPE tmp = zero;
                TYPE tmp1 = zero; // Used to eliminate temporary storage value during leveling
                TYPE tmp2 = zero; // Used to eliminate temporary storage value during leveling
                // Step-by-step mapping
                for (int idx2 = 0; idx2 < SAT_POS_MSG + 1; idx2++) {
                    __mmask32 mask = VECTOR_EQ_MASK(vAbs, VECTOR_SET1(idx2)); // Filter out the part equal to idx2 and map it
#if EF_ELIMINATION >= 1
                    __mmask32 mask_iter = 0;
                    if (nombre_iterations <= floor_iter_thresh) {
                        mask_iter = UINT32_MAX;
                    } else {
                        mask_iter = 0;
                    }

                    __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 2:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 3:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 4:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 5:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    default:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    }
#else
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        break;
                    case 2:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        break;
                    case 3:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        break;
                    case 4:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        break;
                    case 5:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        break;
                    default:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        break;
                    }
#endif
                }
                // Handle overflow, because the magnitude of abs may overflow the bit width
                __mmask32 mask = VECTOR_GE_MASK(vAbs, VECTOR_SET1(SAT_POS_MSG + 1));
#if EF_ELIMINATION >= 1
                __mmask32 mask_iter = 0;
                if (nombre_iterations <= floor_iter_thresh) {
                    mask_iter = UINT32_MAX;
                } else {
                    mask_iter = 0;
                }

                __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 2:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 3:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 4:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 5:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                default:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                }
#else
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    break;
                case 2:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    break;
                case 3:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    break;
                case 4:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    break;
                case 5:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    break;
                default:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    break;
                }
#endif

                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(tmp, min1);
                min2 = VECTOR_MIN_2(tmp, vTemp, min2); // the second minimum value of the Lnm

                tab_vContr[j] = tmp; // Just need to store abs, sign already exists in _sign array
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

            // WE SATURATE DIRECTLY IN MSG FORMAT
#if OMS_MODE == 0 // simple OMS
            min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2 = VECTOR_MIN(min2, max_msg);
            TYPE cste_1 = VECTOR_MIN(VECTOR_SUB(min2, VECTOR_SET1(offset)), max_msg);
            TYPE cste_2 = VECTOR_MIN(VECTOR_SUB(min1, VECTOR_SET1(offset)), max_msg);
#elif OMS_MODE == 1 // selective OMS

            TYPE min1_offed, min2_offed;
            min1_offed = min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2_offed = min2 = VECTOR_MIN(min2, max_msg);

            if (nombre_iterations <= floor_iter_thresh) {
                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min1, factor_2);
                min1_offed = VECTOR_ADD_MASK(msk_lt6_1, min1, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min1, factor_1);
                min1_offed = VECTOR_ADD_MASK(msk_le1_1, min1_offed, ones); // min_offed = min_offed + 1 = min + 2

                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min2, factor_2);
                min2_offed = VECTOR_ADD_MASK(msk_lt6_2, min2, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min2, factor_1);
                min2_offed = VECTOR_ADD_MASK(msk_le1_2, min2_offed, ones); // min_offed = min_offed + 1 = min + 2
            }
            TYPE cste_1 = VECTOR_MIN(min2_offed, max_msg); // Limiting quantization width
            TYPE cste_2 = VECTOR_MIN(min1_offed, max_msg);

#endif

#if PETIT == 1
#if MANUAL_PREFETCH == 1
            for (j = 0; j < DEG_1; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_1]), _MM_HINT_T0);
#endif
#endif

#if (DEG_1 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8); // xor operation
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif
            // immediately update En the operations in Algorithm from line 18 to 21
#pragma unroll(DEG_1)
            for (j = 0; j < DEG_1; j++) {
                TYPE vContr = tab_vContr[j]; // Lnm in the ith iteration
                TYPE vAbs = VECTOR_ABS(vContr); // Lmn_new test if saturates
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, _sign[j]); // sign is in Algorithm Line 15
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig); // compute the Lmn in the next iteration
                TYPE v2Sr
                    = VECTOR_ADD_AND_SATURATE_VAR_8bits(temp_vContr[j], v2St, min_var); // update the En in Algorithm 21
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St); // Store Lmn
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr); // Store En
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
#if PETIT == 1
#if MANUAL_PREFETCH == 1
            _mm_prefetch((const char*)(*p_indice_nod2), _MM_HINT_T0);
#endif
#endif
            cn_count++;
        }

        ////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////

#if NB_DEGRES >= 2
        /********************************************DEG_2******************************************************************/
        for (int i = 0; i < DEG_2_COMPUTATIONS; i++) {

#if (DEG_2 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_2];
            TYPE temp_vContr[DEG_2];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;
            TYPE _sign[DEG_2];

#pragma unroll(DEG_2)
            for (int j = 0; j < DEG_2; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                unsigned short col = PosNoeudsVariable[p_indice_nod1 - p_vn_adr];

                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                vContr = VECTOR_MIN(vContr, max_var);

#if FAID2_SIGN_BACKTRACK == 1
                TYPE cSign = VECTOR_GET_SIGN_BIT(VECTOR_ADD_MASK(VECTOR_EQ_MASK(vContr, zero), vContr, vNoeud), msign8);
#else
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
#endif
                _sign[j] = cSign;
                sign = VECTOR_XOR(sign, cSign); // 0: pos, 1: neg
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                temp_vContr[j] = vContr;

                int idx1 = 0;
                switch (VN_weight_[col]) {
                case 3:
                    idx1 = 0;
                    break;
                case 6:
                    idx1 = 1;
                    break;
                case 11:
                    idx1 = 2;
                    break;
                default:
                    idx1 = 3;
                    break;
                }
                TYPE tmp = zero;
                TYPE tmp1 = zero; // Used to eliminate temporary storage value during leveling
                TYPE tmp2 = zero; // Used to eliminate temporary storage value during leveling
                // Step-by-step mapping
                for (int idx2 = 0; idx2 < SAT_POS_MSG + 1; idx2++) {
                    __mmask32 mask = VECTOR_EQ_MASK(vAbs, VECTOR_SET1(idx2)); // Filter out the part equal to idx2 and map it
#if EF_ELIMINATION >= 1
                    __mmask32 mask_iter = 0;
                    if (nombre_iterations <= floor_iter_thresh) {
                        mask_iter = UINT32_MAX;
                    } else {
                        mask_iter = 0;
                    }

                    __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 2:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 3:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 4:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 5:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    default:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    }
#else
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        break;
                    case 2:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        break;
                    case 3:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        break;
                    case 4:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        break;
                    case 5:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        break;
                    default:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        break;
                    }
#endif
                }
                // Handle overflow, because the magnitude of abs may overflow the bit width
                __mmask32 mask = VECTOR_GE_MASK(vAbs, VECTOR_SET1(SAT_POS_MSG + 1));
#if EF_ELIMINATION >= 1
                __mmask32 mask_iter = 0;
                if (nombre_iterations <= floor_iter_thresh) {
                    mask_iter = UINT32_MAX;
                } else {
                    mask_iter = 0;
                }

                __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 2:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 3:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 4:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 5:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                default:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                }
#else
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    break;
                case 2:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    break;
                case 3:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    break;
                case 4:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    break;
                case 5:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    break;
                default:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    break;
                }
#endif

                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(tmp, min1);
                min2 = VECTOR_MIN_2(tmp, vTemp, min2); // the second minimum value of the Lnm

                tab_vContr[j] = tmp;
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_2; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_2]), _MM_HINT_T0);
#endif
            // WE SATURATE DIRECTLY IN MSG FORMAT
#if OMS_MODE == 0 // simple OMS
            min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2 = VECTOR_MIN(min2, max_msg);
            TYPE cste_1 = VECTOR_MIN(VECTOR_SUB(min2, VECTOR_SET1(offset)), max_msg);
            TYPE cste_2 = VECTOR_MIN(VECTOR_SUB(min1, VECTOR_SET1(offset)), max_msg);
#elif OMS_MODE == 1 // selective OMS

            TYPE min1_offed, min2_offed;
            min1_offed = min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2_offed = min2 = VECTOR_MIN(min2, max_msg);

            if (nombre_iterations <= floor_iter_thresh) {
                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min1, factor_2);
                min1_offed = VECTOR_ADD_MASK(msk_lt6_1, min1, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min1, factor_1);
                min1_offed = VECTOR_ADD_MASK(msk_le1_1, min1_offed, ones); // min_offed = min_offed + 1 = min + 2

                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min2, factor_2);
                min2_offed = VECTOR_ADD_MASK(msk_lt6_2, min2, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min2, factor_1);
                min2_offed = VECTOR_ADD_MASK(msk_le1_2, min2_offed, ones); // min_offed = min_offed + 1 = min + 2
            }
            TYPE cste_1 = VECTOR_MIN(min2_offed, max_msg);
            TYPE cste_2 = VECTOR_MIN(min1_offed, max_msg);

#endif

#if (DEG_2 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_2)
            for (int j = 0; j < DEG_2; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                                             // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE vSig = VECTOR_XOR(sign, _sign[j]);
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(temp_vContr[j], v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }

            cn_count++;
        }
#endif
#if NB_DEGRES >= 3
        /****************************************************DEG_3****************************************************************/
        for (int i = 0; i < DEG_3_COMPUTATIONS; i++) {

#if (DEG_3 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_3];
            TYPE temp_vContr[DEG_3];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;
            TYPE _sign[DEG_3];

#pragma unroll(DEG_3)
            for (int j = 0; j < DEG_3; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                unsigned short col = PosNoeudsVariable[p_indice_nod1 - p_vn_adr];

                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                vContr = VECTOR_MIN(vContr, max_var);

#if FAID2_SIGN_BACKTRACK == 1
                TYPE cSign = VECTOR_GET_SIGN_BIT(VECTOR_ADD_MASK(VECTOR_EQ_MASK(vContr, zero), vContr, vNoeud), msign8);
#else
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
#endif
                _sign[j] = cSign;
                sign = VECTOR_XOR(sign, cSign); // 0: pos, 1: neg
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                temp_vContr[j] = vContr;

                int idx1 = 0;
                switch (VN_weight_[col]) {
                case 3:
                    idx1 = 0;
                    break;
                case 6:
                    idx1 = 1;
                    break;
                case 11:
                    idx1 = 2;
                    break;
                default:
                    idx1 = 3;
                    break;
                }
                TYPE tmp = zero;
                TYPE tmp1 = zero; // Used to eliminate temporary storage value during leveling
                TYPE tmp2 = zero; // Used to eliminate temporary storage value during leveling
                // Step-by-step mapping
                for (int idx2 = 0; idx2 < SAT_POS_MSG + 1; idx2++) {
                    __mmask32 mask = VECTOR_EQ_MASK(vAbs, VECTOR_SET1(idx2)); // Filter out the part equal to idx2 and map it
#if EF_ELIMINATION >= 1
                    __mmask32 mask_iter = 0;
                    if (nombre_iterations <= floor_iter_thresh) {
                        mask_iter = UINT32_MAX;
                    } else {
                        mask_iter = 0;
                    }

                    __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 2:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 3:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 4:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 5:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    default:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    }
#else
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        break;
                    case 2:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        break;
                    case 3:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        break;
                    case 4:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        break;
                    case 5:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        break;
                    default:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        break;
                    }
#endif
                }
                // Handle overflow, because the magnitude of abs may overflow the bit width
                __mmask32 mask = VECTOR_GE_MASK(vAbs, VECTOR_SET1(SAT_POS_MSG + 1));
#if EF_ELIMINATION >= 1
                __mmask32 mask_iter = 0;
                if (nombre_iterations <= floor_iter_thresh) {
                    mask_iter = UINT32_MAX;
                } else {
                    mask_iter = 0;
                }

                __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 2:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 3:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 4:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 5:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                default:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                }
#else
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    break;
                case 2:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    break;
                case 3:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    break;
                case 4:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    break;
                case 5:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    break;
                default:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    break;
                }
#endif

                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(tmp, min1);
                min2 = VECTOR_MIN_2(tmp, vTemp, min2); // the second minimum value of the Lnm

                tab_vContr[j] = tmp;
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_3; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_3]), _MM_HINT_T0);
#endif
            // WE SATURATE DIRECTLY IN MSG FORMAT
#if OMS_MODE == 0 // simple OMS
            min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2 = VECTOR_MIN(min2, max_msg);
            TYPE cste_1 = VECTOR_MIN(VECTOR_SUB(min2, VECTOR_SET1(offset)), max_msg);
            TYPE cste_2 = VECTOR_MIN(VECTOR_SUB(min1, VECTOR_SET1(offset)), max_msg);
#elif OMS_MODE == 1 // selective OMS

            TYPE min1_offed, min2_offed;
            min1_offed = min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2_offed = min2 = VECTOR_MIN(min2, max_msg);

            if (nombre_iterations <= floor_iter_thresh) {
                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min1, factor_2);
                min1_offed = VECTOR_ADD_MASK(msk_lt6_1, min1, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min1, factor_1);
                min1_offed = VECTOR_ADD_MASK(msk_le1_1, min1_offed, ones); // min_offed = min_offed + 1 = min + 2

                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min2, factor_2);
                min2_offed = VECTOR_ADD_MASK(msk_lt6_2, min2, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min2, factor_1);
                min2_offed = VECTOR_ADD_MASK(msk_le1_2, min2_offed, ones); // min_offed = min_offed + 1 = min + 2
            }
            TYPE cste_1 = VECTOR_MIN(min2_offed, max_msg);
            TYPE cste_2 = VECTOR_MIN(min1_offed, max_msg);

#endif
#if (DEG_3 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_3)
            for (int j = 0; j < DEG_3; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, _sign[j]);
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(temp_vContr[j], v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }

            cn_count++;
        }
#endif
#if NB_DEGRES >= 4
        /*************************************DEG_4*************************************************/
        for (int i = 0; i < DEG_4_COMPUTATIONS; i++) {

#if (DEG_4 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_4];
            TYPE temp_vContr[DEG_4];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_4)
            for (int j = 0; j < DEG_4; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                unsigned short col = PosNoeudsVariable[p_indice_nod1 - p_vn_adr];

                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                vContr = VECTOR_MIN(vContr, max_var);

                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
                sign = VECTOR_XOR(sign, cSign); // 0: pos, 1: neg
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                temp_vContr[j] = vContr;

                int idx1 = 0;
                switch (VN_weight_[col]) {
                case 3:
                    idx1 = 0;
                    break;
                case 6:
                    idx1 = 1;
                    break;
                case 11:
                    idx1 = 2;
                    break;
                default:
                    idx1 = 3;
                    break;
                }
                TYPE tmp = zero;
                TYPE tmp1 = zero; // Used to eliminate temporary storage value during leveling
                TYPE tmp2 = zero; // Used to eliminate temporary storage value during leveling
                // Step-by-step mapping
                for (int idx2 = 0; idx2 < SAT_POS_MSG + 1; idx2++) {
                    __mmask32 mask = VECTOR_EQ_MASK(vAbs, VECTOR_SET1(idx2)); // Filter out the part equal to idx2 and map it
#if EF_ELIMINATION >= 1
                    __mmask32 mask_iter = 0;
                    if (nombre_iterations <= floor_iter_thresh) {
                        mask_iter = UINT32_MAX;
                    } else {
                        mask_iter = 0;
                    }

                    __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 2:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 3:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 4:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 5:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    default:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    }
#else
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        break;
                    case 2:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        break;
                    case 3:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        break;
                    case 4:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        break;
                    case 5:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        break;
                    default:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        break;
                    }
#endif
                }
                // Handle overflow, because the magnitude of abs may overflow the bit width
                __mmask32 mask = VECTOR_GE_MASK(vAbs, VECTOR_SET1(SAT_POS_MSG + 1));
#if EF_ELIMINATION >= 1
                __mmask32 mask_iter = 0;
                if (nombre_iterations <= floor_iter_thresh) {
                    mask_iter = UINT32_MAX;
                } else {
                    mask_iter = 0;
                }

                __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 2:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 3:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 4:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 5:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                default:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                }
#else
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    break;
                case 2:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    break;
                case 3:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    break;
                case 4:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    break;
                case 5:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    break;
                default:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    break;
                }
#endif

                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(tmp, min1);
                min2 = VECTOR_MIN_2(tmp, vTemp, min2); // the second minimum value of the Lnm

                tmp = VECTOR_ADD_MASK(VECTOR_EQ_MASK(vContr, zero), VECTOR_SIGN(tmp, vContr), tmp);
                tab_vContr[j] = tmp;
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_4; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_4]), _MM_HINT_T0);
#endif
            // WE SATURATE DIRECTLY IN MSG FORMAT
#if OMS_MODE == 0 // simple OMS
            min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2 = VECTOR_MIN(min2, max_msg);
            TYPE cste_1 = VECTOR_MIN(VECTOR_SUB(min2, VECTOR_SET1(offset)), max_msg);
            TYPE cste_2 = VECTOR_MIN(VECTOR_SUB(min1, VECTOR_SET1(offset)), max_msg);
#elif OMS_MODE == 1 // selective OMS

            TYPE min1_offed, min2_offed;
            min1_offed = min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2_offed = min2 = VECTOR_MIN(min2, max_msg);

            if (nombre_iterations <= floor_iter_thresh) {
                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min1, factor_2);
                min1_offed = VECTOR_ADD_MASK(msk_lt6_1, min1, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min1, factor_1);
                min1_offed = VECTOR_ADD_MASK(msk_le1_1, min1_offed, ones); // min_offed = min_offed + 1 = min + 2

                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min2, factor_2);
                min2_offed = VECTOR_ADD_MASK(msk_lt6_2, min2, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min2, factor_1);
                min2_offed = VECTOR_ADD_MASK(msk_le1_2, min2_offed, ones); // min_offed = min_offed + 1 = min + 2
            }
            TYPE cste_1 = VECTOR_MIN(min2_offed, max_msg);
            TYPE cste_2 = VECTOR_MIN(min1_offed, max_msg);

#endif
#if (DEG_4 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_4)
            for (int j = 0; j < DEG_4; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(temp_vContr[j], v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
            cn_count++;
        }
#endif
#if NB_DEGRES >= 5
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_5_COMPUTATIONS; i++) {

#if (DEG_5 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_5];
            TYPE temp_vContr[DEG_5];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_5)
            for (int j = 0; j < DEG_5; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                unsigned short col = PosNoeudsVariable[p_indice_nod1 - p_vn_adr];

                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                vContr = VECTOR_MIN(vContr, max_var);

                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
                sign = VECTOR_XOR(sign, cSign); // 0: pos, 1: neg
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                temp_vContr[j] = vContr;

                int idx1 = 0;
                switch (VN_weight_[col]) {
                case 3:
                    idx1 = 0;
                    break;
                case 6:
                    idx1 = 1;
                    break;
                case 11:
                    idx1 = 2;
                    break;
                default:
                    idx1 = 3;
                    break;
                }
                TYPE tmp = zero;
                TYPE tmp1 = zero; // Used to eliminate temporary storage value during leveling
                TYPE tmp2 = zero; // Used to eliminate temporary storage value during leveling
                // Step-by-step mapping
                for (int idx2 = 0; idx2 < SAT_POS_MSG + 1; idx2++) {
                    __mmask32 mask = VECTOR_EQ_MASK(vAbs, VECTOR_SET1(idx2)); // Filter out the part equal to idx2 and map it
#if EF_ELIMINATION >= 1
                    __mmask32 mask_iter = 0;
                    if (nombre_iterations <= floor_iter_thresh) {
                        mask_iter = UINT32_MAX;
                    } else {
                        mask_iter = 0;
                    }

                    __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 2:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 3:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 4:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 5:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    default:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    }
#else
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        break;
                    case 2:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        break;
                    case 3:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        break;
                    case 4:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        break;
                    case 5:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        break;
                    default:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        break;
                    }
#endif
                }
                // Handle overflow, because the magnitude of abs may overflow the bit width
                __mmask32 mask = VECTOR_GE_MASK(vAbs, VECTOR_SET1(SAT_POS_MSG + 1));
#if EF_ELIMINATION >= 1
                __mmask32 mask_iter = 0;
                if (nombre_iterations <= floor_iter_thresh) {
                    mask_iter = UINT32_MAX;
                } else {
                    mask_iter = 0;
                }

                __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 2:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 3:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 4:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 5:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                default:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                }
#else
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    break;
                case 2:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    break;
                case 3:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    break;
                case 4:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    break;
                case 5:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    break;
                default:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    break;
                }
#endif

                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(tmp, min1);
                min2 = VECTOR_MIN_2(tmp, vTemp, min2); // the second minimum value of the Lnm

                tmp = VECTOR_ADD_MASK(VECTOR_EQ_MASK(vContr, zero), VECTOR_SIGN(tmp, vContr), tmp);
                tab_vContr[j] = tmp;
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_5; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_5]), _MM_HINT_T0);
#endif
            // WE SATURATE DIRECTLY IN MSG FORMAT
#if OMS_MODE == 0 // simple OMS
            min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2 = VECTOR_MIN(min2, max_msg);
            TYPE cste_1 = VECTOR_MIN(VECTOR_SUB(min2, VECTOR_SET1(offset)), max_msg);
            TYPE cste_2 = VECTOR_MIN(VECTOR_SUB(min1, VECTOR_SET1(offset)), max_msg);
#elif OMS_MODE == 1 // selective OMS

            TYPE min1_offed, min2_offed;
            min1_offed = min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2_offed = min2 = VECTOR_MIN(min2, max_msg);

            if (nombre_iterations <= floor_iter_thresh) {
                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min1, factor_2);
                min1_offed = VECTOR_ADD_MASK(msk_lt6_1, min1, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min1, factor_1);
                min1_offed = VECTOR_ADD_MASK(msk_le1_1, min1_offed, ones); // min_offed = min_offed + 1 = min + 2

                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min2, factor_2);
                min2_offed = VECTOR_ADD_MASK(msk_lt6_2, min2, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min2, factor_1);
                min2_offed = VECTOR_ADD_MASK(msk_le1_2, min2_offed, ones); // min_offed = min_offed + 1 = min + 2
            }
            TYPE cste_1 = VECTOR_MIN(min2_offed, max_msg);
            TYPE cste_2 = VECTOR_MIN(min1_offed, max_msg);

#endif
#if (DEG_5 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_5)
            for (int j = 0; j < DEG_5; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(temp_vContr[j], v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
            cn_count++;
        }
#endif
#if NB_DEGRES >= 6
        /*************************************DEG_6*************************************************/
        for (int i = 0; i < DEG_6_COMPUTATIONS; i++) {

#if (DEG_6 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_6];
            TYPE temp_vContr[DEG_6];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_6)
            for (int j = 0; j < DEG_6; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                unsigned short col = PosNoeudsVariable[p_indice_nod1 - p_vn_adr];

                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                vContr = VECTOR_MIN(vContr, max_var);

                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
                sign = VECTOR_XOR(sign, cSign); // 0: pos, 1: neg
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                temp_vContr[j] = vContr;

                int idx1 = 0;
                switch (VN_weight_[col]) {
                case 3:
                    idx1 = 0;
                    break;
                case 6:
                    idx1 = 1;
                    break;
                case 11:
                    idx1 = 2;
                    break;
                default:
                    idx1 = 3;
                    break;
                }
                TYPE tmp = zero;
                TYPE tmp1 = zero; // Used to eliminate temporary storage value during leveling
                TYPE tmp2 = zero; // Used to eliminate temporary storage value during leveling
                // Step-by-step mapping
                for (int idx2 = 0; idx2 < SAT_POS_MSG + 1; idx2++) {
                    __mmask32 mask = VECTOR_EQ_MASK(vAbs, VECTOR_SET1(idx2)); // Filter out the part equal to idx2 and map it
#if EF_ELIMINATION >= 1
                    __mmask32 mask_iter = 0;
                    if (nombre_iterations <= floor_iter_thresh) {
                        mask_iter = UINT32_MAX;
                    } else {
                        mask_iter = 0;
                    }

                    __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 2:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 3:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 4:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 5:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    default:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    }
#else
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        break;
                    case 2:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        break;
                    case 3:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        break;
                    case 4:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        break;
                    case 5:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        break;
                    default:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        break;
                    }
#endif
                }
                // Handle overflow, because the magnitude of abs may overflow the bit width
                __mmask32 mask = VECTOR_GE_MASK(vAbs, VECTOR_SET1(SAT_POS_MSG + 1));
#if EF_ELIMINATION >= 1
                __mmask32 mask_iter = 0;
                if (nombre_iterations <= floor_iter_thresh) {
                    mask_iter = UINT32_MAX;
                } else {
                    mask_iter = 0;
                }

                __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 2:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 3:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 4:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 5:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                default:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                }
#else
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    break;
                case 2:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    break;
                case 3:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    break;
                case 4:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    break;
                case 5:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    break;
                default:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    break;
                }
#endif

                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(tmp, min1);
                min2 = VECTOR_MIN_2(tmp, vTemp, min2); // the second minimum value of the Lnm

                tmp = VECTOR_ADD_MASK(VECTOR_EQ_MASK(vContr, zero), VECTOR_SIGN(tmp, vContr), tmp);
                tab_vContr[j] = tmp;
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_6; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_6]), _MM_HINT_T0);
#endif
            // WE SATURATE DIRECTLY IN MSG FORMAT
#if OMS_MODE == 0 // simple OMS
            min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2 = VECTOR_MIN(min2, max_msg);
            TYPE cste_1 = VECTOR_MIN(VECTOR_SUB(min2, VECTOR_SET1(offset)), max_msg);
            TYPE cste_2 = VECTOR_MIN(VECTOR_SUB(min1, VECTOR_SET1(offset)), max_msg);
#elif OMS_MODE == 1 // selective OMS

            TYPE min1_offed, min2_offed;
            min1_offed = min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2_offed = min2 = VECTOR_MIN(min2, max_msg);

            if (nombre_iterations <= floor_iter_thresh) {
                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min1, factor_2);
                min1_offed = VECTOR_ADD_MASK(msk_lt6_1, min1, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min1, factor_1);
                min1_offed = VECTOR_ADD_MASK(msk_le1_1, min1_offed, ones); // min_offed = min_offed + 1 = min + 2

                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min2, factor_2);
                min2_offed = VECTOR_ADD_MASK(msk_lt6_2, min2, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min2, factor_1);
                min2_offed = VECTOR_ADD_MASK(msk_le1_2, min2_offed, ones); // min_offed = min_offed + 1 = min + 2
            }
            TYPE cste_1 = VECTOR_MIN(min2_offed, max_msg);
            TYPE cste_2 = VECTOR_MIN(min1_offed, max_msg);

#endif
#if (DEG_6 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_6)
            for (int j = 0; j < DEG_6; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(temp_vContr[j], v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
            cn_count++;
        }
#endif
#if NB_DEGRES >= 7
        /*************************************DEG_7*************************************************/
        for (int i = 0; i < DEG_7_COMPUTATIONS; i++) {

#if (DEG_7 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_7];
            TYPE temp_vContr[DEG_7];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_7)
            for (int j = 0; j < DEG_7; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                unsigned short col = PosNoeudsVariable[p_indice_nod1 - p_vn_adr];

                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                vContr = VECTOR_MIN(vContr, max_var);

                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
                sign = VECTOR_XOR(sign, cSign); // 0: pos, 1: neg
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                temp_vContr[j] = vContr;

                int idx1 = 0;
                switch (VN_weight_[col]) {
                case 3:
                    idx1 = 0;
                    break;
                case 6:
                    idx1 = 1;
                    break;
                case 11:
                    idx1 = 2;
                    break;
                default:
                    idx1 = 3;
                    break;
                }
                TYPE tmp = zero;
                TYPE tmp1 = zero; // Used to eliminate temporary storage value during leveling
                TYPE tmp2 = zero; // Used to eliminate temporary storage value during leveling
                // Step-by-step mapping
                for (int idx2 = 0; idx2 < SAT_POS_MSG + 1; idx2++) {
                    __mmask32 mask = VECTOR_EQ_MASK(vAbs, VECTOR_SET1(idx2)); // Filter out the part equal to idx2 and map it
#if EF_ELIMINATION >= 1
                    __mmask32 mask_iter = 0;
                    if (nombre_iterations <= floor_iter_thresh) {
                        mask_iter = UINT32_MAX;
                    } else {
                        mask_iter = 0;
                    }

                    __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 2:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 3:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 4:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 5:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    default:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    }
#else
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        break;
                    case 2:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        break;
                    case 3:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        break;
                    case 4:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        break;
                    case 5:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        break;
                    default:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        break;
                    }
#endif
                }
                // Handle overflow, because the magnitude of abs may overflow the bit width
                __mmask32 mask = VECTOR_GE_MASK(vAbs, VECTOR_SET1(SAT_POS_MSG + 1));
#if EF_ELIMINATION >= 1
                __mmask32 mask_iter = 0;
                if (nombre_iterations <= floor_iter_thresh) {
                    mask_iter = UINT32_MAX;
                } else {
                    mask_iter = 0;
                }

                __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 2:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 3:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 4:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 5:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                default:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                }
#else
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    break;
                case 2:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    break;
                case 3:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    break;
                case 4:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    break;
                case 5:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    break;
                default:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    break;
                }
#endif

                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(tmp, min1);
                min2 = VECTOR_MIN_2(tmp, vTemp, min2); // the second minimum value of the Lnm

                tmp = VECTOR_ADD_MASK(VECTOR_EQ_MASK(vContr, zero), VECTOR_SIGN(tmp, vContr), tmp);
                tab_vContr[j] = tmp;
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_7; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_7]), _MM_HINT_T0);
#endif
            // WE SATURATE DIRECTLY IN MSG FORMAT
#if OMS_MODE == 0 // simple OMS
            min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2 = VECTOR_MIN(min2, max_msg);
            TYPE cste_1 = VECTOR_MIN(VECTOR_SUB(min2, VECTOR_SET1(offset)), max_msg);
            TYPE cste_2 = VECTOR_MIN(VECTOR_SUB(min1, VECTOR_SET1(offset)), max_msg);
#elif OMS_MODE == 1 // selective OMS

            TYPE min1_offed, min2_offed;
            min1_offed = min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2_offed = min2 = VECTOR_MIN(min2, max_msg);

            if (nombre_iterations <= floor_iter_thresh) {
                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min1, factor_2);
                min1_offed = VECTOR_ADD_MASK(msk_lt6_1, min1, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min1, factor_1);
                min1_offed = VECTOR_ADD_MASK(msk_le1_1, min1_offed, ones); // min_offed = min_offed + 1 = min + 2

                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min2, factor_2);
                min2_offed = VECTOR_ADD_MASK(msk_lt6_2, min2, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min2, factor_1);
                min2_offed = VECTOR_ADD_MASK(msk_le1_2, min2_offed, ones); // min_offed = min_offed + 1 = min + 2
            }
            TYPE cste_1 = VECTOR_MIN(min2_offed, max_msg);
            TYPE cste_2 = VECTOR_MIN(min1_offed, max_msg);

#endif
#if (DEG_7 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_7)
            for (int j = 0; j < DEG_7; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(temp_vContr[j], v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
            cn_count++;
        }
#endif
#if NB_DEGRES >= 8
        /*************************************DEG_8*************************************************/
        for (int i = 0; i < DEG_8_COMPUTATIONS; i++) {

#if (DEG_8 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_8];
            TYPE temp_vContr[DEG_8];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_8)
            for (int j = 0; j < DEG_8; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                unsigned short col = PosNoeudsVariable[p_indice_nod1 - p_vn_adr];

                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                vContr = VECTOR_MIN(vContr, max_var);

                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
                sign = VECTOR_XOR(sign, cSign); // 0: pos, 1: neg
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                temp_vContr[j] = vContr;

                int idx1 = 0;
                switch (VN_weight_[col]) {
                case 3:
                    idx1 = 0;
                    break;
                case 6:
                    idx1 = 1;
                    break;
                case 11:
                    idx1 = 2;
                    break;
                default:
                    idx1 = 3;
                    break;
                }
                TYPE tmp = zero;
                TYPE tmp1 = zero; // Used to eliminate temporary storage value during leveling
                TYPE tmp2 = zero; // Used to eliminate temporary storage value during leveling
                // Step-by-step mapping
                for (int idx2 = 0; idx2 < SAT_POS_MSG + 1; idx2++) {
                    __mmask32 mask = VECTOR_EQ_MASK(vAbs, VECTOR_SET1(idx2)); // Filter out the part equal to idx2 and map it
#if EF_ELIMINATION >= 1
                    __mmask32 mask_iter = 0;
                    if (nombre_iterations <= floor_iter_thresh) {
                        mask_iter = UINT32_MAX;
                    } else {
                        mask_iter = 0;
                    }

                    __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 2:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 3:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 4:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 5:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    default:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    }
#else
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        break;
                    case 2:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        break;
                    case 3:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        break;
                    case 4:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        break;
                    case 5:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        break;
                    default:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        break;
                    }
#endif
                }
                // Handle overflow, because the magnitude of abs may overflow the bit width
                __mmask32 mask = VECTOR_GE_MASK(vAbs, VECTOR_SET1(SAT_POS_MSG + 1));
#if EF_ELIMINATION >= 1
                __mmask32 mask_iter = 0;
                if (nombre_iterations <= floor_iter_thresh) {
                    mask_iter = UINT32_MAX;
                } else {
                    mask_iter = 0;
                }

                __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 2:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 3:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 4:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 5:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                default:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                }
#else
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    break;
                case 2:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    break;
                case 3:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    break;
                case 4:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    break;
                case 5:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    break;
                default:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    break;
                }
#endif

                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(tmp, min1);
                min2 = VECTOR_MIN_2(tmp, vTemp, min2); // the second minimum value of the Lnm

                tmp = VECTOR_ADD_MASK(VECTOR_EQ_MASK(vContr, zero), VECTOR_SIGN(tmp, vContr), tmp);
                tab_vContr[j] = tmp;
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_8; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_8]), _MM_HINT_T0);
#endif
            // WE SATURATE DIRECTLY IN MSG FORMAT
#if OMS_MODE == 0 // simple OMS
            min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2 = VECTOR_MIN(min2, max_msg);
            TYPE cste_1 = VECTOR_MIN(VECTOR_SUB(min2, VECTOR_SET1(offset)), max_msg);
            TYPE cste_2 = VECTOR_MIN(VECTOR_SUB(min1, VECTOR_SET1(offset)), max_msg);
#elif OMS_MODE == 1 // selective OMS

            TYPE min1_offed, min2_offed;
            min1_offed = min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2_offed = min2 = VECTOR_MIN(min2, max_msg);

            if (nombre_iterations <= floor_iter_thresh) {
                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min1, factor_2);
                min1_offed = VECTOR_ADD_MASK(msk_lt6_1, min1, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min1, factor_1);
                min1_offed = VECTOR_ADD_MASK(msk_le1_1, min1_offed, ones); // min_offed = min_offed + 1 = min + 2

                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min2, factor_2);
                min2_offed = VECTOR_ADD_MASK(msk_lt6_2, min2, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min2, factor_1);
                min2_offed = VECTOR_ADD_MASK(msk_le1_2, min2_offed, ones); // min_offed = min_offed + 1 = min + 2
            }
            TYPE cste_1 = VECTOR_MIN(min2_offed, max_msg);
            TYPE cste_2 = VECTOR_MIN(min1_offed, max_msg);

#endif
#if (DEG_8 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_8)
            for (int j = 0; j < DEG_8; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(temp_vContr[j], v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
            cn_count++;
        }
#endif
#if NB_DEGRES >= 9
        /*************************************DEG_9*************************************************/
        for (int i = 0; i < DEG_9_COMPUTATIONS; i++) {

#if (DEG_9 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_9];
            TYPE temp_vContr[DEG_9];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_9)
            for (int j = 0; j < DEG_9; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                unsigned short col = PosNoeudsVariable[p_indice_nod1 - p_vn_adr];

                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                vContr = VECTOR_MIN(vContr, max_var);

                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
                sign = VECTOR_XOR(sign, cSign); // 0: pos, 1: neg
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                temp_vContr[j] = vContr;

                int idx1 = 0;
                switch (VN_weight_[col]) {
                case 3:
                    idx1 = 0;
                    break;
                case 6:
                    idx1 = 1;
                    break;
                case 11:
                    idx1 = 2;
                    break;
                default:
                    idx1 = 3;
                    break;
                }
                TYPE tmp = zero;
                TYPE tmp1 = zero; // Used to eliminate temporary storage value during leveling
                TYPE tmp2 = zero; // Used to eliminate temporary storage value during leveling
                // Step-by-step mapping
                for (int idx2 = 0; idx2 < SAT_POS_MSG + 1; idx2++) {
                    __mmask32 mask = VECTOR_EQ_MASK(vAbs, VECTOR_SET1(idx2)); // Filter out the part equal to idx2 and map it
#if EF_ELIMINATION >= 1
                    __mmask32 mask_iter = 0;
                    if (nombre_iterations <= floor_iter_thresh) {
                        mask_iter = UINT32_MAX;
                    } else {
                        mask_iter = 0;
                    }

                    __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 2:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 3:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 4:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 5:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    default:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    }
#else
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        break;
                    case 2:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        break;
                    case 3:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        break;
                    case 4:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        break;
                    case 5:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        break;
                    default:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        break;
                    }
#endif
                }
                // Handle overflow, because the magnitude of abs may overflow the bit width
                __mmask32 mask = VECTOR_GE_MASK(vAbs, VECTOR_SET1(SAT_POS_MSG + 1));
#if EF_ELIMINATION >= 1
                __mmask32 mask_iter = 0;
                if (nombre_iterations <= floor_iter_thresh) {
                    mask_iter = UINT32_MAX;
                } else {
                    mask_iter = 0;
                }

                __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 2:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 3:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 4:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 5:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                default:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                }
#else
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    break;
                case 2:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    break;
                case 3:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    break;
                case 4:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    break;
                case 5:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    break;
                default:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    break;
                }
#endif

                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(tmp, min1);
                min2 = VECTOR_MIN_2(tmp, vTemp, min2); // the second minimum value of the Lnm

                tmp = VECTOR_ADD_MASK(VECTOR_EQ_MASK(vContr, zero), VECTOR_SIGN(tmp, vContr), tmp);
                tab_vContr[j] = tmp;
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_9; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_9]), _MM_HINT_T0);
#endif
            // WE SATURATE DIRECTLY IN MSG FORMAT
#if OMS_MODE == 0 // simple OMS
            min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2 = VECTOR_MIN(min2, max_msg);
            TYPE cste_1 = VECTOR_MIN(VECTOR_SUB(min2, VECTOR_SET1(offset)), max_msg);
            TYPE cste_2 = VECTOR_MIN(VECTOR_SUB(min1, VECTOR_SET1(offset)), max_msg);
#elif OMS_MODE == 1 // selective OMS

            TYPE min1_offed, min2_offed;
            min1_offed = min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2_offed = min2 = VECTOR_MIN(min2, max_msg);

            if (nombre_iterations <= floor_iter_thresh) {
                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min1, factor_2);
                min1_offed = VECTOR_ADD_MASK(msk_lt6_1, min1, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min1, factor_1);
                min1_offed = VECTOR_ADD_MASK(msk_le1_1, min1_offed, ones); // min_offed = min_offed + 1 = min + 2

                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min2, factor_2);
                min2_offed = VECTOR_ADD_MASK(msk_lt6_2, min2, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min2, factor_1);
                min2_offed = VECTOR_ADD_MASK(msk_le1_2, min2_offed, ones); // min_offed = min_offed + 1 = min + 2
            }
            TYPE cste_1 = VECTOR_MIN(min2_offed, max_msg);
            TYPE cste_2 = VECTOR_MIN(min1_offed, max_msg);

#endif
#if (DEG_9 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_9)
            for (int j = 0; j < DEG_9; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(temp_vContr[j], v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
            cn_count++;
        }
#endif
#if NB_DEGRES >= 10
        /*************************************DEG_10*************************************************/
        for (int i = 0; i < DEG_10_COMPUTATIONS; i++) {

#if (DEG_10 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_10];
            TYPE temp_vContr[DEG_10];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_10)
            for (int j = 0; j < DEG_10; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                unsigned short col = PosNoeudsVariable[p_indice_nod1 - p_vn_adr];

                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                vContr = VECTOR_MIN(vContr, max_var);

                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
                sign = VECTOR_XOR(sign, cSign); // 0: pos, 1: neg
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                temp_vContr[j] = vContr;

                int idx1 = 0;
                switch (VN_weight_[col]) {
                case 3:
                    idx1 = 0;
                    break;
                case 6:
                    idx1 = 1;
                    break;
                case 11:
                    idx1 = 2;
                    break;
                default:
                    idx1 = 3;
                    break;
                }
                TYPE tmp = zero;
                TYPE tmp1 = zero; // Used to eliminate temporary storage value during leveling
                TYPE tmp2 = zero; // Used to eliminate temporary storage value during leveling
                // Step-by-step mapping
                for (int idx2 = 0; idx2 < SAT_POS_MSG + 1; idx2++) {
                    __mmask32 mask = VECTOR_EQ_MASK(vAbs, VECTOR_SET1(idx2)); // Filter out the part equal to idx2 and map it
#if EF_ELIMINATION >= 1
                    __mmask32 mask_iter = 0;
                    if (nombre_iterations <= floor_iter_thresh) {
                        mask_iter = UINT32_MAX;
                    } else {
                        mask_iter = 0;
                    }

                    __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 2:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 3:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 4:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 5:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    default:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    }
#else
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        break;
                    case 2:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        break;
                    case 3:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        break;
                    case 4:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        break;
                    case 5:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        break;
                    default:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        break;
                    }
#endif
                }
                // Handle overflow, because the magnitude of abs may overflow the bit width
                __mmask32 mask = VECTOR_GE_MASK(vAbs, VECTOR_SET1(SAT_POS_MSG + 1));
#if EF_ELIMINATION >= 1
                __mmask32 mask_iter = 0;
                if (nombre_iterations <= floor_iter_thresh) {
                    mask_iter = UINT32_MAX;
                } else {
                    mask_iter = 0;
                }

                __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 2:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 3:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 4:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 5:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                default:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                }
#else
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    break;
                case 2:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    break;
                case 3:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    break;
                case 4:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    break;
                case 5:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    break;
                default:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    break;
                }
#endif

                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(tmp, min1);
                min2 = VECTOR_MIN_2(tmp, vTemp, min2); // the second minimum value of the Lnm

                tmp = VECTOR_ADD_MASK(VECTOR_EQ_MASK(vContr, zero), VECTOR_SIGN(tmp, vContr), tmp);
                tab_vContr[j] = tmp;
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_10; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_10]), _MM_HINT_T0);
#endif
            // WE SATURATE DIRECTLY IN MSG FORMAT
#if OMS_MODE == 0 // simple OMS
            min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2 = VECTOR_MIN(min2, max_msg);
            TYPE cste_1 = VECTOR_MIN(VECTOR_SUB(min2, VECTOR_SET1(offset)), max_msg);
            TYPE cste_2 = VECTOR_MIN(VECTOR_SUB(min1, VECTOR_SET1(offset)), max_msg);
#elif OMS_MODE == 1 // selective OMS

            TYPE min1_offed, min2_offed;
            min1_offed = min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2_offed = min2 = VECTOR_MIN(min2, max_msg);

            if (nombre_iterations <= floor_iter_thresh) {
                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min1, factor_2);
                min1_offed = VECTOR_ADD_MASK(msk_lt6_1, min1, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min1, factor_1);
                min1_offed = VECTOR_ADD_MASK(msk_le1_1, min1_offed, ones); // min_offed = min_offed + 1 = min + 2

                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min2, factor_2);
                min2_offed = VECTOR_ADD_MASK(msk_lt6_2, min2, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min2, factor_1);
                min2_offed = VECTOR_ADD_MASK(msk_le1_2, min2_offed, ones); // min_offed = min_offed + 1 = min + 2
            }
            TYPE cste_1 = VECTOR_MIN(min2_offed, max_msg);
            TYPE cste_2 = VECTOR_MIN(min1_offed, max_msg);

#endif
#if (DEG_10 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_10)
            for (int j = 0; j < DEG_10; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(temp_vContr[j], v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
            cn_count++;
        }
#endif
#if NB_DEGRES >= 11
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_11_COMPUTATIONS; i++) {

#if (DEG_11 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_11];
            TYPE temp_vContr[DEG_11];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_11)
            for (int j = 0; j < DEG_11; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                unsigned short col = PosNoeudsVariable[p_indice_nod1 - p_vn_adr];

                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                vContr = VECTOR_MIN(vContr, max_var);

                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
                sign = VECTOR_XOR(sign, cSign); // 0: pos, 1: neg
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                temp_vContr[j] = vContr;

                int idx1 = 0;
                switch (VN_weight_[col]) {
                case 3:
                    idx1 = 0;
                    break;
                case 6:
                    idx1 = 1;
                    break;
                case 11:
                    idx1 = 2;
                    break;
                default:
                    idx1 = 3;
                    break;
                }
                TYPE tmp = zero;
                TYPE tmp1 = zero; // Used to eliminate temporary storage value during leveling
                TYPE tmp2 = zero; // Used to eliminate temporary storage value during leveling
                // Step-by-step mapping
                for (int idx2 = 0; idx2 < SAT_POS_MSG + 1; idx2++) {
                    __mmask32 mask = VECTOR_EQ_MASK(vAbs, VECTOR_SET1(idx2)); // Filter out the part equal to idx2 and map it
#if EF_ELIMINATION >= 1
                    __mmask32 mask_iter = 0;
                    if (nombre_iterations <= floor_iter_thresh) {
                        mask_iter = UINT32_MAX;
                    } else {
                        mask_iter = 0;
                    }

                    __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 2:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 3:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 4:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 5:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    default:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    }
#else
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        break;
                    case 2:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        break;
                    case 3:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        break;
                    case 4:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        break;
                    case 5:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        break;
                    default:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        break;
                    }
#endif
                }
                // Handle overflow, because the magnitude of abs may overflow the bit width
                __mmask32 mask = VECTOR_GE_MASK(vAbs, VECTOR_SET1(SAT_POS_MSG + 1));
#if EF_ELIMINATION >= 1
                __mmask32 mask_iter = 0;
                if (nombre_iterations <= floor_iter_thresh) {
                    mask_iter = UINT32_MAX;
                } else {
                    mask_iter = 0;
                }

                __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 2:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 3:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 4:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 5:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                default:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                }
#else
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    break;
                case 2:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    break;
                case 3:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    break;
                case 4:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    break;
                case 5:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    break;
                default:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    break;
                }
#endif

                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(tmp, min1);
                min2 = VECTOR_MIN_2(tmp, vTemp, min2); // the second minimum value of the Lnm

                tmp = VECTOR_ADD_MASK(VECTOR_EQ_MASK(vContr, zero), VECTOR_SIGN(tmp, vContr), tmp);
                tab_vContr[j] = tmp;
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_11; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_11]), _MM_HINT_T0);
#endif
            // WE SATURATE DIRECTLY IN MSG FORMAT
#if OMS_MODE == 0 // simple OMS
            min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2 = VECTOR_MIN(min2, max_msg);
            TYPE cste_1 = VECTOR_MIN(VECTOR_SUB(min2, VECTOR_SET1(offset)), max_msg);
            TYPE cste_2 = VECTOR_MIN(VECTOR_SUB(min1, VECTOR_SET1(offset)), max_msg);
#elif OMS_MODE == 1 // selective OMS

            TYPE min1_offed, min2_offed;
            min1_offed = min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2_offed = min2 = VECTOR_MIN(min2, max_msg);

            if (nombre_iterations <= floor_iter_thresh) {
                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min1, factor_2);
                min1_offed = VECTOR_ADD_MASK(msk_lt6_1, min1, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min1, factor_1);
                min1_offed = VECTOR_ADD_MASK(msk_le1_1, min1_offed, ones); // min_offed = min_offed + 1 = min + 2

                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min2, factor_2);
                min2_offed = VECTOR_ADD_MASK(msk_lt6_2, min2, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min2, factor_1);
                min2_offed = VECTOR_ADD_MASK(msk_le1_2, min2_offed, ones); // min_offed = min_offed + 1 = min + 2
            }
            TYPE cste_1 = VECTOR_MIN(min2_offed, max_msg);
            TYPE cste_2 = VECTOR_MIN(min1_offed, max_msg);

#endif
#if (DEG_11 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_11)
            for (int j = 0; j < DEG_11; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(temp_vContr[j], v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
            cn_count++;
        }
#endif
#if NB_DEGRES >= 12
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_12_COMPUTATIONS; i++) {

#if (DEG_12 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_12];
            TYPE temp_vContr[DEG_12];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_12)
            for (int j = 0; j < DEG_12; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                unsigned short col = PosNoeudsVariable[p_indice_nod1 - p_vn_adr];

                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                vContr = VECTOR_MIN(vContr, max_var);

                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
                sign = VECTOR_XOR(sign, cSign); // 0: pos, 1: neg
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                temp_vContr[j] = vContr;

                int idx1 = 0;
                switch (VN_weight_[col]) {
                case 3:
                    idx1 = 0;
                    break;
                case 6:
                    idx1 = 1;
                    break;
                case 11:
                    idx1 = 2;
                    break;
                default:
                    idx1 = 3;
                    break;
                }
                TYPE tmp = zero;
                TYPE tmp1 = zero; // Used to eliminate temporary storage value during leveling
                TYPE tmp2 = zero; // Used to eliminate temporary storage value during leveling
                // Step-by-step mapping
                for (int idx2 = 0; idx2 < SAT_POS_MSG + 1; idx2++) {
                    __mmask32 mask = VECTOR_EQ_MASK(vAbs, VECTOR_SET1(idx2)); // Filter out the part equal to idx2 and map it
#if EF_ELIMINATION >= 1
                    __mmask32 mask_iter = 0;
                    if (nombre_iterations <= floor_iter_thresh) {
                        mask_iter = UINT32_MAX;
                    } else {
                        mask_iter = 0;
                    }

                    __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 2:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 3:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 4:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 5:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    default:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    }
#else
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        break;
                    case 2:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        break;
                    case 3:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        break;
                    case 4:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        break;
                    case 5:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        break;
                    default:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        break;
                    }
#endif
                }
                // Handle overflow, because the magnitude of abs may overflow the bit width
                __mmask32 mask = VECTOR_GE_MASK(vAbs, VECTOR_SET1(SAT_POS_MSG + 1));
#if EF_ELIMINATION >= 1
                __mmask32 mask_iter = 0;
                if (nombre_iterations <= floor_iter_thresh) {
                    mask_iter = UINT32_MAX;
                } else {
                    mask_iter = 0;
                }

                __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 2:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 3:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 4:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 5:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                default:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                }
#else
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    break;
                case 2:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    break;
                case 3:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    break;
                case 4:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    break;
                case 5:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    break;
                default:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    break;
                }
#endif

                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(tmp, min1);
                min2 = VECTOR_MIN_2(tmp, vTemp, min2); // the second minimum value of the Lnm

                tmp = VECTOR_ADD_MASK(VECTOR_EQ_MASK(vContr, zero), VECTOR_SIGN(tmp, vContr), tmp);
                tab_vContr[j] = tmp;
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_12; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_12]), _MM_HINT_T0);
#endif
            // WE SATURATE DIRECTLY IN MSG FORMAT
#if OMS_MODE == 0 // simple OMS
            min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2 = VECTOR_MIN(min2, max_msg);
            TYPE cste_1 = VECTOR_MIN(VECTOR_SUB(min2, VECTOR_SET1(offset)), max_msg);
            TYPE cste_2 = VECTOR_MIN(VECTOR_SUB(min1, VECTOR_SET1(offset)), max_msg);
#elif OMS_MODE == 1 // selective OMS

            TYPE min1_offed, min2_offed;
            min1_offed = min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2_offed = min2 = VECTOR_MIN(min2, max_msg);

            if (nombre_iterations <= floor_iter_thresh) {
                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min1, factor_2);
                min1_offed = VECTOR_ADD_MASK(msk_lt6_1, min1, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min1, factor_1);
                min1_offed = VECTOR_ADD_MASK(msk_le1_1, min1_offed, ones); // min_offed = min_offed + 1 = min + 2

                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min2, factor_2);
                min2_offed = VECTOR_ADD_MASK(msk_lt6_2, min2, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min2, factor_1);
                min2_offed = VECTOR_ADD_MASK(msk_le1_2, min2_offed, ones); // min_offed = min_offed + 1 = min + 2
            }
            TYPE cste_1 = VECTOR_MIN(min2_offed, max_msg);
            TYPE cste_2 = VECTOR_MIN(min1_offed, max_msg);

#endif
#if (DEG_12 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_12)
            for (int j = 0; j < DEG_12; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(temp_vContr[j], v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
            cn_count++;
        }
#endif

#if NB_DEGRES >= 13
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_13_COMPUTATIONS; i++) {

#if (DEG_13 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_13];
            TYPE temp_vContr[DEG_13];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_13)
            for (int j = 0; j < DEG_13; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                unsigned short col = PosNoeudsVariable[p_indice_nod1 - p_vn_adr];

                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                vContr = VECTOR_MIN(vContr, max_var);

                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
                sign = VECTOR_XOR(sign, cSign); // 0: pos, 1: neg
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                temp_vContr[j] = vContr;

                int idx1 = 0;
                switch (VN_weight_[col]) {
                case 3:
                    idx1 = 0;
                    break;
                case 6:
                    idx1 = 1;
                    break;
                case 11:
                    idx1 = 2;
                    break;
                default:
                    idx1 = 3;
                    break;
                }
                TYPE tmp = zero;
                TYPE tmp1 = zero; // Used to eliminate temporary storage value during leveling
                TYPE tmp2 = zero; // Used to eliminate temporary storage value during leveling
                // Step-by-step mapping
                for (int idx2 = 0; idx2 < SAT_POS_MSG + 1; idx2++) {
                    __mmask32 mask = VECTOR_EQ_MASK(vAbs, VECTOR_SET1(idx2)); // Filter out the part equal to idx2 and map it
#if EF_ELIMINATION >= 1
                    __mmask32 mask_iter = 0;
                    if (nombre_iterations <= floor_iter_thresh) {
                        mask_iter = UINT32_MAX;
                    } else {
                        mask_iter = 0;
                    }

                    __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 2:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 3:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 4:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 5:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    default:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    }
#else
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        break;
                    case 2:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        break;
                    case 3:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        break;
                    case 4:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        break;
                    case 5:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        break;
                    default:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        break;
                    }
#endif
                }
                // Handle overflow, because the magnitude of abs may overflow the bit width
                __mmask32 mask = VECTOR_GE_MASK(vAbs, VECTOR_SET1(SAT_POS_MSG + 1));
#if EF_ELIMINATION >= 1
                __mmask32 mask_iter = 0;
                if (nombre_iterations <= floor_iter_thresh) {
                    mask_iter = UINT32_MAX;
                } else {
                    mask_iter = 0;
                }

                __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 2:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 3:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 4:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 5:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                default:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                }
#else
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    break;
                case 2:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    break;
                case 3:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    break;
                case 4:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    break;
                case 5:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    break;
                default:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    break;
                }
#endif

                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(tmp, min1);
                min2 = VECTOR_MIN_2(tmp, vTemp, min2); // the second minimum value of the Lnm

                tmp = VECTOR_ADD_MASK(VECTOR_EQ_MASK(vContr, zero), VECTOR_SIGN(tmp, vContr), tmp);
                tab_vContr[j] = tmp;
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_13; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_13]), _MM_HINT_T0);
#endif
            // WE SATURATE DIRECTLY IN MSG FORMAT
#if OMS_MODE == 0 // simple OMS
            min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2 = VECTOR_MIN(min2, max_msg);
            TYPE cste_1 = VECTOR_MIN(VECTOR_SUB(min2, VECTOR_SET1(offset)), max_msg);
            TYPE cste_2 = VECTOR_MIN(VECTOR_SUB(min1, VECTOR_SET1(offset)), max_msg);
#elif OMS_MODE == 1 // selective OMS

            TYPE min1_offed, min2_offed;
            min1_offed = min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2_offed = min2 = VECTOR_MIN(min2, max_msg);

            if (nombre_iterations <= floor_iter_thresh) {
                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min1, factor_2);
                min1_offed = VECTOR_ADD_MASK(msk_lt6_1, min1, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min1, factor_1);
                min1_offed = VECTOR_ADD_MASK(msk_le1_1, min1_offed, ones); // min_offed = min_offed + 1 = min + 2

                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min2, factor_2);
                min2_offed = VECTOR_ADD_MASK(msk_lt6_2, min2, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min2, factor_1);
                min2_offed = VECTOR_ADD_MASK(msk_le1_2, min2_offed, ones); // min_offed = min_offed + 1 = min + 2
            }
            TYPE cste_1 = VECTOR_MIN(min2_offed, max_msg);
            TYPE cste_2 = VECTOR_MIN(min1_offed, max_msg);

#endif
#if (DEG_13 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_13)
            for (int j = 0; j < DEG_13; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(temp_vContr[j], v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
            cn_count++;
        }
#endif
#if NB_DEGRES >= 14
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_14_COMPUTATIONS; i++) {

#if (DEG_14 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_14];
            TYPE temp_vContr[DEG_14];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_14)
            for (int j = 0; j < DEG_14; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                unsigned short col = PosNoeudsVariable[p_indice_nod1 - p_vn_adr];

                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                vContr = VECTOR_MIN(vContr, max_var);

                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
                sign = VECTOR_XOR(sign, cSign); // 0: pos, 1: neg
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                temp_vContr[j] = vContr;

                int idx1 = 0;
                switch (VN_weight_[col]) {
                case 3:
                    idx1 = 0;
                    break;
                case 6:
                    idx1 = 1;
                    break;
                case 11:
                    idx1 = 2;
                    break;
                default:
                    idx1 = 3;
                    break;
                }
                TYPE tmp = zero;
                TYPE tmp1 = zero; // Used to eliminate temporary storage value during leveling
                TYPE tmp2 = zero; // Used to eliminate temporary storage value during leveling
                // Step-by-step mapping
                for (int idx2 = 0; idx2 < SAT_POS_MSG + 1; idx2++) {
                    __mmask32 mask = VECTOR_EQ_MASK(vAbs, VECTOR_SET1(idx2)); // Filter out the part equal to idx2 and map it
#if EF_ELIMINATION >= 1
                    __mmask32 mask_iter = 0;
                    if (nombre_iterations <= floor_iter_thresh) {
                        mask_iter = UINT32_MAX;
                    } else {
                        mask_iter = 0;
                    }

                    __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 2:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 3:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 4:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 5:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    default:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    }
#else
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        break;
                    case 2:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        break;
                    case 3:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        break;
                    case 4:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        break;
                    case 5:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        break;
                    default:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        break;
                    }
#endif
                }
                // Handle overflow, because the magnitude of abs may overflow the bit width
                __mmask32 mask = VECTOR_GE_MASK(vAbs, VECTOR_SET1(SAT_POS_MSG + 1));
#if EF_ELIMINATION >= 1
                __mmask32 mask_iter = 0;
                if (nombre_iterations <= floor_iter_thresh) {
                    mask_iter = UINT32_MAX;
                } else {
                    mask_iter = 0;
                }

                __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 2:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 3:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 4:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 5:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                default:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                }
#else
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    break;
                case 2:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    break;
                case 3:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    break;
                case 4:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    break;
                case 5:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    break;
                default:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    break;
                }
#endif

                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(tmp, min1);
                min2 = VECTOR_MIN_2(tmp, vTemp, min2); // the second minimum value of the Lnm

                tmp = VECTOR_ADD_MASK(VECTOR_EQ_MASK(vContr, zero), VECTOR_SIGN(tmp, vContr), tmp);
                tab_vContr[j] = tmp;
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_14; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_14]), _MM_HINT_T0);
#endif
            // WE SATURATE DIRECTLY IN MSG FORMAT
#if OMS_MODE == 0 // simple OMS
            min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2 = VECTOR_MIN(min2, max_msg);
            TYPE cste_1 = VECTOR_MIN(VECTOR_SUB(min2, VECTOR_SET1(offset)), max_msg);
            TYPE cste_2 = VECTOR_MIN(VECTOR_SUB(min1, VECTOR_SET1(offset)), max_msg);
#elif OMS_MODE == 1 // selective OMS

            TYPE min1_offed, min2_offed;
            min1_offed = min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2_offed = min2 = VECTOR_MIN(min2, max_msg);

            if (nombre_iterations <= floor_iter_thresh) {
                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min1, factor_2);
                min1_offed = VECTOR_ADD_MASK(msk_lt6_1, min1, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min1, factor_1);
                min1_offed = VECTOR_ADD_MASK(msk_le1_1, min1_offed, ones); // min_offed = min_offed + 1 = min + 2

                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min2, factor_2);
                min2_offed = VECTOR_ADD_MASK(msk_lt6_2, min2, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min2, factor_1);
                min2_offed = VECTOR_ADD_MASK(msk_le1_2, min2_offed, ones); // min_offed = min_offed + 1 = min + 2
            }
            TYPE cste_1 = VECTOR_MIN(min2_offed, max_msg);
            TYPE cste_2 = VECTOR_MIN(min1_offed, max_msg);

#endif
#if (DEG_14 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_14)
            for (int j = 0; j < DEG_14; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(temp_vContr[j], v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
            cn_count++;
        }
#endif
#if NB_DEGRES >= 15
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_15_COMPUTATIONS; i++) {

#if (DEG_15 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_15];
            TYPE temp_vContr[DEG_15];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_15)
            for (int j = 0; j < DEG_15; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                unsigned short col = PosNoeudsVariable[p_indice_nod1 - p_vn_adr];

                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                vContr = VECTOR_MIN(vContr, max_var);

                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
                sign = VECTOR_XOR(sign, cSign); // 0: pos, 1: neg
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                temp_vContr[j] = vContr;

                int idx1 = 0;
                switch (VN_weight_[col]) {
                case 3:
                    idx1 = 0;
                    break;
                case 6:
                    idx1 = 1;
                    break;
                case 11:
                    idx1 = 2;
                    break;
                default:
                    idx1 = 3;
                    break;
                }
                TYPE tmp = zero;
                TYPE tmp1 = zero; // Used to eliminate temporary storage value during leveling
                TYPE tmp2 = zero; // Used to eliminate temporary storage value during leveling
                // Step-by-step mapping
                for (int idx2 = 0; idx2 < SAT_POS_MSG + 1; idx2++) {
                    __mmask32 mask = VECTOR_EQ_MASK(vAbs, VECTOR_SET1(idx2)); // Filter out the part equal to idx2 and map it
#if EF_ELIMINATION >= 1
                    __mmask32 mask_iter = 0;
                    if (nombre_iterations <= floor_iter_thresh) {
                        mask_iter = UINT32_MAX;
                    } else {
                        mask_iter = 0;
                    }

                    __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 2:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 3:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 4:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 5:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    default:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    }
#else
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        break;
                    case 2:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        break;
                    case 3:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        break;
                    case 4:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        break;
                    case 5:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        break;
                    default:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        break;
                    }
#endif
                }
                // Handle overflow, because the magnitude of abs may overflow the bit width
                __mmask32 mask = VECTOR_GE_MASK(vAbs, VECTOR_SET1(SAT_POS_MSG + 1));
#if EF_ELIMINATION >= 1
                __mmask32 mask_iter = 0;
                if (nombre_iterations <= floor_iter_thresh) {
                    mask_iter = UINT32_MAX;
                } else {
                    mask_iter = 0;
                }

                __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 2:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 3:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 4:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 5:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                default:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                }
#else
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    break;
                case 2:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    break;
                case 3:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    break;
                case 4:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    break;
                case 5:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    break;
                default:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    break;
                }
#endif

                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(tmp, min1);
                min2 = VECTOR_MIN_2(tmp, vTemp, min2); // the second minimum value of the Lnm

                tmp = VECTOR_ADD_MASK(VECTOR_EQ_MASK(vContr, zero), VECTOR_SIGN(tmp, vContr), tmp);
                tab_vContr[j] = tmp;
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_15; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_15]), _MM_HINT_T0);
#endif
            // WE SATURATE DIRECTLY IN MSG FORMAT
#if OMS_MODE == 0 // simple OMS
            min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2 = VECTOR_MIN(min2, max_msg);
            TYPE cste_1 = VECTOR_MIN(VECTOR_SUB(min2, VECTOR_SET1(offset)), max_msg);
            TYPE cste_2 = VECTOR_MIN(VECTOR_SUB(min1, VECTOR_SET1(offset)), max_msg);
#elif OMS_MODE == 1 // selective OMS

            TYPE min1_offed, min2_offed;
            min1_offed = min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2_offed = min2 = VECTOR_MIN(min2, max_msg);

            if (nombre_iterations <= floor_iter_thresh) {
                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min1, factor_2);
                min1_offed = VECTOR_ADD_MASK(msk_lt6_1, min1, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min1, factor_1);
                min1_offed = VECTOR_ADD_MASK(msk_le1_1, min1_offed, ones); // min_offed = min_offed + 1 = min + 2

                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min2, factor_2);
                min2_offed = VECTOR_ADD_MASK(msk_lt6_2, min2, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min2, factor_1);
                min2_offed = VECTOR_ADD_MASK(msk_le1_2, min2_offed, ones); // min_offed = min_offed + 1 = min + 2
            }
            TYPE cste_1 = VECTOR_MIN(min2_offed, max_msg);
            TYPE cste_2 = VECTOR_MIN(min1_offed, max_msg);

#endif
#if (DEG_15 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_15)
            for (int j = 0; j < DEG_15; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(temp_vContr[j], v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
            cn_count++;
        }
#endif
#if NB_DEGRES >= 16
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_16_COMPUTATIONS; i++) {

#if (DEG_16 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_16];
            TYPE temp_vContr[DEG_16];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_16)
            for (int j = 0; j < DEG_16; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                unsigned short col = PosNoeudsVariable[p_indice_nod1 - p_vn_adr];

                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                vContr = VECTOR_MIN(vContr, max_var);

                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
                sign = VECTOR_XOR(sign, cSign); // 0: pos, 1: neg
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                temp_vContr[j] = vContr;

                int idx1 = 0;
                switch (VN_weight_[col]) {
                case 3:
                    idx1 = 0;
                    break;
                case 6:
                    idx1 = 1;
                    break;
                case 11:
                    idx1 = 2;
                    break;
                default:
                    idx1 = 3;
                    break;
                }
                TYPE tmp = zero;
                TYPE tmp1 = zero; // Used to eliminate temporary storage value during leveling
                TYPE tmp2 = zero; // Used to eliminate temporary storage value during leveling
                // Step-by-step mapping
                for (int idx2 = 0; idx2 < SAT_POS_MSG + 1; idx2++) {
                    __mmask32 mask = VECTOR_EQ_MASK(vAbs, VECTOR_SET1(idx2)); // Filter out the part equal to idx2 and map it
#if EF_ELIMINATION >= 1
                    __mmask32 mask_iter = 0;
                    if (nombre_iterations <= floor_iter_thresh) {
                        mask_iter = UINT32_MAX;
                    } else {
                        mask_iter = 0;
                    }

                    __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 2:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 3:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 4:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 5:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    default:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    }
#else
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        break;
                    case 2:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        break;
                    case 3:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        break;
                    case 4:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        break;
                    case 5:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        break;
                    default:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        break;
                    }
#endif
                }
                // Handle overflow, because the magnitude of abs may overflow the bit width
                __mmask32 mask = VECTOR_GE_MASK(vAbs, VECTOR_SET1(SAT_POS_MSG + 1));
#if EF_ELIMINATION >= 1
                __mmask32 mask_iter = 0;
                if (nombre_iterations <= floor_iter_thresh) {
                    mask_iter = UINT32_MAX;
                } else {
                    mask_iter = 0;
                }

                __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 2:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 3:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 4:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 5:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                default:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                }
#else
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    break;
                case 2:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    break;
                case 3:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    break;
                case 4:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    break;
                case 5:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    break;
                default:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    break;
                }
#endif

                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(tmp, min1);
                min2 = VECTOR_MIN_2(tmp, vTemp, min2); // the second minimum value of the Lnm

                tmp = VECTOR_ADD_MASK(VECTOR_EQ_MASK(vContr, zero), VECTOR_SIGN(tmp, vContr), tmp);
                tab_vContr[j] = tmp;
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_16; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_16]), _MM_HINT_T0);
#endif
            // WE SATURATE DIRECTLY IN MSG FORMAT
#if OMS_MODE == 0 // simple OMS
            min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2 = VECTOR_MIN(min2, max_msg);
            TYPE cste_1 = VECTOR_MIN(VECTOR_SUB(min2, VECTOR_SET1(offset)), max_msg);
            TYPE cste_2 = VECTOR_MIN(VECTOR_SUB(min1, VECTOR_SET1(offset)), max_msg);
#elif OMS_MODE == 1 // selective OMS

            TYPE min1_offed, min2_offed;
            min1_offed = min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2_offed = min2 = VECTOR_MIN(min2, max_msg);

            if (nombre_iterations <= floor_iter_thresh) {
                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min1, factor_2);
                min1_offed = VECTOR_ADD_MASK(msk_lt6_1, min1, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min1, factor_1);
                min1_offed = VECTOR_ADD_MASK(msk_le1_1, min1_offed, ones); // min_offed = min_offed + 1 = min + 2

                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min2, factor_2);
                min2_offed = VECTOR_ADD_MASK(msk_lt6_2, min2, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min2, factor_1);
                min2_offed = VECTOR_ADD_MASK(msk_le1_2, min2_offed, ones); // min_offed = min_offed + 1 = min + 2
            }
            TYPE cste_1 = VECTOR_MIN(min2_offed, max_msg);
            TYPE cste_2 = VECTOR_MIN(min1_offed, max_msg);

#endif
#if (DEG_16 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_16)
            for (int j = 0; j < DEG_16; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(temp_vContr[j], v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
            cn_count++;
        }
#endif
#if NB_DEGRES >= 17
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_17_COMPUTATIONS; i++) {

#if (DEG_17 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_17];
            TYPE temp_vContr[DEG_17];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_17)
            for (int j = 0; j < DEG_17; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                unsigned short col = PosNoeudsVariable[p_indice_nod1 - p_vn_adr];

                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                vContr = VECTOR_MIN(vContr, max_var);

                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
                sign = VECTOR_XOR(sign, cSign); // 0: pos, 1: neg
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                temp_vContr[j] = vContr;

                int idx1 = 0;
                switch (VN_weight_[col]) {
                case 3:
                    idx1 = 0;
                    break;
                case 6:
                    idx1 = 1;
                    break;
                case 11:
                    idx1 = 2;
                    break;
                default:
                    idx1 = 3;
                    break;
                }
                TYPE tmp = zero;
                / TYPE tmp1 = zero; // Used to eliminate temporary storage value during leveling
                TYPE tmp2 = zero; // Used to eliminate temporary storage value during leveling
                // Step-by-step mapping
                for (int idx2 = 0; idx2 < SAT_POS_MSG + 1; idx2++) {
                    __mmask32 mask = VECTOR_EQ_MASK(vAbs, VECTOR_SET1(idx2)); // Filter out the part equal to idx2 and map it
#if EF_ELIMINATION >= 1
                    __mmask32 mask_iter = 0;
                    if (nombre_iterations <= floor_iter_thresh) {
                        mask_iter = UINT32_MAX;
                    } else {
                        mask_iter = 0;
                    }

                    __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 2:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 3:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 4:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 5:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    default:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    }
#else
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        break;
                    case 2:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        break;
                    case 3:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        break;
                    case 4:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        break;
                    case 5:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        break;
                    default:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        break;
                    }
#endif
                }
                // Handle overflow, because the magnitude of abs may overflow the bit width
                __mmask32 mask = VECTOR_GE_MASK(vAbs, VECTOR_SET1(SAT_POS_MSG + 1));
#if EF_ELIMINATION >= 1
                __mmask32 mask_iter = 0;
                if (nombre_iterations <= floor_iter_thresh) {
                    mask_iter = UINT32_MAX;
                } else {
                    mask_iter = 0;
                }

                __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 2:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 3:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 4:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 5:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                default:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                }
#else
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    break;
                case 2:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    break;
                case 3:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    break;
                case 4:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    break;
                case 5:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    break;
                default:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    break;
                }
#endif

                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(tmp, min1);
                min2 = VECTOR_MIN_2(tmp, vTemp, min2); // the second minimum value of the Lnm

                tmp = VECTOR_ADD_MASK(VECTOR_EQ_MASK(vContr, zero), VECTOR_SIGN(tmp, vContr), tmp);
                tab_vContr[j] = tmp;
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_17; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_17]), _MM_HINT_T0);
#endif
            // WE SATURATE DIRECTLY IN MSG FORMAT
#if OMS_MODE == 0 // simple OMS
            min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2 = VECTOR_MIN(min2, max_msg);
            TYPE cste_1 = VECTOR_MIN(VECTOR_SUB(min2, VECTOR_SET1(offset)), max_msg);
            TYPE cste_2 = VECTOR_MIN(VECTOR_SUB(min1, VECTOR_SET1(offset)), max_msg);
#elif OMS_MODE == 1 // selective OMS

            TYPE min1_offed, min2_offed;
            min1_offed = min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2_offed = min2 = VECTOR_MIN(min2, max_msg);

            if (nombre_iterations <= floor_iter_thresh) {
                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min1, factor_2);
                min1_offed = VECTOR_ADD_MASK(msk_lt6_1, min1, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min1, factor_1);
                min1_offed = VECTOR_ADD_MASK(msk_le1_1, min1_offed, ones); // min_offed = min_offed + 1 = min + 2

                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min2, factor_2);
                min2_offed = VECTOR_ADD_MASK(msk_lt6_2, min2, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min2, factor_1);
                min2_offed = VECTOR_ADD_MASK(msk_le1_2, min2_offed, ones); // min_offed = min_offed + 1 = min + 2
            }
            TYPE cste_1 = VECTOR_MIN(min2_offed, max_msg);
            TYPE cste_2 = VECTOR_MIN(min1_offed, max_msg);

#endif
#if (DEG_17 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_17)
            for (int j = 0; j < DEG_17; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(temp_vContr[j], v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
            cn_count++;
        }
#endif
#if NB_DEGRES >= 18
        /*************************************DEG_18*************************************************/
        for (int i = 0; i < DEG_18_COMPUTATIONS; i++) {

#if (DEG_18 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_18];
            TYPE temp_vContr[DEG_18];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_18)
            for (int j = 0; j < DEG_18; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                unsigned short col = PosNoeudsVariable[p_indice_nod1 - p_vn_adr];

                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                vContr = VECTOR_MIN(vContr, max_var);

                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
                sign = VECTOR_XOR(sign, cSign); // 0: pos, 1: neg
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                temp_vContr[j] = vContr;

                int idx1 = 0;
                switch (VN_weight_[col]) {
                case 3:
                    idx1 = 0;
                    break;
                case 6:
                    idx1 = 1;
                    break;
                case 11:
                    idx1 = 2;
                    break;
                default:
                    idx1 = 3;
                    break;
                }
                TYPE tmp = zero;
                TYPE tmp1 = zero; // Used to eliminate temporary storage value during leveling
                TYPE tmp2 = zero; // Used to eliminate temporary storage value during leveling
                // Step-by-step mapping
                for (int idx2 = 0; idx2 < SAT_POS_MSG + 1; idx2++) {
                    __mmask32 mask = VECTOR_EQ_MASK(vAbs, VECTOR_SET1(idx2)); // Filter out the part equal to idx2 and map it
#if EF_ELIMINATION >= 1
                    __mmask32 mask_iter = 0;
                    if (nombre_iterations <= floor_iter_thresh) {
                        mask_iter = UINT32_MAX;
                    } else {
                        mask_iter = 0;
                    }

                    __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 2:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 3:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 4:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 5:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    default:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    }
#else
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        break;
                    case 2:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        break;
                    case 3:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        break;
                    case 4:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        break;
                    case 5:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        break;
                    default:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        break;
                    }
#endif
                }
                // Handle overflow, because the magnitude of abs may overflow the bit width
                __mmask32 mask = VECTOR_GE_MASK(vAbs, VECTOR_SET1(SAT_POS_MSG + 1));
#if EF_ELIMINATION >= 1
                __mmask32 mask_iter = 0;
                if (nombre_iterations <= floor_iter_thresh) {
                    mask_iter = UINT32_MAX;
                } else {
                    mask_iter = 0;
                }

                __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 2:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 3:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 4:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 5:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                default:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                }
#else
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    break;
                case 2:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    break;
                case 3:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    break;
                case 4:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    break;
                case 5:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    break;
                default:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    break;
                }
#endif

                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(tmp, min1);
                min2 = VECTOR_MIN_2(tmp, vTemp, min2); // the second minimum value of the Lnm

                tmp = VECTOR_ADD_MASK(VECTOR_EQ_MASK(vContr, zero), VECTOR_SIGN(tmp, vContr), tmp);
                tab_vContr[j] = tmp;
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_18; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_18]), _MM_HINT_T0);
#endif
            // WE SATURATE DIRECTLY IN MSG FORMAT
#if OMS_MODE == 0 // simple OMS
            min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2 = VECTOR_MIN(min2, max_msg);
            TYPE cste_1 = VECTOR_MIN(VECTOR_SUB(min2, VECTOR_SET1(offset)), max_msg);
            TYPE cste_2 = VECTOR_MIN(VECTOR_SUB(min1, VECTOR_SET1(offset)), max_msg);
#elif OMS_MODE == 1 // selective OMS

            TYPE min1_offed, min2_offed;
            min1_offed = min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2_offed = min2 = VECTOR_MIN(min2, max_msg);

            if (nombre_iterations <= floor_iter_thresh) {
                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min1, factor_2);
                min1_offed = VECTOR_ADD_MASK(msk_lt6_1, min1, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min1, factor_1);
                min1_offed = VECTOR_ADD_MASK(msk_le1_1, min1_offed, ones); // min_offed = min_offed + 1 = min + 2

                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min2, factor_2);
                min2_offed = VECTOR_ADD_MASK(msk_lt6_2, min2, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min2, factor_1);
                min2_offed = VECTOR_ADD_MASK(msk_le1_2, min2_offed, ones); // min_offed = min_offed + 1 = min + 2
            }
            TYPE cste_1 = VECTOR_MIN(min2_offed, max_msg);
            TYPE cste_2 = VECTOR_MIN(min1_offed, max_msg);

#endif
#if (DEG_18 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_18)
            for (int j = 0; j < DEG_18; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(temp_vContr[j], v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
            cn_count++;
        }
#endif
#if NB_DEGRES >= 19
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_19_COMPUTATIONS; i++) {

#if (DEG_19 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_19];
            TYPE temp_vContr[DEG_19];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_19)
            for (int j = 0; j < DEG_19; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                unsigned short col = PosNoeudsVariable[p_indice_nod1 - p_vn_adr];

                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                vContr = VECTOR_MIN(vContr, max_var);

                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
                sign = VECTOR_XOR(sign, cSign); // 0: pos, 1: neg
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                temp_vContr[j] = vContr;

                int idx1 = 0;
                switch (VN_weight_[col]) {
                case 3:
                    idx1 = 0;
                    break;
                case 6:
                    idx1 = 1;
                    break;
                case 11:
                    idx1 = 2;
                    break;
                default:
                    idx1 = 3;
                    break;
                }
                TYPE tmp = zero;
                TYPE tmp1 = zero; // Used to eliminate temporary storage value during leveling
                TYPE tmp2 = zero; // Used to eliminate temporary storage value during leveling
                // Step-by-step mapping
                for (int idx2 = 0; idx2 < SAT_POS_MSG + 1; idx2++) {
                    __mmask32 mask = VECTOR_EQ_MASK(vAbs, VECTOR_SET1(idx2)); // Filter out the part equal to idx2 and map it
#if EF_ELIMINATION >= 1
                    __mmask32 mask_iter = 0;
                    if (nombre_iterations <= floor_iter_thresh) {
                        mask_iter = UINT32_MAX;
                    } else {
                        mask_iter = 0;
                    }

                    __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 2:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 3:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 4:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 5:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    default:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    }
#else
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        break;
                    case 2:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        break;
                    case 3:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        break;
                    case 4:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        break;
                    case 5:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        break;
                    default:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        break;
                    }
#endif
                }
                // Handle overflow, because the magnitude of abs may overflow the bit width
                __mmask32 mask = VECTOR_GE_MASK(vAbs, VECTOR_SET1(SAT_POS_MSG + 1));
#if EF_ELIMINATION >= 1
                __mmask32 mask_iter = 0;
                if (nombre_iterations <= floor_iter_thresh) {
                    mask_iter = UINT32_MAX;
                } else {
                    mask_iter = 0;
                }

                __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 2:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 3:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 4:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 5:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                default:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                }
#else
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    break;
                case 2:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    break;
                case 3:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    break;
                case 4:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    break;
                case 5:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    break;
                default:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    break;
                }
#endif

                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(tmp, min1);
                min2 = VECTOR_MIN_2(tmp, vTemp, min2); // the second minimum value of the Lnm

                tmp = VECTOR_ADD_MASK(VECTOR_EQ_MASK(vContr, zero), VECTOR_SIGN(tmp, vContr), tmp);
                tab_vContr[j] = tmp;
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_19; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_19]), _MM_HINT_T0);
#endif
            // WE SATURATE DIRECTLY IN MSG FORMAT
#if OMS_MODE == 0 // simple OMS
            min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2 = VECTOR_MIN(min2, max_msg);
            TYPE cste_1 = VECTOR_MIN(VECTOR_SUB(min2, VECTOR_SET1(offset)), max_msg);
            TYPE cste_2 = VECTOR_MIN(VECTOR_SUB(min1, VECTOR_SET1(offset)), max_msg);
#elif OMS_MODE == 1 // selective OMS

            TYPE min1_offed, min2_offed;
            min1_offed = min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2_offed = min2 = VECTOR_MIN(min2, max_msg);

            if (nombre_iterations <= floor_iter_thresh) {
                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min1, factor_2);
                min1_offed = VECTOR_ADD_MASK(msk_lt6_1, min1, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min1, factor_1);
                min1_offed = VECTOR_ADD_MASK(msk_le1_1, min1_offed, ones); // min_offed = min_offed + 1 = min + 2

                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min2, factor_2);
                min2_offed = VECTOR_ADD_MASK(msk_lt6_2, min2, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min2, factor_1);
                min2_offed = VECTOR_ADD_MASK(msk_le1_2, min2_offed, ones); // min_offed = min_offed + 1 = min + 2
            }
            TYPE cste_1 = VECTOR_MIN(min2_offed, max_msg);
            TYPE cste_2 = VECTOR_MIN(min1_offed, max_msg);

#endif
#if (DEG_19 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_19)
            for (int j = 0; j < DEG_19; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(temp_vContr[j], v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
            cn_count++;
        }
#endif
#if NB_DEGRES >= 20
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_20_COMPUTATIONS; i++) {

#if (DEG_20 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_20];
            TYPE temp_vContr[DEG_20];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_20)
            for (int j = 0; j < DEG_20; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                unsigned short col = PosNoeudsVariable[p_indice_nod1 - p_vn_adr];

                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                vContr = VECTOR_MIN(vContr, max_var);

                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
                sign = VECTOR_XOR(sign, cSign); // 0: pos, 1: neg
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                temp_vContr[j] = vContr;

                int idx1 = 0;
                switch (VN_weight_[col]) {
                case 3:
                    idx1 = 0;
                    break;
                case 6:
                    idx1 = 1;
                    break;
                case 11:
                    idx1 = 2;
                    break;
                default:
                    idx1 = 3;
                    break;
                }
                TYPE tmp = zero;
                TYPE tmp1 = zero; // Used to eliminate temporary storage value during leveling
                TYPE tmp2 = zero; // Used to eliminate temporary storage value during leveling
                // Step-by-step mapping
                for (int idx2 = 0; idx2 < SAT_POS_MSG + 1; idx2++) {
                    __mmask32 mask = VECTOR_EQ_MASK(vAbs, VECTOR_SET1(idx2)); // Filter out the part equal to idx2 and map it
#if EF_ELIMINATION >= 1
                    __mmask32 mask_iter = 0;
                    if (nombre_iterations <= floor_iter_thresh) {
                        mask_iter = UINT32_MAX;
                    } else {
                        mask_iter = 0;
                    }

                    __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 2:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 3:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 4:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    case 5:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    default:
                        tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][idx2]));
                        tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                        tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                        break;
                    }
#else
                    switch (nb_iteration - nombre_iterations) {
                    case 1:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][idx2]));
                        break;
                    case 2:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][idx2]));
                        break;
                    case 3:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][idx2]));
                        break;
                    case 4:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][idx2]));
                        break;
                    case 5:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][idx2]));
                        break;
                    default:
                        tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][idx2]));
                        break;
                    }
#endif
                }
                // Handle overflow, because the magnitude of abs may overflow the bit width
                __mmask32 mask = VECTOR_GE_MASK(vAbs, VECTOR_SET1(SAT_POS_MSG + 1));
#if EF_ELIMINATION >= 1
                __mmask32 mask_iter = 0;
                if (nombre_iterations <= floor_iter_thresh) {
                    mask_iter = UINT32_MAX;
                } else {
                    mask_iter = 0;
                }

                __mmask32 mask_eef = mask_iter & l_m_error_sum & l_checksum_[cn_count];
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it1_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 2:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it2_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 3:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it3_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 4:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it4_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                case 5:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it5_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                default:
                    tmp1 = VECTOR_ADD_MASK(mask_eef, zero, VECTOR_SET1(V2C_map_it6_ef[idx1][SAT_POS_MSG]));
                    tmp2 = VECTOR_ADD_MASK(~mask_eef, zero, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp1);
                    tmp = VECTOR_ADD_MASK(mask, tmp, tmp2);
                    break;
                }
#else
                switch (nb_iteration - nombre_iterations) {
                case 1:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it1_[idx1][SAT_POS_MSG]));
                    break;
                case 2:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it2_[idx1][SAT_POS_MSG]));
                    break;
                case 3:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it3_[idx1][SAT_POS_MSG]));
                    break;
                case 4:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it4_[idx1][SAT_POS_MSG]));
                    break;
                case 5:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it5_[idx1][SAT_POS_MSG]));
                    break;
                default:
                    tmp = VECTOR_ADD_MASK(mask, tmp, VECTOR_SET1(V2C_map_it6_[idx1][SAT_POS_MSG]));
                    break;
                }
#endif

                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(tmp, min1);
                min2 = VECTOR_MIN_2(tmp, vTemp, min2); // the second minimum value of the Lnm

                tmp = VECTOR_ADD_MASK(VECTOR_EQ_MASK(vContr, zero), VECTOR_SIGN(tmp, vContr), tmp);
                tab_vContr[j] = tmp;
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_20; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_20]), _MM_HINT_T0);
#endif
            // WE SATURATE DIRECTLY IN MSG FORMAT
#if OMS_MODE == 0 // simple OMS
            min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2 = VECTOR_MIN(min2, max_msg);
            TYPE cste_1 = VECTOR_MIN(VECTOR_SUB(min2, VECTOR_SET1(offset)), max_msg);
            TYPE cste_2 = VECTOR_MIN(VECTOR_SUB(min1, VECTOR_SET1(offset)), max_msg);
#elif OMS_MODE == 1 // selective OMS

            TYPE min1_offed, min2_offed;
            min1_offed = min1 = VECTOR_MIN(min1, max_msg); // Limiting quantization width
            min2_offed = min2 = VECTOR_MIN(min2, max_msg);

            if (nombre_iterations <= floor_iter_thresh) {
                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min1, factor_2);
                min1_offed = VECTOR_ADD_MASK(msk_lt6_1, min1, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_1 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min1, factor_1);
                min1_offed = VECTOR_ADD_MASK(msk_le1_1, min1_offed, ones); // min_offed = min_offed + 1 = min + 2

                // where checksum is wrong and error_sum < floor_err_count and min < factor_2
                __mmask32 msk_lt6_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LT_MASK(min2, factor_2);
                min2_offed = VECTOR_ADD_MASK(msk_lt6_2, min2, twos); // min_offed = min + 1
                // where checksum is wrong and error_sum < floor_err_count and min <= factor_1
                __mmask32 msk_le1_2 = l_checksum_[cn_count] & l_m_error_sum & VECTOR_LE_MASK(min2, factor_1);
                min2_offed = VECTOR_ADD_MASK(msk_le1_2, min2_offed, ones); // min_offed = min_offed + 1 = min + 2
            }
            TYPE cste_1 = VECTOR_MIN(min2_offed, max_msg);
            TYPE cste_2 = VECTOR_MIN(min1_offed, max_msg);

#endif
#if (DEG_20 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_20)
            for (int j = 0; j < DEG_20; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(temp_vContr[j], v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
            cn_count++;
        }
#endif
    }

    // ?? BF
    __mmask32 hard_llr[_NoVar] = { 0 }; // Hard LLR, highest bit
    __mmask32 hard2_llr[_NoVar] = { 0 }; // Hard LLR, second highest position
    __mmask32 hard_ch[_NoVar] = { 0 }; // Hard LLR, highest bit
    // Record the position to be flipped in a BF iteration (no need to manually reset to 0 in each BF iteration, it will be automatically reset during the traversal of non-zero elements)
    __mmask32 mask_flip_record[_NoVar] = { 0 };
    const TYPE highest2 = VECTOR_SET1(13); // LLR is ANDed with it to get the next highest position
    const TYPE highest2_n = VECTOR_SET1(-13); // LLR is AND with it to get the next highest position
    for (i = 0; i < _NoVar; i++) {
        hard_llr[i] = VECTOR_GT_MASK(var_nodes[i], zero);
        hard2_llr[i] = VECTOR_GE_MASK(var_nodes[i], highest2) | VECTOR_LE_MASK(var_nodes[i], highest2_n);
        hard_ch[i] = hard_llr[i];
    }
    int BFiter = 0;
    __mmask32 t = UINT32_MAX; // Whether it was flipped in the previous iteration, initialized to 1
    TYPE Th = VECTOR_SET1(REGULAR_COL_WEIGHT); // Flip Threshold
    TYPE l0 = zero; // The maximum threshold number
    TYPE l1 = zero; // The number of times the second largest threshold
    const TYPE L0 = VECTOR_SET1(_L0); // Maximum threshold and maximum number of iterations
    const TYPE L1 = VECTOR_SET1(_L1); // Maximum number of iterations for the next largest threshold
    const TYPE alpha = VECTOR_SET1(_alpha);
    const TYPE delta = VECTOR_SET1(_delta); // Th Decrease step length

    while (BFiter < _maxBFiter) {
        //  The weight of each variable node participating in this verification equation is increased by 1, with a maximum of 255 (although the data type is epi8, epu8 is used for calculation)
        TYPE flip_vote[_NoVar];
        for (i = 0; i < _NoVar; i++) {
            flip_vote[i] = zero;
        }
        TYPE max_vote = ones; // Cannot be set to 0, otherwise all correct answers will be flipped

        int cn_count = 0;
        const unsigned short* pCN = PosNoeudsVariable;
        const unsigned short* pCN2 = PosNoeudsVariable; // Used for flip vote writing
        __mmask32 mask_sum; // == : 0; != : 1, initially both are correct
        TYPE error_sum = zero; //Add 1 to the wrong values, up to 255 (although the data type is epi8, epu8 is used for calculation)
        for (i = 0; i < DEG_1_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_1; j++) {
                // > 0 ? (== 0 ? ch > 0)
                mask_sum ^= hard_llr[*pCN];
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADDU_MASK(mask_sum, error_sum, ones);

            for (j = 0; j < DEG_1; j++) {
                flip_vote[*pCN2] = VECTOR_ADDU_MASK(mask_sum, flip_vote[*pCN2], ones);
                max_vote = VECTOR_MAX(max_vote, flip_vote[*pCN2]);
                pCN2++;
            }
        }
#if NB_DEGRES >= 2
        for (i = 0; i < DEG_2_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_2; j++) {
                mask_sum ^= hard_llr[*pCN];
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADDU_MASK(mask_sum, error_sum, ones);

            for (j = 0; j < DEG_2; j++) {
                flip_vote[*pCN2] = VECTOR_ADDU_MASK(mask_sum, flip_vote[*pCN2], ones);
                max_vote = VECTOR_MAX(max_vote, flip_vote[*pCN2]);
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 3
        for (i = 0; i < DEG_3_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_3; j++) {
                mask_sum ^= hard_llr[*pCN];
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADDU_MASK(mask_sum, error_sum, ones);

            for (j = 0; j < DEG_3; j++) {
                flip_vote[*pCN2] = VECTOR_ADDU_MASK(mask_sum, flip_vote[*pCN2], ones);
                max_vote = VECTOR_MAX(max_vote, flip_vote[*pCN2]);
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 4
        for (i = 0; i < DEG_4_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_4; j++) {
                mask_sum ^= hard_llr[*pCN];
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADDU_MASK(mask_sum, error_sum, ones);

            for (j = 0; j < DEG_4; j++) {
                flip_vote[*pCN2] = VECTOR_ADDU_MASK(mask_sum, flip_vote[*pCN2], ones);
                max_vote = VECTOR_MAX(max_vote, flip_vote[*pCN2]);
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 5
        for (i = 0; i < DEG_5_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_5; j++) {
                mask_sum ^= hard_llr[*pCN];
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADDU_MASK(mask_sum, error_sum, ones);

            for (j = 0; j < DEG_5; j++) {
                flip_vote[*pCN2] = VECTOR_ADDU_MASK(mask_sum, flip_vote[*pCN2], ones);
                max_vote = VECTOR_MAX(max_vote, flip_vote[*pCN2]);
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 6
        for (i = 0; i < DEG_6_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_6; j++) {
                mask_sum ^= hard_llr[*pCN];
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADDU_MASK(mask_sum, error_sum, ones);

            for (j = 0; j < DEG_6; j++) {
                flip_vote[*pCN2] = VECTOR_ADDU_MASK(mask_sum, flip_vote[*pCN2], ones);
                max_vote = VECTOR_MAX(max_vote, flip_vote[*pCN2]);
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 7
        for (i = 0; i < DEG_7_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_7; j++) {
                mask_sum ^= hard_llr[*pCN];
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADDU_MASK(mask_sum, error_sum, ones);

            for (j = 0; j < DEG_7; j++) {
                flip_vote[*pCN2] = VECTOR_ADDU_MASK(mask_sum, flip_vote[*pCN2], ones);
                max_vote = VECTOR_MAX(max_vote, flip_vote[*pCN2]);
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 8
        for (i = 0; i < DEG_8_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_8; j++) {
                mask_sum ^= hard_llr[*pCN];
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADDU_MASK(mask_sum, error_sum, ones);

            for (j = 0; j < DEG_8; j++) {
                flip_vote[*pCN2] = VECTOR_ADDU_MASK(mask_sum, flip_vote[*pCN2], ones);
                max_vote = VECTOR_MAX(max_vote, flip_vote[*pCN2]);
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 9
        for (i = 0; i < DEG_9_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_9; j++) {
                mask_sum ^= hard_llr[*pCN];
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADDU_MASK(mask_sum, error_sum, ones);

            for (j = 0; j < DEG_9; j++) {
                flip_vote[*pCN2] = VECTOR_ADDU_MASK(mask_sum, flip_vote[*pCN2], ones);
                max_vote = VECTOR_MAX(max_vote, flip_vote[*pCN2]);
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 10
        for (i = 0; i < DEG_10_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_10; j++) {
                mask_sum ^= hard_llr[*pCN];
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADDU_MASK(mask_sum, error_sum, ones);

            for (j = 0; j < DEG_10; j++) {
                flip_vote[*pCN2] = VECTOR_ADDU_MASK(mask_sum, flip_vote[*pCN2], ones);
                max_vote = VECTOR_MAX(max_vote, flip_vote[*pCN2]);
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 11
        for (i = 0; i < DEG_11_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_11; j++) {
                mask_sum ^= hard_llr[*pCN];
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADDU_MASK(mask_sum, error_sum, ones);

            for (j = 0; j < DEG_11; j++) {
                flip_vote[*pCN2] = VECTOR_ADDU_MASK(mask_sum, flip_vote[*pCN2], ones);
                max_vote = VECTOR_MAX(max_vote, flip_vote[*pCN2]);
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 12
        for (i = 0; i < DEG_12_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_12; j++) {
                mask_sum ^= hard_llr[*pCN];
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADDU_MASK(mask_sum, error_sum, ones);

            for (j = 0; j < DEG_12; j++) {
                flip_vote[*pCN2] = VECTOR_ADDU_MASK(mask_sum, flip_vote[*pCN2], ones);
                max_vote = VECTOR_MAX(max_vote, flip_vote[*pCN2]);
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 13
        for (i = 0; i < DEG_13_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_13; j++) {
                mask_sum ^= hard_llr[*pCN];
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADDU_MASK(mask_sum, error_sum, ones);

            for (j = 0; j < DEG_13; j++) {
                flip_vote[*pCN2] = VECTOR_ADDU_MASK(mask_sum, flip_vote[*pCN2], ones);
                max_vote = VECTOR_MAX(max_vote, flip_vote[*pCN2]);
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 14
        for (i = 0; i < DEG_14_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_14; j++) {
                mask_sum ^= hard_llr[*pCN];
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADDU_MASK(mask_sum, error_sum, ones);

            for (j = 0; j < DEG_14; j++) {
                flip_vote[*pCN2] = VECTOR_ADDU_MASK(mask_sum, flip_vote[*pCN2], ones);
                max_vote = VECTOR_MAX(max_vote, flip_vote[*pCN2]);
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 15
        for (i = 0; i < DEG_15_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_15; j++) {
                mask_sum ^= hard_llr[*pCN];
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADDU_MASK(mask_sum, error_sum, ones);

            for (j = 0; j < DEG_15; j++) {
                flip_vote[*pCN2] = VECTOR_ADDU_MASK(mask_sum, flip_vote[*pCN2], ones);
                max_vote = VECTOR_MAX(max_vote, flip_vote[*pCN2]);
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 16
        for (i = 0; i < DEG_16_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_16; j++) {
                mask_sum ^= hard_llr[*pCN];
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADDU_MASK(mask_sum, error_sum, ones);

            for (j = 0; j < DEG_16; j++) {
                flip_vote[*pCN2] = VECTOR_ADDU_MASK(mask_sum, flip_vote[*pCN2], ones);
                max_vote = VECTOR_MAX(max_vote, flip_vote[*pCN2]);
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 17
        for (i = 0; i < DEG_17_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_17; j++) {
                mask_sum ^= hard_llr[*pCN];
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADDU_MASK(mask_sum, error_sum, ones);

            for (j = 0; j < DEG_17; j++) {
                flip_vote[*pCN2] = VECTOR_ADDU_MASK(mask_sum, flip_vote[*pCN2], ones);
                max_vote = VECTOR_MAX(max_vote, flip_vote[*pCN2]);
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 18
        for (i = 0; i < DEG_18_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_18; j++) {
                mask_sum ^= hard_llr[*pCN];
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADDU_MASK(mask_sum, error_sum, ones);

            for (j = 0; j < DEG_18; j++) {
                flip_vote[*pCN2] = VECTOR_ADDU_MASK(mask_sum, flip_vote[*pCN2], ones);
                max_vote = VECTOR_MAX(max_vote, flip_vote[*pCN2]);
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 19
        for (i = 0; i < DEG_19_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_19; j++) {
                mask_sum ^= hard_llr[*pCN];
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADDU_MASK(mask_sum, error_sum, ones);

            for (j = 0; j < DEG_19; j++) {
                flip_vote[*pCN2] = VECTOR_ADDU_MASK(mask_sum, flip_vote[*pCN2], ones);
                max_vote = VECTOR_MAX(max_vote, flip_vote[*pCN2]);
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 20
        for (i = 0; i < DEG_20_COMPUTATIONS; i++) {
            mask_sum = 0;
            for (j = 0; j < DEG_20; j++) {
                mask_sum ^= hard_llr[*pCN];
                pCN++;
            }
            l_checksum_[cn_count++] = mask_sum;
            error_sum = VECTOR_ADDU_MASK(mask_sum, error_sum, ones);

            for (j = 0; j < DEG_20; j++) {
                flip_vote[*pCN2] = VECTOR_ADDU_MASK(mask_sum, flip_vote[*pCN2], ones);
                max_vote = VECTOR_MAX(max_vote, flip_vote[*pCN2]);
                pCN2++;
            }
        }
#endif
        if (VECTOR_GTU_MASK(error_sum, zero) == 0) { // all the frames are right
            break;
        }
        l_m_error_sum = VECTOR_LTU_MASK(error_sum, VECTOR_SET1(floor_err_count));

        Th = VECTOR_SUB_MASK(~t, Th, delta); // If it didn?t turn over last time, then lower the threshold
        //t is 1 and l0 < L0
        __mmask32 max_Th = t & VECTOR_LT_MASK(l0, L0);
        Th = VECTOR_MOV_MASK(Th, max_Th, VECTOR_SET1(REGULAR_COL_WEIGHT + _alpha));
        l0 = VECTOR_ADD_MASK(max_Th, l0, ones);
        // t is 1 and l0 >= L0 and l1 < L1
        __mmask32 submax_Th = t & (~max_Th) & VECTOR_LT_MASK(l1, L1);
        Th = VECTOR_MOV_MASK(Th, submax_Th, VECTOR_SET1(REGULAR_COL_WEIGHT + _alpha - _delta));
        l1 = VECTOR_ADD_MASK(submax_Th, l1, ones);
        // t is 1 and l0 >= L0 and l1 >= L1
        __mmask32 ssubmax_Th = t & (~max_Th) & (~submax_Th);
        Th = VECTOR_MOV_MASK(Th, ssubmax_Th, VECTOR_SET1(REGULAR_COL_WEIGHT + _alpha - 2 * _delta));
        Th = VECTOR_MAX(Th, ones); // The threshold is at least 1, otherwise the correct answer will be flipped.

        t = 0; // Reset rollover statistics

        cn_count = 0;
        __mmask32 mask_flip = 0; // vote Set to 1 if the threshold is met. The current threshold is min(5, max_vote)
        pCN2 = PosNoeudsVariable;
        for (i = 0; i < DEG_1_COMPUTATIONS; i++) {
            for (j = 0; j < DEG_1; j++) {
                if (VN_weight_[*pCN2] == REGULAR_COL_WEIGHT) {
                    mask_flip = VECTOR_GE_MASK(
                        VECTOR_ADD_MASK(hard_llr[*pCN2] ^ hard_ch[*pCN2], flip_vote[*pCN2], alpha), Th);
                    mask_flip_record[*pCN2] = mask_flip; // The records will be repeated during the traversal, but it will not affect the results.
                    // hard_llr[*pCN2] ^= mask_flip; // Flip if the flip condition is met
                    t |= mask_flip; // The flip position is set to 1
                }
                pCN2++;
            }
        }
#if NB_DEGRES >= 2
        for (i = 0; i < DEG_2_COMPUTATIONS; i++) {
            for (j = 0; j < DEG_2; j++) {
                if (VN_weight_[*pCN2] == REGULAR_COL_WEIGHT) {
                    mask_flip = VECTOR_GE_MASK(
                        VECTOR_ADD_MASK(hard_llr[*pCN2] ^ hard_ch[*pCN2], flip_vote[*pCN2], alpha), Th);
                    mask_flip_record[*pCN2] = mask_flip; // The records will be repeated during the traversal, but it will not affect the results.
                    // hard_llr[*pCN2] ^= mask_flip; // Flip if the flip condition is met
                    t |= mask_flip; // The flip position is set to 1
                }
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 3
        for (i = 0; i < DEG_3_COMPUTATIONS; i++) {
            for (j = 0; j < DEG_3; j++) {
                if (VN_weight_[*pCN2] == REGULAR_COL_WEIGHT) {
                    mask_flip = VECTOR_GE_MASK(
                        VECTOR_ADD_MASK(hard_llr[*pCN2] ^ hard_ch[*pCN2], flip_vote[*pCN2], alpha), Th);
                    mask_flip_record[*pCN2] = mask_flip; // The records will be repeated during the traversal, but it will not affect the results.
                    // hard_llr[*pCN2] ^= mask_flip; // Flip if the flip condition is met
                    t |= mask_flip; // The flip position is set to 1
                }
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 4
        for (i = 0; i < DEG_4_COMPUTATIONS; i++) {
            for (j = 0; j < DEG_4; j++) {
                if (VN_weight_[*pCN2] == REGULAR_COL_WEIGHT) {
                    mask_flip = VECTOR_GE_MASK(
                        VECTOR_ADD_MASK(hard_llr[*pCN2] ^ hard_ch[*pCN2], flip_vote[*pCN2], alpha), Th);
                    mask_flip_record[*pCN2] = mask_flip; // The records will be repeated during the traversal, but it will not affect the results.
                    // hard_llr[*pCN2] ^= mask_flip; // Flip if the flip condition is met
                    t |= mask_flip; // The flip position is set to 1
                }
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 5
        for (i = 0; i < DEG_5_COMPUTATIONS; i++) {
            for (j = 0; j < DEG_5; j++) {
                if (VN_weight_[*pCN2] == REGULAR_COL_WEIGHT) {
                    mask_flip = VECTOR_GE_MASK(
                        VECTOR_ADD_MASK(hard_llr[*pCN2] ^ hard_ch[*pCN2], flip_vote[*pCN2], alpha), Th);
                    mask_flip_record[*pCN2] = mask_flip; // The records will be repeated during the traversal, but it will not affect the results.
                    // hard_llr[*pCN2] ^= mask_flip; // Flip if the flip condition is met
                    t |= mask_flip; // The flip position is set to 1
                }
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 6
        for (i = 0; i < DEG_6_COMPUTATIONS; i++) {
            for (j = 0; j < DEG_6; j++) {
                if (VN_weight_[*pCN2] == REGULAR_COL_WEIGHT) {
                    mask_flip = VECTOR_GE_MASK(
                        VECTOR_ADD_MASK(hard_llr[*pCN2] ^ hard_ch[*pCN2], flip_vote[*pCN2], alpha), Th);
                    mask_flip_record[*pCN2] = mask_flip; // The records will be repeated during the traversal, but it will not affect the results.
                    // hard_llr[*pCN2] ^= mask_flip; // Flip if the flip condition is met
                    t |= mask_flip; // The flip position is set to 1
                }
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 7
        for (i = 0; i < DEG_7_COMPUTATIONS; i++) {
            for (j = 0; j < DEG_7; j++) {
                if (VN_weight_[*pCN2] == REGULAR_COL_WEIGHT) {
                    mask_flip = VECTOR_GE_MASK(
                        VECTOR_ADD_MASK(hard_llr[*pCN2] ^ hard_ch[*pCN2], flip_vote[*pCN2], alpha), Th);
                    mask_flip_record[*pCN2] = mask_flip; // The records will be repeated during the traversal, but it will not affect the results.
                    // hard_llr[*pCN2] ^= mask_flip; // Flip if the flip condition is met
                    t |= mask_flip; // The flip position is set to 1
                }
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 8
        for (i = 0; i < DEG_8_COMPUTATIONS; i++) {
            for (j = 0; j < DEG_8; j++) {
                if (VN_weight_[*pCN2] == REGULAR_COL_WEIGHT) {
                    mask_flip = VECTOR_GE_MASK(
                        VECTOR_ADD_MASK(hard_llr[*pCN2] ^ hard_ch[*pCN2], flip_vote[*pCN2], alpha), Th);
                    mask_flip_record[*pCN2] = mask_flip; // The records will be repeated during the traversal, but it will not affect the results.
                    // hard_llr[*pCN2] ^= mask_flip; // Flip if the flip condition is met
                    t |= mask_flip; // The flip position is set to 1
                }
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 9
        for (i = 0; i < DEG_9_COMPUTATIONS; i++) {
            for (j = 0; j < DEG_9; j++) {
                if (VN_weight_[*pCN2] == REGULAR_COL_WEIGHT) {
                    mask_flip = VECTOR_GE_MASK(
                        VECTOR_ADD_MASK(hard_llr[*pCN2] ^ hard_ch[*pCN2], flip_vote[*pCN2], alpha), Th);
                    mask_flip_record[*pCN2] = mask_flip; // The records will be repeated during the traversal, but it will not affect the results.
                    // hard_llr[*pCN2] ^= mask_flip; // Flip if the flip condition is met
                    t |= mask_flip; // The flip position is set to 1
                }
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 10
        for (i = 0; i < DEG_10_COMPUTATIONS; i++) {
            for (j = 0; j < DEG_10; j++) {
                if (VN_weight_[*pCN2] == REGULAR_COL_WEIGHT) {
                    mask_flip = VECTOR_GE_MASK(
                        VECTOR_ADD_MASK(hard_llr[*pCN2] ^ hard_ch[*pCN2], flip_vote[*pCN2], alpha), Th);
                    mask_flip_record[*pCN2] = mask_flip; // The records will be repeated during the traversal, but it will not affect the results.
                    // hard_llr[*pCN2] ^= mask_flip; // Flip if the flip condition is met
                    t |= mask_flip; // The flip position is set to 1
                }
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 11
        for (i = 0; i < DEG_11_COMPUTATIONS; i++) {
            for (j = 0; j < DEG_11; j++) {
                if (VN_weight_[*pCN2] == REGULAR_COL_WEIGHT) {
                    mask_flip = VECTOR_GE_MASK(
                        VECTOR_ADD_MASK(hard_llr[*pCN2] ^ hard_ch[*pCN2], flip_vote[*pCN2], alpha), Th);
                    mask_flip_record[*pCN2] = mask_flip; // The records will be repeated during the traversal, but it will not affect the results.
                    // hard_llr[*pCN2] ^= mask_flip; // Flip if the flip condition is met
                    t |= mask_flip; // The flip position is set to 1
                }
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 12
        for (i = 0; i < DEG_12_COMPUTATIONS; i++) {
            for (j = 0; j < DEG_12; j++) {
                if (VN_weight_[*pCN2] == REGULAR_COL_WEIGHT) {
                    mask_flip = VECTOR_GE_MASK(
                        VECTOR_ADD_MASK(hard_llr[*pCN2] ^ hard_ch[*pCN2], flip_vote[*pCN2], alpha), Th);
                    mask_flip_record[*pCN2] = mask_flip; // The records will be repeated during the traversal, but it will not affect the results.
                    // hard_llr[*pCN2] ^= mask_flip; // Flip if the flip condition is met
                    t |= mask_flip; // The flip position is set to 1
                }
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 13
        for (i = 0; i < DEG_13_COMPUTATIONS; i++) {
            for (j = 0; j < DEG_13; j++) {
                if (VN_weight_[*pCN2] == REGULAR_COL_WEIGHT) {
                    mask_flip = VECTOR_GE_MASK(
                        VECTOR_ADD_MASK(hard_llr[*pCN2] ^ hard_ch[*pCN2], flip_vote[*pCN2], alpha), Th);
                    mask_flip_record[*pCN2] = mask_flip; // The records will be repeated during the traversal, but it will not affect the results.
                    // hard_llr[*pCN2] ^= mask_flip; // Flip if the flip condition is met
                    t |= mask_flip; // The flip position is set to 1
                }
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 14
        for (i = 0; i < DEG_14_COMPUTATIONS; i++) {
            for (j = 0; j < DEG_14; j++) {
                if (VN_weight_[*pCN2] == REGULAR_COL_WEIGHT) {
                    mask_flip = VECTOR_GE_MASK(
                        VECTOR_ADD_MASK(hard_llr[*pCN2] ^ hard_ch[*pCN2], flip_vote[*pCN2], alpha), Th);
                    mask_flip_record[*pCN2] = mask_flip; // The records will be repeated during the traversal, but it will not affect the results.
                    // hard_llr[*pCN2] ^= mask_flip; // Flip if the flip condition is met
                    t |= mask_flip; // The flip position is set to 1
                }
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 15
        for (i = 0; i < DEG_15_COMPUTATIONS; i++) {
            for (j = 0; j < DEG_15; j++) {
                if (VN_weight_[*pCN2] == REGULAR_COL_WEIGHT) {
                    mask_flip = VECTOR_GE_MASK(
                        VECTOR_ADD_MASK(hard_llr[*pCN2] ^ hard_ch[*pCN2], flip_vote[*pCN2], alpha), Th);
                    mask_flip_record[*pCN2] = mask_flip; // The records will be repeated during the traversal, but it will not affect the results.
                    // hard_llr[*pCN2] ^= mask_flip; // Flip if the flip condition is met
                    t |= mask_flip; // The flip position is set to 1
                }
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 16
        for (i = 0; i < DEG_16_COMPUTATIONS; i++) {
            for (j = 0; j < DEG_16; j++) {
                if (VN_weight_[*pCN2] == REGULAR_COL_WEIGHT) {
                    mask_flip = VECTOR_GE_MASK(
                        VECTOR_ADD_MASK(hard_llr[*pCN2] ^ hard_ch[*pCN2], flip_vote[*pCN2], alpha), Th);
                    mask_flip_record[*pCN2] = mask_flip; // The records will be repeated during the traversal, but it will not affect the results.
                    // hard_llr[*pCN2] ^= mask_flip; // Flip if the flip condition is met
                    t |= mask_flip; // The flip position is set to 1
                }
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 17
        for (i = 0; i < DEG_17_COMPUTATIONS; i++) {
            for (j = 0; j < DEG_17; j++) {
                if (VN_weight_[*pCN2] == REGULAR_COL_WEIGHT) {
                    mask_flip = VECTOR_GE_MASK(
                        VECTOR_ADD_MASK(hard_llr[*pCN2] ^ hard_ch[*pCN2], flip_vote[*pCN2], alpha), Th);
                    mask_flip_record[*pCN2] = mask_flip; // The records will be repeated during the traversal, but it will not affect the results.
                    // hard_llr[*pCN2] ^= mask_flip; // Flip if the flip condition is met
                    t |= mask_flip; // The flip position is set to 1
                }
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 18
        for (i = 0; i < DEG_18_COMPUTATIONS; i++) {
            for (j = 0; j < DEG_18; j++) {
                if (VN_weight_[*pCN2] == REGULAR_COL_WEIGHT) {
                    mask_flip = VECTOR_GE_MASK(
                        VECTOR_ADD_MASK(hard_llr[*pCN2] ^ hard_ch[*pCN2], flip_vote[*pCN2], alpha), Th);
                    mask_flip_record[*pCN2] = mask_flip; // The records will be repeated during the traversal, but it will not affect the results.
                    // hard_llr[*pCN2] ^= mask_flip; // Flip if the flip condition is met
                    t |= mask_flip; // The flip position is set to 1
                }
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 19
        for (i = 0; i < DEG_19_COMPUTATIONS; i++) {
            for (j = 0; j < DEG_19; j++) {
                if (VN_weight_[*pCN2] == REGULAR_COL_WEIGHT) {
                    mask_flip = VECTOR_GE_MASK(
                        VECTOR_ADD_MASK(hard_llr[*pCN2] ^ hard_ch[*pCN2], flip_vote[*pCN2], alpha), Th);
                    mask_flip_record[*pCN2] = mask_flip; // The records will be repeated during the traversal, but it will not affect the results.
                    // hard_llr[*pCN2] ^= mask_flip; // Flip if the flip condition is met
                    t |= mask_flip; // The flip position is set to 1
                }
                pCN2++;
            }
        }
#endif
#if NB_DEGRES >= 20
        for (i = 0; i < DEG_20_COMPUTATIONS; i++) {
            for (j = 0; j < DEG_20; j++) {
                if (VN_weight_[*pCN2] == REGULAR_COL_WEIGHT) {
                    mask_flip = VECTOR_GE_MASK(
                        VECTOR_ADD_MASK(hard_llr[*pCN2] ^ hard_ch[*pCN2], flip_vote[*pCN2], alpha), Th);
                    mask_flip_record[*pCN2] = mask_flip; // The records will be repeated during the traversal, but it will not affect the results.
                    // hard_llr[*pCN2] ^= mask_flip; // Flip if the flip condition is met
                    t |= mask_flip; // The flip position is set to 1
                }
                pCN2++;
            }
        }
#endif
        // 1 is a big jump, 0 is a small jump (each bit should have the same value)
        __mmask32 mask_big_jump = VECTOR_GE_MASK(Th, VECTOR_SET1(REGULAR_COL_WEIGHT));
        // Flipping is performed based on the statistical flipping position
        for (i = 0; i < _NoVar; i++) {
            __mmask32 xor3 = mask_big_jump & mask_flip_record[i]; // XOR position with 3
            // xor3 flips the position that is 1 (i.e. XOR with 1), and the other positions remain unchanged (i.e. XOR with 0)
            // Both bits are XORed with 1, which is equivalent to XORing with 3.
            hard_llr[i] ^= xor3;
            hard2_llr[i] ^= xor3;

            // When the jump is small, if the second highest is 1, only the second highest will be flipped, otherwise only the highest will be flipped
            hard_llr[i] ^= (~mask_big_jump) & mask_flip_record[i] & (~hard2_llr[i]);
            hard2_llr[i] ^= (~mask_big_jump) & mask_flip_record[i] & (hard2_llr[i]);
        }
        // for (i = 0; i < _NoVar; i++) {
        //     hard_llr[i] ^= mask_flip_record[i];
        // }
        BFiter++;
    }

    // Assign the value back to llr, assign 1 to the 1 place and -1 to the 0 place
    for (i = 0; i < _NoVar; i++) {
        var_nodes[i] = VECTOR_MOV_MASK(-ones, hard_llr[i], ones);
    }
    // DATE:20181028
    /*
                The Program is modified for LDPC codes without puncture and shorten and the decoder output
                is the whole codeword,and the statistic result of BER and FER is the whole codeword
                it is different from the Program for 5G platform
        */
    // Deinterleaving
    if ((NOEUD) % 32 == 0) {
        uchar_itranspose_avx((TYPE*)var_nodes, (TYPE*)decodedBits, (NOEUD));
    } else {
        char* ptr = (char*)var_nodes;
        for (i = 0; i < (NOEUD); i += 1) {
            for (int j = 0; j < 32; j += 1) {
                decodedBits[j * (NOEUD) + i] = (ptr[32 * i + j] > 0);
            }
        }
    }
    ////DATE:20190306
    ///*
    // The program is modified for 5G LDPC codes, the decodedBits are the information bits
    // and the FER and BER are compared with the information bits
    //*/
    //
    // if ((NmoinsK - _ShortenBits) % 32 == 0)
    //{
    //	uchar_itranspose_avx((TYPE*)var_nodes, (TYPE*)decodedBits, (NmoinsK - _ShortenBits));
    //}
    // else
    //{
    //	char* ptr = (char*)var_nodes;
    //	for (i = 0; i < (NmoinsK - _ShortenBits); i += 1)
    //	{
    //		for (int j = 0; j < 32; j += 1)
    //		{
    //			decodedBits[j *(NmoinsK - _ShortenBits) + i] = (ptr[32 * i + j] >
    // 0);//varnode0 stores the first information of the first frame, varnode1 stores the first information of the second frame, and varnode32 stores the second information of the first frame.
    //		}
    //	}
    //}
    //
    // return BFiter;
}
