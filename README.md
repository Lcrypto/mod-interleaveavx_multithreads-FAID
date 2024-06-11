# LNS FAID

Source codes for the paper "Optimized LNS-FAID of LDPC Codes: A Hybrid Precision Decoding Approach for 50G-PON".

## Environment Requirements

Requires installation of icc and mkl, supporting AVX512 instruction set.

Before running the code, you need to copy the content in 'Constants/50GPON-dc-original/encode_matrix2.txt' to the back of 'Constants/50GPON-dc-original/encode_matrix1.txt', and then use the combined content to replace the definition of 'GenMatrix' at the last few lines in 'Constants/50GPON-dc-original/Constants_SSE.h'.

### Profile.txt Parameter Explanation

- `DecodeMethod`: 0: NMS, 1: OMS, 2: FAID+DTBF, 3: OMS+BF, 4: OMS+DTBF, 5: FAID+2B1C
  - BF and DTBF iteration numbers need to be modified in the corresponding decode function's `_maxBFiter`.
- `modType`: 1: BPSK, 2: QPSK, 4: 16-QAM, 6: 64-QAM
- `Factor_1`: For NMS, it acts as the normalize parameter for the smallest message; for OMS, it's used as the threshold in the clipping layer strategy.
  - Under NMS, try 26 (`26 / 32 = 0.8125`); usually 1 for OMS.
- `Factor_2`: For NMS, it acts as the normalize parameter for the smallest message; for OMS, it's used as the threshold in the clipping layer strategy.
  - Under NMS, try 26 (`26 / 32 = 0.8125`); usually 6 for OMS.
- `scale`: Quantization coefficient. In the letter, results for QPSK under quantized LNMS, 3-bit LNS-FAID(LUTs), 3-bit LNS-FAID[7], 2-bit LNS-FAID(LUTs) are all 13, for hybrid precision LNS-FAID it's 12.5, and for 2-bit LNS-FAID[7] it's 14.

Different scenarios of LNS-FAID's LUTS are defined in 'CDecoder_FAID.cpp', switched by means like '#define FAID3'.


