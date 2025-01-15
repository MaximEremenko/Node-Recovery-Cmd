#pragma once

#define INTEL_SSSE3 1
#include <mmintrin.h>
#include <immintrin.h>
#include "gf_complete.h"

#define GF_W16_INLINE_MULT_LOG(log, a, b) ((uint32_t)log[a]+(uint32_t)log[b])

struct HiLoTableData
{
	__m128i highTable[4];
	__m128i lowTable[4];
};

extern struct HiLoTableData highLowTable[65536];

void gf_multby_one_ex_128(__m256i* src, __m256i* dest);

void GF_multiply_region_w32_prepared_128(gf_t* gf, __m256i* src, __m256i* dest, struct HiLoTableData* table);
#define GF_multiply_region_w32_dispatch_128(gf, src, dest, val, table)  \
	if(val == 1) gf_multby_one_ex_128(src, dest);	\
	else if(val > 1) GF_multiply_region_w32_prepared_128(gf, src, dest, table)


void gf_multby_one_ex_512(__m256i* src, __m256i* dest);

void GF_multiply_region_w32_prepared_512(gf_t* gf, __m256i* src, __m256i* dest, struct HiLoTableData* table);
#define GF_multiply_region_w32_dispatch_512(gf, src, dest, val, table)  \
	if(val == 1) gf_multby_one_ex_512(src, dest);	\
	else if(val > 1) GF_multiply_region_w32_prepared_512(gf, src, dest, table)


void calcHiLoTables(gf_t* gf);
