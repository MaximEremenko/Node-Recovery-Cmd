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

void GF_multiply_region_w32_prepared(gf_t* gf, __m256i* src, __m256i* dest, gf_val_32_t val, struct HiLoTableData* table);

void calcHiLoTables(gf_t* gf);
