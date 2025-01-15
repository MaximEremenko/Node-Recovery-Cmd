#pragma once

#include "gf_complete.h"

#define GF_W16_INLINE_MULT_LOG(log, a, b) ((uint32_t)log[a]+(uint32_t)log[b])

void GF_multiply_region_w32(gf_t* gf, uint8_t* src, uint8_t* dest, gf_val_32_t val);

void calcHiLoTables(gf_t* gf);
