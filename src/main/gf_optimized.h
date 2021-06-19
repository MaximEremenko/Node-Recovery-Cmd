#pragma once

#include "gf_complete.h"


void
GF_multiply_region_w32(gf_t* gf, uint8_t* src, uint8_t* dest, gf_val_32_t val, int bytes, int xor);
