#define INTEL_SSSE3 1
#include <mmintrin.h>
#include <immintrin.h>
#include "gf_optimized.h"

#include "gf_int.h"

#define GF_FIELD_WIDTH (16)
#define GF_FIELD_SIZE (1 << GF_FIELD_WIDTH)


struct gf_w16_logtable_data {
	uint16_t      log_tbl[GF_FIELD_SIZE];
	uint16_t      antilog_tbl[GF_FIELD_SIZE * 2];
	uint16_t      inv_tbl[GF_FIELD_SIZE];
	uint16_t* d_antilog;
};


inline
gf_val_32_t
GF_w16_log_multiply(gf_t* gf, gf_val_32_t a, gf_val_32_t b)
{
	struct gf_w16_logtable_data* ltd;

	ltd = (struct gf_w16_logtable_data*)((gf_internal_t*)gf->scratch)->private;
	return (a == 0 || b == 0) ? 0 : ltd->antilog_tbl[(int)ltd->log_tbl[a] + (int)ltd->log_tbl[b]];
}


inline
gf_val_32_t
GF_w16_log_multiply_by_log(struct gf_w16_logtable_data* ltd, gf_val_32_t a, int log_b)
{
	// Нам гарантировано, что b != 0 , так как это случай обрабатывается отдельно
	// для a=0 эта функция тоже не вызывается
	return ltd->antilog_tbl[(int)ltd->log_tbl[a] + log_b];
}

struct HiLoTableData
{
	__m128i highTable[4];
	__m128i lowTable[4];
};

struct HiLoTableData highLowTable[65536];

void calcHiLoTables(gf_t* gf)
{
	struct gf_w16_logtable_data* ltd;
	ltd = (struct gf_w16_logtable_data*)((gf_internal_t*)gf->scratch)->private;

	for (int val = 0; val < 65536; ++val)
	{
		int logVal = (int)ltd->log_tbl[val];

		uint8_t low[4][16];
		uint8_t high[4][16];

		low[0][0] = 0;
		low[1][0] = 0;
		low[2][0] = 0;
		low[3][0] = 0;

		high[0][0] = 0;
		high[1][0] = 0;
		high[2][0] = 0;
		high[3][0] = 0;

		for (int j = 1; j < 16; j++) {
			for (int i = 0; i < 4; i++) {
				const uint64_t c = (j << (i * 4));
				const uint64_t prod = GF_w16_log_multiply_by_log(ltd, c, logVal);
				low[i][j] = (prod & 0xff);
				high[i][j] = (prod >> 8);
			}
		}

		for (int i = 0; i < 4; i++) {
			highLowTable[val].lowTable[i] = _mm_loadu_si128((__m128i*)low[i]);
			highLowTable[val].highTable[i] = _mm_loadu_si128((__m128i*)high[i]);
		}
	}
}

void gf_multby_one_ex(uint8_t* src, uint8_t* dest)
{
#ifdef   INTEL_SSE2
	__m128i ms, md;
#endif
	uint8_t * send;
	send = src + 512;

	while (src < send) {
		ms = _mm_load_si128((__m128i*)(src));
		md = _mm_load_si128((__m128i*)(dest));
		md = _mm_xor_si128(md, ms);
		_mm_store_si128((__m128i*)(dest), md);
		src += 16;
		dest += 16;
	}
	return;
}

void GF_multiply_region_w32(gf_t* gf, uint8_t* src, uint8_t* dest, gf_val_32_t val)
{
	uint64_t i, j, * s64, * d64, * top64;
	uint64_t a, c, prod;
	uint8_t low[4][16];
	uint8_t high[4][16];
//	gf_region_data rd;

	__m256i  mask, ta, tb, ti, tpl, tph, tta, ttb, shuffler, unshuffler, lmask;

	if (val == 0) { return; }
	if (val == 1) { gf_multby_one_ex(src, dest); return; }



//	gf_set_region_data(&rd, gf, src, dest, bytes, val, 1, 32);
	//	gf_do_initial_region_alignment(&rd);

	struct gf_w16_logtable_data* ltd;
	ltd = (struct gf_w16_logtable_data*)((gf_internal_t*)gf->scratch)->private;

	// Нам гарантировано, что val != 0 , так как это случай обрабатывается отдельно
	int logVal = (int)ltd->log_tbl[val];

	s64 = (uint64_t*)src;
	d64 = (uint64_t*)dest;
	top64 = (uint64_t*)(dest + 512);

	mask = _mm256_set1_epi8(0x0f);
	lmask = _mm256_set1_epi16(0xff);


	const __m256i thigh_0 = _mm256_loadu2_m128i(&highLowTable[val].highTable[0], &highLowTable[val].highTable[0]);
	const __m256i thigh_1 = _mm256_loadu2_m128i(&highLowTable[val].highTable[1], &highLowTable[val].highTable[1]);
	const __m256i thigh_2 = _mm256_loadu2_m128i(&highLowTable[val].highTable[2], &highLowTable[val].highTable[2]);
	const __m256i thigh_3 = _mm256_loadu2_m128i(&highLowTable[val].highTable[3], &highLowTable[val].highTable[3]);
	const __m256i tlow_0 = _mm256_loadu2_m128i(&highLowTable[val].lowTable[0], &highLowTable[val].lowTable[0]);
	const __m256i tlow_1 = _mm256_loadu2_m128i(&highLowTable[val].lowTable[1], &highLowTable[val].lowTable[1]);
	const __m256i tlow_2 = _mm256_loadu2_m128i(&highLowTable[val].lowTable[2], &highLowTable[val].lowTable[2]);
	const __m256i tlow_3 = _mm256_loadu2_m128i(&highLowTable[val].lowTable[3], &highLowTable[val].lowTable[3]);

	while (d64 != top64) {

		ta = _mm256_loadu_epi16(s64);
		tb = _mm256_loadu_epi16(s64 + 4);

		tta = _mm256_srli_epi16(ta, 8);
		ttb = _mm256_srli_epi16(tb, 8);

		tpl = _mm256_and_si256(tb, lmask);
		tph = _mm256_and_si256(ta, lmask);

		tb = _mm256_packus_epi16(tpl, tph);
		ta = _mm256_packus_epi16(ttb, tta);

		ti = _mm256_and_si256(mask, tb);
		tph = _mm256_shuffle_epi8(thigh_0, ti);
		tpl = _mm256_shuffle_epi8(tlow_0, ti);

		tb = _mm256_srli_epi16(tb, 4);
		ti = _mm256_and_si256(mask, tb);
		tpl = _mm256_xor_si256(_mm256_shuffle_epi8(tlow_1, ti), tpl);
		tph = _mm256_xor_si256(_mm256_shuffle_epi8(thigh_1, ti), tph);

		ti = _mm256_and_si256(mask, ta);
		tpl = _mm256_xor_si256(_mm256_shuffle_epi8(tlow_2, ti), tpl);
		tph = _mm256_xor_si256(_mm256_shuffle_epi8(thigh_2, ti), tph);

		ta = _mm256_srli_epi16(ta, 4);
		ti = _mm256_and_si256(mask, ta);
		tpl = _mm256_xor_si256(_mm256_shuffle_epi8(tlow_3, ti), tpl);
		tph = _mm256_xor_si256(_mm256_shuffle_epi8(thigh_3, ti), tph);

		ta = _mm256_unpackhi_epi8(tpl, tph);
		tb = _mm256_unpacklo_epi8(tpl, tph);

		tta = _mm256_loadu_epi16((__m128i*) d64);
		ta = _mm256_xor_si256(ta, tta);
		ttb = _mm256_loadu_epi16((__m128i*) (d64 + 4));
		tb = _mm256_xor_si256(tb, ttb);
		_mm256_storeu_epi16((__m128i*)d64, ta);
		_mm256_storeu_epi16((__m128i*)(d64 + 4), tb);

		d64 += 8;
		s64 += 8;

	}


	//	gf_do_final_region_alignment(&rd);
}