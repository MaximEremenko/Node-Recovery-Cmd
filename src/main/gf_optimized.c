#define INTEL_SSSE3 1
#include <mmintrin.h>
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

__m128i lowTable[65536][4];
__m128i highTable[65536][4];

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
			lowTable[val][i] = _mm_loadu_si128((__m128i*)low[i]);
			highTable[val][i] = _mm_loadu_si128((__m128i*)high[i]);
		}
	}
}

void GF_multiply_region_w32(gf_t* gf, uint8_t* src, uint8_t* dest, gf_val_32_t val, int bytes, int xor)
{
	uint64_t i, j, * s64, * d64, * top64;;
	uint64_t a, c, prod;
	uint8_t low[4][16];
	uint8_t high[4][16];
	gf_region_data rd;

	__m128i  mask, ta, tb, ti, tpl, tph, tta, ttb, shuffler, unshuffler, lmask;

	if (val == 0) { gf_multby_zero(dest, bytes, xor); return; }
	if (val == 1) { gf_multby_one(src, dest, bytes, xor); return; }

	gf_set_region_data(&rd, gf, src, dest, bytes, val, xor, 32);
	gf_do_initial_region_alignment(&rd);

	struct gf_w16_logtable_data* ltd;
	ltd = (struct gf_w16_logtable_data*)((gf_internal_t*)gf->scratch)->private;

	// Нам гарантировано, что val != 0 , так как это случай обрабатывается отдельно
	int logVal = (int)ltd->log_tbl[val];

	// Пропускаем 0 итерацию цикла ниже, так как в её результате гарантировано получается результат 0 
	/*low[0][0] = 0;
	low[1][0] = 0;
	low[2][0] = 0;
	low[3][0] = 0;

	high[0][0] = 0;
	high[1][0] = 0;
	high[2][0] = 0;
	high[3][0] = 0;
	
	for (j = 1; j < 16; j++) {
		for (i = 0; i < 4; i++) {
			c = (j << (i * 4));
			prod = GF_w16_log_multiply_by_log(ltd, c, logVal);
			low[i][j] = (prod & 0xff);
			high[i][j] = (prod >> 8);
		}
	}

	for (i = 0; i < 4; i++) {
		tlow[i] = _mm_loadu_si128((__m128i*)low[i]);
		thigh[i] = _mm_loadu_si128((__m128i*)high[i]);
	}

	for (i = 0; i < 4; i++) {
		tlow[i] = _mm_loadu_si128(&lowTable[val][i]);
		thigh[i] = _mm_loadu_si128(&highTable[val][i]);
	}
	*/


	s64 = (uint64_t*)rd.s_start;
	d64 = (uint64_t*)rd.d_start;
	top64 = (uint64_t*)rd.d_top;

	mask = _mm_set1_epi8(0x0f);
	lmask = _mm_set1_epi16(0xff);

	if (xor) {
		while (d64 != top64) {

			ta = _mm_load_si128((__m128i*) s64);
			tb = _mm_load_si128((__m128i*) (s64 + 2));

			tta = _mm_srli_epi16(ta, 8);
			ttb = _mm_srli_epi16(tb, 8);
			tpl = _mm_and_si128(tb, lmask);
			tph = _mm_and_si128(ta, lmask);

			tb = _mm_packus_epi16(tpl, tph);
			ta = _mm_packus_epi16(ttb, tta);

			ti = _mm_and_si128(mask, tb);
			tph = _mm_shuffle_epi8(_mm_loadu_si128(&highTable[val][0]), ti);
			tpl = _mm_shuffle_epi8(_mm_loadu_si128(&lowTable[val][0]), ti);

			tb = _mm_srli_epi16(tb, 4);
			ti = _mm_and_si128(mask, tb);
			tpl = _mm_xor_si128(_mm_shuffle_epi8(_mm_loadu_si128(&lowTable[val][1]), ti), tpl);
			tph = _mm_xor_si128(_mm_shuffle_epi8(_mm_loadu_si128(&highTable[val][1]), ti), tph);

			ti = _mm_and_si128(mask, ta);
			tpl = _mm_xor_si128(_mm_shuffle_epi8(_mm_loadu_si128(&lowTable[val][2]), ti), tpl);
			tph = _mm_xor_si128(_mm_shuffle_epi8(_mm_loadu_si128(&highTable[val][2]), ti), tph);

			ta = _mm_srli_epi16(ta, 4);
			ti = _mm_and_si128(mask, ta);
			tpl = _mm_xor_si128(_mm_shuffle_epi8(_mm_loadu_si128(&lowTable[val][3]), ti), tpl);
			tph = _mm_xor_si128(_mm_shuffle_epi8(_mm_loadu_si128(&highTable[val][3]), ti), tph);

			ta = _mm_unpackhi_epi8(tpl, tph);
			tb = _mm_unpacklo_epi8(tpl, tph);

			tta = _mm_load_si128((__m128i*) d64);
			ta = _mm_xor_si128(ta, tta);
			ttb = _mm_load_si128((__m128i*) (d64 + 2));
			tb = _mm_xor_si128(tb, ttb);
			_mm_store_si128((__m128i*)d64, ta);
			_mm_store_si128((__m128i*)(d64 + 2), tb);

			d64 += 4;
			s64 += 4;

		}
	}
	else {
		while (d64 != top64) {

			ta = _mm_load_si128((__m128i*) s64);
			tb = _mm_load_si128((__m128i*) (s64 + 2));

			tta = _mm_srli_epi16(ta, 8);
			ttb = _mm_srli_epi16(tb, 8);
			tpl = _mm_and_si128(tb, lmask);
			tph = _mm_and_si128(ta, lmask);

			tb = _mm_packus_epi16(tpl, tph);
			ta = _mm_packus_epi16(ttb, tta);

			ti = _mm_and_si128(mask, tb);
			tph = _mm_shuffle_epi8(_mm_loadu_si128(&highTable[val][0]), ti);
			tpl = _mm_shuffle_epi8(_mm_loadu_si128(&lowTable[val][0]), ti);

			tb = _mm_srli_epi16(tb, 4);
			ti = _mm_and_si128(mask, tb);
			tpl = _mm_xor_si128(_mm_shuffle_epi8(_mm_loadu_si128(&lowTable[val][1]), ti), tpl);
			tph = _mm_xor_si128(_mm_shuffle_epi8(_mm_loadu_si128(&highTable[val][1]), ti), tph);

			ti = _mm_and_si128(mask, ta);
			tpl = _mm_xor_si128(_mm_shuffle_epi8(_mm_loadu_si128(&lowTable[val][2]), ti), tpl);
			tph = _mm_xor_si128(_mm_shuffle_epi8(_mm_loadu_si128(&highTable[val][2]), ti), tph);

			ta = _mm_srli_epi16(ta, 4);
			ti = _mm_and_si128(mask, ta);
			tpl = _mm_xor_si128(_mm_shuffle_epi8(_mm_loadu_si128(&lowTable[val][3]), ti), tpl);
			tph = _mm_xor_si128(_mm_shuffle_epi8(_mm_loadu_si128(&highTable[val][3]), ti), tph);

			ta = _mm_unpackhi_epi8(tpl, tph);
			tb = _mm_unpacklo_epi8(tpl, tph);

			_mm_store_si128((__m128i*)d64, ta);
			_mm_store_si128((__m128i*)(d64 + 2), tb);

			d64 += 4;
			s64 += 4;
		}
	}

	gf_do_final_region_alignment(&rd);
}