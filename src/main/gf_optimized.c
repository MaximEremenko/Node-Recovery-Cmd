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

struct HiLoTableData highLowTable[65536];

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
	// ��� �������������, ��� b != 0 , ��� ��� ��� ������ �������������� ��������
	// ��� a=0 ��� ������� ���� �� ����������
	return ltd->antilog_tbl[(int)ltd->log_tbl[a] + log_b];
}

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

void gf_multby_one_ex_512(__m256i* src, __m256i* dest)
{
	__m256i v;

	for(int i = 0; i < 512 / 32; ++i)
	{
		v = _mm256_loadu_epi16(dest + i);
		v = _mm256_xor_si256(v, _mm256_loadu_epi16(src + i));
		_mm256_storeu_epi16(dest + i, v);
	}
	return;
}

void gf_multby_one_ex_128(__m256i* src, __m256i* dest)
{
	__m256i v;

	for (int i = 0; i < 128 / 32; ++i)
	{
		v = _mm256_loadu_epi16(dest + i);
		v = _mm256_xor_si256(v, _mm256_loadu_epi16(src + i));
		_mm256_storeu_epi16(dest + i, v);
	}
	return;
}


void GF_multiply_region_w32_prepared_512(gf_t* gf, __m256i* src, __m256i* dest, struct HiLoTableData *table)
{
	uint64_t i;

	__m256i  mask, ta, tb, ti, tpl, tph, tta, ttb, lmask;

//	if (val == 0) { return; }
//	if (val == 1) { gf_multby_one_ex(src, dest); return; }

	mask = _mm256_set1_epi8(0x0f);
	lmask = _mm256_set1_epi16(0xff);


	const __m256i thigh_0 = _mm256_loadu2_m128i(&table->highTable[0], &table->highTable[0]);
	const __m256i thigh_1 = _mm256_loadu2_m128i(&table->highTable[1], &table->highTable[1]);
	const __m256i thigh_2 = _mm256_loadu2_m128i(&table->highTable[2], &table->highTable[2]);
	const __m256i thigh_3 = _mm256_loadu2_m128i(&table->highTable[3], &table->highTable[3]);
	const __m256i tlow_0 = _mm256_loadu2_m128i(&table->lowTable[0], &table->lowTable[0]);
	const __m256i tlow_1 = _mm256_loadu2_m128i(&table->lowTable[1], &table->lowTable[1]);
	const __m256i tlow_2 = _mm256_loadu2_m128i(&table->lowTable[2], &table->lowTable[2]);
	const __m256i tlow_3 = _mm256_loadu2_m128i(&table->lowTable[3], &table->lowTable[3]);

	for(int i = 0; i < 512 / 32; i+=2) {

		ta = _mm256_loadu_epi16(src + i);
		tb = _mm256_loadu_epi16(src + i + 1);

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

		tta = _mm256_loadu_epi16(dest + i);
		ta = _mm256_xor_si256(ta, tta);
		ttb = _mm256_loadu_epi16(dest + i + 1);
		tb = _mm256_xor_si256(tb, ttb);
		_mm256_storeu_epi16(dest + i, ta);
		_mm256_storeu_epi16(dest + i + 1, tb);
	}
}

void GF_multiply_region_w32_prepared_128(gf_t* gf, __m256i* src, __m256i* dest, struct HiLoTableData* table)
{
	uint64_t i;

	__m256i  mask, ta, tb, ti, tpl, tph, tta, ttb, lmask;

	//	if (val == 0) { return; }
	//	if (val == 1) { gf_multby_one_ex(src, dest); return; }

	mask = _mm256_set1_epi8(0x0f);
	lmask = _mm256_set1_epi16(0xff);


	const __m256i thigh_0 = _mm256_loadu2_m128i(&table->highTable[0], &table->highTable[0]);
	const __m256i thigh_1 = _mm256_loadu2_m128i(&table->highTable[1], &table->highTable[1]);
	const __m256i thigh_2 = _mm256_loadu2_m128i(&table->highTable[2], &table->highTable[2]);
	const __m256i thigh_3 = _mm256_loadu2_m128i(&table->highTable[3], &table->highTable[3]);
	const __m256i tlow_0 = _mm256_loadu2_m128i(&table->lowTable[0], &table->lowTable[0]);
	const __m256i tlow_1 = _mm256_loadu2_m128i(&table->lowTable[1], &table->lowTable[1]);
	const __m256i tlow_2 = _mm256_loadu2_m128i(&table->lowTable[2], &table->lowTable[2]);
	const __m256i tlow_3 = _mm256_loadu2_m128i(&table->lowTable[3], &table->lowTable[3]);

	for (int i = 0; i < 128 / 32; i += 2) {

		ta = _mm256_loadu_epi16(src + i);
		tb = _mm256_loadu_epi16(src + i + 1);

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

		tta = _mm256_loadu_epi16(dest + i);
		ta = _mm256_xor_si256(ta, tta);
		ttb = _mm256_loadu_epi16(dest + i + 1);
		tb = _mm256_xor_si256(tb, ttb);
		_mm256_storeu_epi16(dest + i, ta);
		_mm256_storeu_epi16(dest + i + 1, tb);
	}
}