#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <errno.h>
#include <vector>
#include <list>

#ifdef _WIN32
#include "sys/getopt.h"
#else
#include <getopt.h>
#endif

#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include <chrono>
//#ifdef _WIN32
//  #include "../sys/time.h"
//#else
//  #include <sys/time.h>
//#endif



#include "Node.h"
#include "FieldDefs.h"
#include "FieldElement.h"
#include "ConstructionDefs.h"

extern "C"
{
#include "gf_optimized.h"
}

using namespace std;

//=======================================================================================================================================
//=============================================================== Auxiliary functions ===================================================
//=======================================================================================================================================

// Time handling procedures
//void timer_start(double* t)
//{
//    struct timeval  tv;
//
//    gettimeofday(&tv, NULL);
//    *t = (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
//}
//
//double timer_split(const double* t)
//{
//    struct timeval  tv;
//    double  cur_t;
//
//    gettimeofday(&tv, NULL);
//    cur_t = (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
//    return (cur_t - *t);
//}

void ext_vand3_inv(fe_type B[3][3], fe_type a, fe_type b, fe_type c)
{
    fe_type ab = a ^ b;
    fe_type ac = a ^ c;
    fe_type bc = b ^ c;

    fe_type w0 = GF_W16_INLINE_MULT(LOG16, ALOG16, ab, ac);
    fe_type w1 = GF_W16_INLINE_MULT(LOG16, ALOG16, ab, bc);
    fe_type w2 = GF_W16_INLINE_MULT(LOG16, ALOG16, ac, bc);

    B[0][0] = GF_W16_INLINE_MULT(LOG16, ALOG16, b, c);
    B[0][0] = GF_W16_INLINE_DIV(LOG16, DALOG16, B[0][0], w0);
    B[0][1] = GF_W16_INLINE_DIV(LOG16, DALOG16, bc, w0);
    B[0][2] = GF_W16_INLINE_DIV(LOG16, DALOG16, 1, w0);

    B[1][0] = GF_W16_INLINE_MULT(LOG16, ALOG16, a, c);
    B[1][0] = GF_W16_INLINE_DIV(LOG16, DALOG16, B[1][0], w1);
    B[1][1] = GF_W16_INLINE_DIV(LOG16, DALOG16, ac, w1);
    B[1][2] = GF_W16_INLINE_DIV(LOG16, DALOG16, 1, w1);

    B[2][0] = GF_W16_INLINE_MULT(LOG16, ALOG16, a, b);
    B[2][0] = GF_W16_INLINE_DIV(LOG16, DALOG16, B[2][0], w2);
    B[2][1] = GF_W16_INLINE_DIV(LOG16, DALOG16, ab, w2);
    B[2][2] = GF_W16_INLINE_DIV(LOG16, DALOG16, 1, w2);
}

void ext_vand2_inv(fe_type B[2][2], fe_type a, fe_type b)
{
    fe_type w0 = a ^ b;

    B[0][0] = GF_W16_INLINE_DIV(LOG16, DALOG16, b, w0);
    B[0][1] = GF_W16_INLINE_DIV(LOG16, DALOG16, 1, w0);

    B[1][0] = GF_W16_INLINE_DIV(LOG16, DALOG16, a, w0);
    B[1][1] = GF_W16_INLINE_DIV(LOG16, DALOG16, 1, w0);
}


void ext_vand4_inv(fe_type B[4][4], fe_type a, fe_type b, fe_type c, fe_type d)
{
    fe_type ab = a ^ b;
    fe_type ac = a ^ c;
    fe_type ad = a ^ d;
    fe_type bc = b ^ c;
    fe_type bd = b ^ d;
    fe_type cd = c ^ d;
    fe_type abcd = ab ^ cd;

    fe_type a_b = GF_W16_INLINE_MULT(LOG16, ALOG16, a, b);
    fe_type c_d = GF_W16_INLINE_MULT(LOG16, ALOG16, c, d);
    fe_type a_b_c_d = GF_W16_INLINE_MULT(LOG16, ALOG16, a_b, c_d);
    fe_type c_c = GF_W16_INLINE_MULT(LOG16, ALOG16, c, c);
    fe_type d_d = GF_W16_INLINE_MULT(LOG16, ALOG16, d, d);

    fe_type w0 = GF_W16_INLINE_MULT(LOG16, ALOG16, ab, ac);
    w0 = GF_W16_INLINE_MULT(LOG16, ALOG16, w0, ad);

    fe_type w0_a = GF_W16_INLINE_MULT(LOG16, ALOG16, w0, a);

    fe_type w1 = GF_W16_INLINE_MULT(LOG16, ALOG16, ab, bc);
    w1 = GF_W16_INLINE_MULT(LOG16, ALOG16, w1, bd);

    fe_type w1_b = GF_W16_INLINE_MULT(LOG16, ALOG16, w1, b);

    fe_type w2 = GF_W16_INLINE_MULT(LOG16, ALOG16, ac, bc);
    w2 = GF_W16_INLINE_MULT(LOG16, ALOG16, w2, cd);

    fe_type w2_c = GF_W16_INLINE_MULT(LOG16, ALOG16, w2, c);

    fe_type w3 = GF_W16_INLINE_MULT(LOG16, ALOG16, ad, bd);
    w3 = GF_W16_INLINE_MULT(LOG16, ALOG16, w3, cd);

    fe_type w3_d = GF_W16_INLINE_MULT(LOG16, ALOG16, w3, d);

    B[0][0] = GF_W16_INLINE_DIV(LOG16, DALOG16, a_b_c_d, w0_a);

    B[0][1] = GF_W16_INLINE_MULT(LOG16, ALOG16, bd, cd) ^ d_d;
    B[0][1] = GF_W16_INLINE_DIV(LOG16, DALOG16, B[0][1], w0);

    B[0][2] = GF_W16_INLINE_DIV(LOG16, DALOG16, abcd ^ a, w0);

    B[0][3] = GF_W16_INLINE_DIV(LOG16, DALOG16, 1, w0);

    B[1][0] = GF_W16_INLINE_DIV(LOG16, DALOG16, a_b_c_d, w1_b);

    B[1][1] = GF_W16_INLINE_MULT(LOG16, ALOG16, ad, cd) ^ d_d;
    B[1][1] = GF_W16_INLINE_DIV(LOG16, DALOG16, B[1][1], w1);

    B[1][2] = GF_W16_INLINE_DIV(LOG16, DALOG16, abcd ^ b, w1);

    B[1][3] = GF_W16_INLINE_DIV(LOG16, DALOG16, 1, w1);

    B[2][0] = GF_W16_INLINE_DIV(LOG16, DALOG16, a_b_c_d, w2_c);

    B[2][1] = GF_W16_INLINE_MULT(LOG16, ALOG16, ad, bd) ^ d_d;
    B[2][1] = GF_W16_INLINE_DIV(LOG16, DALOG16, B[2][1], w2);

    B[2][2] = GF_W16_INLINE_DIV(LOG16, DALOG16, abcd ^ c, w2);

    B[2][3] = GF_W16_INLINE_DIV(LOG16, DALOG16, 1, w2);

    B[3][0] = GF_W16_INLINE_DIV(LOG16, DALOG16, a_b_c_d, w3_d);

    B[3][1] = GF_W16_INLINE_MULT(LOG16, ALOG16, ac, bc) ^ c_c;
    B[3][1] = GF_W16_INLINE_DIV(LOG16, DALOG16, B[3][1], w3);

    B[3][2] = GF_W16_INLINE_DIV(LOG16, DALOG16, abcd ^ d, w3);

    B[3][3] = GF_W16_INLINE_DIV(LOG16, DALOG16, 1, w3);
}

// Inversion of Vandermond matrix of size 4x4
/*void vand4_inv(FieldElement B[4][4], FieldElement a, FieldElement b, FieldElement c, FieldElement d)
{
    FieldElement ab = a + b;
    FieldElement ac = a + c;
    FieldElement ad = a + d;
    FieldElement bc = b + c;
    FieldElement bd = b + d;
    FieldElement cd = c + d;
    FieldElement abcd = ab + cd;

    FieldElement a_b_c_d = a * b * c * d;
    FieldElement c_c = c * c;
    FieldElement d_d = d * d;

    FieldElement inva = ~a;
    FieldElement invb = ~b;
    FieldElement invc = ~c;
    FieldElement invd = ~d;

    FieldElement w0 = ~(ab * ac * ad);
    FieldElement w1 = ~(ab * bc * bd);
    FieldElement w2 = ~(ac * bc * cd);
    FieldElement w3 = ~(ad * bd * cd);

    B[0][0] = w0 * a_b_c_d * inva;
    B[0][1] = w0 * (bd * cd + d_d);
    B[0][2] = w0 * (abcd + a);
    B[0][3] = w0;

    B[1][0] = w1 * a_b_c_d * invb;
    B[1][1] = w1 * (ad * cd + d_d);
    B[1][2] = w1 * (abcd + b);
    B[1][3] = w1;

    B[2][0] = w2 * a_b_c_d * invc;
    B[2][1] = w2 * (ad * bd + d_d);
    B[2][2] = w2 * (abcd + c);
    B[2][3] = w2;

    B[3][0] = w3 * a_b_c_d * invd;
    B[3][1] = w3 * (ac * bc + c_c);
    B[3][2] = w3 * (abcd + d);
    B[3][3] = w3;
}*/

// Inversion of Vandermond matrix of size 4x4
void vand4_inv2(LogFieldElement B[4][4], FieldElement a, FieldElement b, FieldElement c, FieldElement d)
{
    LogFieldElement abLog = (a + b).toLog();
    LogFieldElement acLog = (a + c).toLog();
    LogFieldElement adLog = (a + d).toLog();
    LogFieldElement bcLog = (b + c).toLog();
    LogFieldElement bdLog = (b + d).toLog();
    LogFieldElement cdLog = (c + d).toLog();
    FieldElement abcd = ((a + b) + (c + d));

    LogFieldElement cLog = c.toLog();
    LogFieldElement dLog = d.toLog();
    LogFieldElement a_b_c_d = a.toLog() * b.toLog() * cLog * dLog;
    LogFieldElement c_c = cLog * cLog;
    LogFieldElement d_d = dLog * dLog;

    LogFieldElement inva = (~a).toLog();
    LogFieldElement invb = (~b).toLog();
    LogFieldElement invc = (~c).toLog();
    LogFieldElement invd = (~d).toLog();

    LogFieldElement w0 = (~(abLog * acLog * adLog).toNormal()).toLog();
    LogFieldElement w1 = (~(abLog * bcLog * bdLog).toNormal()).toLog();
    LogFieldElement w2 = (~(acLog * bcLog * cdLog).toNormal()).toLog();
    LogFieldElement w3 = (~(adLog * bdLog * cdLog).toNormal()).toLog();

    B[0][0] = (w0 * a_b_c_d * inva);
    B[0][1] = (w0 * ((bdLog * cdLog).toNormal() + d_d.toNormal()));
    B[0][2] = (w0 * (abcd + a));
    B[0][3] = w0;

    B[1][0] = (w1 * a_b_c_d * invb);
    B[1][1] = (w1 * ((adLog * cdLog).toNormal() + d_d.toNormal()));
    B[1][2] = (w1 * (abcd + b));
    B[1][3] = w1;

    B[2][0] = ((w2 * a_b_c_d) * invc);
    B[2][1] = (w2 * ((adLog * bdLog).toNormal() + d_d.toNormal()));
    B[2][2] = (w2 * (abcd + c));
    B[2][3] = w2;

    B[3][0] = ((w3 * a_b_c_d) * invd);
    B[3][1] = (w3 * ((acLog * bcLog).toNormal() + c_c.toNormal()));
    B[3][2] = (w3 * (abcd + d));
    B[3][3] = w3;
}


// Inversion of Vandermond matrix of size 3x3
void vand3_inv(LogFieldElement B[3][3], FieldElement a, FieldElement b, FieldElement c)
{
    LogFieldElement ab = (a + b).toLog();
    LogFieldElement ac = (a + c).toLog();
    LogFieldElement bc = (b + c).toLog();

    LogFieldElement w0 = (~((ab * ac).toNormal())).toLog();
    LogFieldElement w1 = (~((ab * bc).toNormal())).toLog();
    LogFieldElement w2 = (~((ac * bc).toNormal())).toLog();

    B[0][0] = w0 * b * c;
    B[0][1] = w0 * bc;
    B[0][2] = w0;

    B[1][0] = w1 * a * c;
    B[1][1] = w1 * ac;
    B[1][2] = w1;

    B[2][0] = w2 * a * b;
    B[2][1] = w2 * ab;
    B[2][2] = w2;
}

// Inversion of Vandermond matrix of size 2x2
void vand2_inv(LogFieldElement B[2][2], FieldElement a, FieldElement b)
{
    LogFieldElement w0 = (~(a + b)).toLog();

    B[0][0] = w0 * b;
    B[0][1] = w0;

    B[1][0] = w0 * a;
    B[1][1] = w0;
}

// Retrieve positions and lambdas ids for current A id and subblock number (for 1 node recovery)
void curr_pos_get(unsigned int* pos, unsigned int* lambda_ids, int curr_matr_A_id, int curr_sub_block)
{
    int coeff[4];
    int temp = 0;

    coeff[0] = (curr_sub_block & 0x3);
    coeff[1] = (curr_sub_block & 0xC) >> 2;
    coeff[2] = (curr_sub_block & 0x30) >> 4;
    coeff[3] = (curr_sub_block & 0xC0) >> 6;

    for (int i = 0; i < 4; ++i)
    {
        if (i < curr_matr_A_id)
        {
            temp = temp + coeff[i] * (1 << (2 * i));
            lambda_ids[i] = coeff[i];
        }
        else
        {
            temp = temp + coeff[i] * (1 << (2 * (i + 1)));
            lambda_ids[i + 1] = coeff[i];
        }
    }

    for (int i = 0; i < 4; ++i)
    {
        pos[i] = temp + i * (1 << (2 * curr_matr_A_id));
    }
}

//=======================================================================================================================================
//=============================================================== Recovery procedures ===================================================
//=======================================================================================================================================

uint8_t sameLambdaData[154][4][512]; //[NodeId][LambdaId][Data]
uint8_t sameLambdaRes[5][4][4][512]; //[AId][LambdaId][Pow][Data]
unsigned int sameLambdaDataCounter[154][4];
unsigned int sameLambdaResCounter[5][4]; //[AId][LambdaId]

//fe_type* pArray2Vector[5][1024]; //[AId][Pointer to Data]

//fe_type resVectors[4][1024]; //[Pow][Data]

fe_type coeffPow[154][4][4]; //[NodeId][LanbdaId][Pow]
fe_type recNodesCoeff[1024][4];

unsigned int availNodes[5][31];
unsigned int availNodesCounter[5];

uint8_t* dataSrc[2400];
uint8_t* resDst[2400];
gf_val_32_t currCoeff[2400];


unsigned int tempNodeIdCounter[2400];
unsigned int tempAIdxCounter[2400];
int tempKCounter[2400];
int tempJCointer[2400];



double ext_recover_4_nodes_core(unsigned int* pNodesToRecoverIdx, Node* pNodes, int lambdasIdx[5], fe_type* pCurrData, double* inneTime1, double* innerTime2)
{
    uint8_t dstTest[512];

    auto start_time = chrono::high_resolution_clock::now();

    fe_type currCol[4] = { 0, 0, 0, 0 };
    LogFieldElement currColLog[4];
    fe_type recData[4] = { 0, 0, 0, 0 };
    LogFieldElement invMatr[4][4];


    unsigned int AIdx = 0, lambdaId = 0, nodeId = 0;
    ;

    for (int i = 0; i < 2400; ++i)
    {
        GF_multiply_region_w32(&GF, dataSrc[i], resDst[i], currCoeff[i]);
//                GF.multiply_region.w32(&GF, dataSrc[i], resDst[i], currCoeff[i], 512, 1);


        /*        memcpy(dstTest, resDst[i], 512);
                GF.multiply_region.w32(&GF, dataSrc[i], resDst[i], currCoeff[i], 512, 1);
                GF_multiply_region_w32(&GF, dataSrc[i], dstTest, currCoeff[i], 512, 1);
                if (memcmp(resDst[i], dstTest, 512) != 0)
                {
                    printf("test failed\n");
                    exit(0);
                }
          */
    }

    auto  innerTime1End = chrono::high_resolution_clock::now();
    auto  inner1Elapsed = chrono::duration_cast<chrono::microseconds>(innerTime1End - start_time);
    *inneTime1 = std::chrono::duration<double>(inner1Elapsed).count();

    auto innerTime2start = chrono::high_resolution_clock::now();

    for (int i = 0; i < 1024; ++i)
    {
        for (int j = 0; j < 4; ++j)
            currCol[j] = 0;

        lambdasIdx[0] = (i & 0x3); //11
        lambdasIdx[1] = (i & 0xC) >> 2; //1100
        lambdasIdx[2] = (i & 0x30) >> 4; //110000
        lambdasIdx[3] = (i & 0xC0) >> 6; //11000000
        lambdasIdx[4] = (i & 0x300) >> 8; //1100000000

        for (AIdx = 0; AIdx < 5; ++AIdx)
        {
            lambdaId = lambdasIdx[AIdx];

            pCurrData = (fe_type*)(&sameLambdaRes[AIdx][lambdaId][0][0]);
            currCol[0] ^= pCurrData[sameLambdaResCounter[AIdx][lambdaId]];

            pCurrData = (fe_type*)(&sameLambdaRes[AIdx][lambdaId][1][0]);
            currCol[1] ^= pCurrData[sameLambdaResCounter[AIdx][lambdaId]];

            pCurrData = (fe_type*)(&sameLambdaRes[AIdx][lambdaId][2][0]);
            currCol[2] ^= pCurrData[sameLambdaResCounter[AIdx][lambdaId]];

            pCurrData = (fe_type*)(&sameLambdaRes[AIdx][lambdaId][3][0]);
            currCol[3] ^= pCurrData[sameLambdaResCounter[AIdx][lambdaId]];


            ++sameLambdaResCounter[AIdx][lambdaId];
        }
        currColLog[0] = LogFieldElement(FieldElement(currCol[0]));
        currColLog[1] = LogFieldElement(FieldElement(currCol[1]));
        currColLog[2] = LogFieldElement(FieldElement(currCol[2]));
        currColLog[3] = LogFieldElement(FieldElement(currCol[3]));

        vand4_inv2(invMatr, recNodesCoeff[i][0], recNodesCoeff[i][1], recNodesCoeff[i][2], recNodesCoeff[i][3]);

        recData[0] = (invMatr[0][0] * currColLog[0]).toNormal().getElement();
        recData[0] ^= (invMatr[0][1] * currColLog[1]).toNormal().getElement();
        recData[0] ^= (invMatr[0][2] * currColLog[2]).toNormal().getElement();
        recData[0] ^= (invMatr[0][3] * currColLog[3]).toNormal().getElement();

        recData[1] = (invMatr[1][0] * currColLog[0]).toNormal().getElement();
        recData[1] ^= (invMatr[1][1] * currColLog[1]).toNormal().getElement();
        recData[1] ^= (invMatr[1][2] * currColLog[2]).toNormal().getElement();
        recData[1] ^= (invMatr[1][3] * currColLog[3]).toNormal().getElement();

        recData[2] = (invMatr[2][0] * currColLog[0]).toNormal().getElement();
        recData[2] ^= (invMatr[2][1] * currColLog[1]).toNormal().getElement();
        recData[2] ^= (invMatr[2][2] * currColLog[2]).toNormal().getElement();
        recData[2] ^= (invMatr[2][3] * currColLog[3]).toNormal().getElement();

        recData[3] = (invMatr[3][0] * currColLog[0]).toNormal().getElement();
        recData[3] ^= (invMatr[3][1] * currColLog[1]).toNormal().getElement();
        recData[3] ^= (invMatr[3][2] * currColLog[2]).toNormal().getElement();
        recData[3] ^= (invMatr[3][3] * currColLog[3]).toNormal().getElement();

        pNodes[pNodesToRecoverIdx[0]].setData(i, recData[0]);
        pNodes[pNodesToRecoverIdx[1]].setData(i, recData[1]);
        pNodes[pNodesToRecoverIdx[2]].setData(i, recData[2]);
        pNodes[pNodesToRecoverIdx[3]].setData(i, recData[3]);

    }
    auto innerTime2End = chrono::high_resolution_clock::now();
    auto inner2Elapsed = chrono::duration_cast<chrono::microseconds>(innerTime2End - innerTime2start);
    *innerTime2 = std::chrono::duration<double>(inner2Elapsed).count();
    // *innerTime2 = timer_split(&start_time2);

    auto  end_time = chrono::high_resolution_clock::now();
    auto  elapsed = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    double elapsed_time = std::chrono::duration<double>(elapsed).count();
    //double elapsed_time = timer_split(&start_time);

    return elapsed_time;
}

// 4 nodes recovery procedure extended
double ext_recover_4_nodes(unsigned int* pNodesToRecoverIdx, Node* pNodes, double* inneTime1, double* innerTime2)
{
    int lambdasIdx[5] = { 0, 0, 0, 0, 0 };
    int i = 0, j = 0, k = 0, v = 0;
    unsigned int AIdx = 0, lambdaId = 0, nodeId = 0;
    fe_type* pCurrData = 0, * pCurrRes = 0;

    for (i = 0; i < 4; ++i)
        for (j = 0; j < 154; ++j)
            sameLambdaDataCounter[j][i] = 0;

    for (i = 0; i < 512; ++i)
        for (j = 0; j < 4; ++j)
            for (k = 0; k < 4; ++k)
                for (v = 0; v < 5; ++v)
                    sameLambdaRes[v][k][j][i] = (uint8_t)0;

    for (i = 0; i < 4; ++i)
        for (j = 0; j < 5; ++j)
            sameLambdaResCounter[j][i] = 0;

    for (i = 0; i < 5; ++i)
        availNodesCounter[i] = 0;

    for (i = 0; i < 1024; ++i)
    {
        lambdasIdx[0] = (i & 0x3); //11
        lambdasIdx[1] = (i & 0xC) >> 2; //1100
        lambdasIdx[2] = (i & 0x30) >> 4; //110000
        lambdasIdx[3] = (i & 0xC0) >> 6; //11000000
        lambdasIdx[4] = (i & 0x300) >> 8; //1100000000


        for (j = 0; j < 154; ++j)
        {
            AIdx = AIdxArray[j];
            lambdaId = lambdasIdx[AIdx];

            pCurrData = (fe_type*)(&sameLambdaData[j][lambdaId][0]);
            pCurrData[sameLambdaDataCounter[j][lambdaId]] = pNodes[j].getData(i).getElement();
            ++sameLambdaDataCounter[j][lambdaId];

            coeffPow[j][lambdaId][0] = (fe_type)1;
            coeffPow[j][lambdaId][1] = GF_W16_INLINE_MULT(LOG16, ALOG16, sigmas[j], lambdas[AIdx][lambdaId]);
            coeffPow[j][lambdaId][2] = GF_W16_INLINE_MULT(LOG16, ALOG16, sigmas_pow2[j], lambdas_pow2[AIdx][lambdaId]);
            coeffPow[j][lambdaId][3] = GF_W16_INLINE_MULT(LOG16, ALOG16, sigmas_pow3[j], lambdas_pow3[AIdx][lambdaId]);

            if (j == pNodesToRecoverIdx[0])
                recNodesCoeff[i][0] = coeffPow[j][lambdaId][1];
            else if (j == pNodesToRecoverIdx[1])
                recNodesCoeff[i][1] = coeffPow[j][lambdaId][1];
            else if (j == pNodesToRecoverIdx[2])
                recNodesCoeff[i][2] = coeffPow[j][lambdaId][1];
            else if (j == pNodesToRecoverIdx[3])
                recNodesCoeff[i][3] = coeffPow[j][lambdaId][1];
        }

        /*for (AIdx = 0; AIdx < 5; ++AIdx)
        {
            lambdaId = lambdasIdx[AIdx];

            pCurrRes = (fe_type*)(&sameLambdaRes[AIdx][lambdaId][0]);
            pArray2Vector[AIdx][i] = &pCurrRes[sameLambdaResCounter[AIdx][lambdaId]];
            ++sameLambdaResCounter[AIdx][lambdaId];
        }*/
    }

    for (i = 0; i < 154; ++i)
    {
        if (i != pNodesToRecoverIdx[0] &&
            i != pNodesToRecoverIdx[1] &&
            i != pNodesToRecoverIdx[2] &&
            i != pNodesToRecoverIdx[3])
        {
            AIdx = AIdxArray[i];

            availNodes[AIdx][availNodesCounter[AIdx]] = i;
            ++availNodesCounter[AIdx];
        }
    }

    unsigned int counter = 0;
    for (j = 0; j < 4; ++j)
    {
        for (AIdx = 0; AIdx < 5; ++AIdx)
        {
            for (i = 0; i < availNodesCounter[AIdx]; ++i)
            {
                nodeId = availNodes[AIdx][i];

                for (k = 0; k < 4; ++k)
                {
                    dataSrc[counter] = &sameLambdaData[nodeId][k][0];
                    resDst[counter] = &sameLambdaRes[AIdx][k][j][0];
                    currCoeff[counter] = coeffPow[nodeId][k][j];
                    ++counter;
                }
            }
        }
    }

//    for (int x = 0; x < 10; ++x)
        //ext_recover_4_nodes_core(pNodesToRecoverIdx, pNodes, lambdasIdx, pCurrData, inneTime1, innerTime2);

    return ext_recover_4_nodes_core(pNodesToRecoverIdx, pNodes, lambdasIdx, pCurrData, inneTime1, innerTime2);
}

// Definitions of global auxiliary variables for 2, 3 and 4 nodes recovery procedures
FieldElement AData[1024][5][31];
FieldElement ACol[1024][5][3];
FieldElement ASigmas[5][31];
FieldElement ASigmasPow2[5][31];
FieldElement ASigmasPow3[5][31];
FieldElement nodesToRecLambdas[1024][4];

// Definitions of global auxiliary variable for 1 node recovery procedure
FieldElement diffAData[256][4][31];
FieldElement diffASigma[4][31];
FieldElement diffASigmaPow2[4][31];
FieldElement diffASigmaPow3[4][31];
FieldElement diffACol[256][4][3];
FieldElement sameAData[256][4][31];
FieldElement sameASigma[31];
FieldElement sameASigmaPow2[31];
FieldElement sameASigmaPow3[31];
FieldElement sameACol[4][3];

double recover_4_nodes_core(unsigned int* pNodesToRecoverIdx, unsigned int ADataNodesNum[5], Node* pNodes)
{
    FieldElement sumCol[4];
    LogFieldElement invMatr[4][4];
    FieldElement currRecData[4];
    FieldElement currSum[5], currSigmaSum[5], currSigmaPow2Sum[5], currSigmaPow3Sum[5];

    auto start_time = chrono::high_resolution_clock::now();
    int i, j, k;

    for (i = 0; i < 1024; ++i)
    {
        for (j = 0; j < 5; ++j)
        {
            currSum[j] = ZERO_ELEMENT;
            currSigmaSum[j] = ZERO_ELEMENT;
            currSigmaPow2Sum[j] = ZERO_ELEMENT;
            currSigmaPow3Sum[j] = ZERO_ELEMENT;

            for (k = 0; k < ADataNodesNum[j]; ++k)
            {
                currSum[j] += AData[i][j][k];
                currSigmaSum[j] += (ASigmas[j][k].toLog() * AData[i][j][k]).toNormal();
                currSigmaPow2Sum[j] += (ASigmasPow2[j][k].toLog() * AData[i][j][k]).toNormal();
                currSigmaPow3Sum[j] += (ASigmasPow3[j][k].toLog() * AData[i][j][k]).toNormal();
            }
        }
        sumCol[0] = currSum[0] + currSum[1] + currSum[2] + currSum[3] + currSum[4];
        sumCol[1] = (currSigmaSum[0].toLog() * ACol[i][0][0]).toNormal()
            + (currSigmaSum[1].toLog() * ACol[i][1][0]).toNormal()
            + (currSigmaSum[2].toLog() * ACol[i][2][0] ).toNormal()
            + (currSigmaSum[3].toLog() * ACol[i][3][0]).toNormal()
            + (currSigmaSum[4].toLog() * ACol[i][4][0]).toNormal();
        sumCol[2] = (currSigmaPow2Sum[0].toLog() * ACol[i][0][1]).toNormal()
            + (currSigmaPow2Sum[1].toLog() * ACol[i][1][1]).toNormal()
            + (currSigmaPow2Sum[2].toLog() * ACol[i][2][1]).toNormal()
            + (currSigmaPow2Sum[3].toLog() * ACol[i][3][1]).toNormal()
            + (currSigmaPow2Sum[4].toLog() * ACol[i][4][1]).toNormal();
        sumCol[3] = (currSigmaPow3Sum[0].toLog() * ACol[i][0][2]).toNormal()
            + (currSigmaPow3Sum[1].toLog() * ACol[i][1][2]).toNormal()
            + (currSigmaPow3Sum[2].toLog() * ACol[i][2][2]).toNormal()
            + (currSigmaPow3Sum[3].toLog() * ACol[i][3][2]).toNormal()
            + (currSigmaPow3Sum[4].toLog() * ACol[i][4][2]).toNormal();

        vand4_inv2(invMatr, nodesToRecLambdas[i][0], nodesToRecLambdas[i][1], nodesToRecLambdas[i][2], nodesToRecLambdas[i][3]);

        currRecData[0] = (invMatr[0][0] * sumCol[0]).toNormal()
            + (invMatr[0][1] * sumCol[1]).toNormal()
            + (invMatr[0][2] * sumCol[2]).toNormal()
            + (invMatr[0][3] * sumCol[3]).toNormal();
        currRecData[1] = (invMatr[1][0] * sumCol[0]).toNormal()
            + (invMatr[1][1] * sumCol[1]).toNormal()
            + (invMatr[1][2] * sumCol[2]).toNormal()
            + (invMatr[1][3] * sumCol[3]).toNormal();
        currRecData[2] = (invMatr[2][0] * sumCol[0]).toNormal()
            + (invMatr[2][1] * sumCol[1]).toNormal()
            + (invMatr[2][2] * sumCol[2]).toNormal()
            + (invMatr[2][3] * sumCol[3]).toNormal();
        currRecData[3] = (invMatr[3][0] * sumCol[0]).toNormal()
            + (invMatr[3][1] * sumCol[1]).toNormal()
            + (invMatr[3][2] * sumCol[2]).toNormal()
            + (invMatr[3][3] * sumCol[3]).toNormal();

        pNodes[pNodesToRecoverIdx[0]].setData(i, currRecData[0]);
        pNodes[pNodesToRecoverIdx[1]].setData(i, currRecData[1]);
        pNodes[pNodesToRecoverIdx[2]].setData(i, currRecData[2]);
        pNodes[pNodesToRecoverIdx[3]].setData(i, currRecData[3]);

    }

    auto  end_time = chrono::high_resolution_clock::now();
    auto  elapsed = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    double elapsed_time = std::chrono::duration<double>(elapsed).count();
    return elapsed_time;

}

// 4 nodes recovery procedure
double recover_4_nodes(unsigned int* pNodesToRecoverIdx, Node* pNodes)
{
    int lambdasIdx[5] = { 0, 0, 0, 0, 0 };
    unsigned int i = 0, j = 0, k = 0;
    int AIdx, currNodeToRec = 0, currLambdaIdx = 0;

    FieldElement currSum[5], currSigmaSum[5], currSigmaPow2Sum[5], currSigmaPow3Sum[5], sumCol[4];
    FieldElement invMatr[4][4], currRecData[4];

    unsigned int ADataNodesNum[5] = { 0, 0, 0, 0, 0 };

    for (i = 0; i < 1024; ++i)
    {
        lambdasIdx[0] = (i & 0x3); //11
        lambdasIdx[1] = (i & 0xC) >> 2; //1100
        lambdasIdx[2] = (i & 0x30) >> 4; //110000
        lambdasIdx[3] = (i & 0xC0) >> 6; //11000000
        lambdasIdx[4] = (i & 0x300) >> 8; //1100000000

        ADataNodesNum[0] = 0;
        ADataNodesNum[1] = 0;
        ADataNodesNum[2] = 0;
        ADataNodesNum[3] = 0;
        ADataNodesNum[4] = 0;

        for (j = 0; j < 154; ++j)
        {
            AIdx = AIdxArray[j];

            if (j != pNodesToRecoverIdx[0] &&
                j != pNodesToRecoverIdx[1] &&
                j != pNodesToRecoverIdx[2] &&
                j != pNodesToRecoverIdx[3])
            {
                AData[i][AIdx][ADataNodesNum[AIdx]] = pNodes[j].getData(i);
                ++ADataNodesNum[AIdx];
            }
        }

        for (AIdx = 0; AIdx < 5; ++AIdx)
        {
            currLambdaIdx = lambdasIdx[AIdx];
            ACol[i][AIdx][0] = lambdas[AIdx][currLambdaIdx];
            ACol[i][AIdx][1] = lambdas_pow2[AIdx][currLambdaIdx];
            ACol[i][AIdx][2] = lambdas_pow3[AIdx][currLambdaIdx];
        }

        currNodeToRec = pNodesToRecoverIdx[0];
        AIdx = AIdxArray[currNodeToRec];
        currLambdaIdx = lambdasIdx[AIdx];
        nodesToRecLambdas[i][0] = (FieldElement(lambdas[AIdx][currLambdaIdx]).toLog() * sigmas[currNodeToRec]).toNormal();
        currNodeToRec = pNodesToRecoverIdx[1];
        AIdx = AIdxArray[currNodeToRec];
        currLambdaIdx = lambdasIdx[AIdx];
        nodesToRecLambdas[i][1] = (FieldElement(lambdas[AIdx][currLambdaIdx]).toLog() * sigmas[currNodeToRec]).toNormal();
        currNodeToRec = pNodesToRecoverIdx[2];
        AIdx = AIdxArray[currNodeToRec];
        currLambdaIdx = lambdasIdx[AIdx];
        nodesToRecLambdas[i][2] = (FieldElement(lambdas[AIdx][currLambdaIdx]).toLog() * sigmas[currNodeToRec]).toNormal();
        currNodeToRec = pNodesToRecoverIdx[3];
        AIdx = AIdxArray[currNodeToRec];
        currLambdaIdx = lambdasIdx[AIdx];
        nodesToRecLambdas[i][3] = (FieldElement(lambdas[AIdx][currLambdaIdx]).toLog() * sigmas[currNodeToRec]).toNormal();
    }

    ADataNodesNum[0] = 0;
    ADataNodesNum[1] = 0;
    ADataNodesNum[2] = 0;
    ADataNodesNum[3] = 0;
    ADataNodesNum[4] = 0;

    for (j = 0; j < 154; ++j)
    {
        AIdx = AIdxArray[j];

        if (j != pNodesToRecoverIdx[0] &&
            j != pNodesToRecoverIdx[1] &&
            j != pNodesToRecoverIdx[2] &&
            j != pNodesToRecoverIdx[3])
        {
            ASigmas[AIdx][ADataNodesNum[AIdx]] = FieldElement(sigmas[j]);
            ASigmasPow2[AIdx][ADataNodesNum[AIdx]] = FieldElement(sigmas_pow2[j]);
            ASigmasPow3[AIdx][ADataNodesNum[AIdx]] = FieldElement(sigmas_pow3[j]);
            ++ADataNodesNum[AIdx];
        }
    }

    // double start_time = 0;
     //timer_start(&start_time);

    return recover_4_nodes_core(pNodesToRecoverIdx, ADataNodesNum, pNodes);
}

// 3 nodes recovery procedure extended
double ext_recover_3_nodes(unsigned int* pNodesToRecoverIdx, Node* pNodes)
{
    int lambdasIdx[5] = { 0, 0, 0, 0, 0 };
    int i = 0, j = 0, k = 0, v = 0;
    unsigned int AIdx = 0, lambdaId = 0, nodeId = 0;
    fe_type* pCurrData = 0, * pCurrRes = 0;
    fe_type currCol[3] = { 0, 0, 0 };
    fe_type recData[3] = { 0, 0, 0 };
    fe_type invMatr[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };

    for (i = 0; i < 4; ++i)
        for (j = 0; j < 154; ++j)
            sameLambdaDataCounter[j][i] = 0;

    for (i = 0; i < 512; ++i)
        for (j = 0; j < 4; ++j)
            for (k = 0; k < 4; ++k)
                for (v = 0; v < 5; ++v)
                    sameLambdaRes[v][k][j][i] = (uint8_t)0;

    for (i = 0; i < 4; ++i)
        for (j = 0; j < 5; ++j)
            sameLambdaResCounter[j][i] = 0;

    for (i = 0; i < 5; ++i)
        availNodesCounter[i] = 0;

    for (i = 0; i < 1024; ++i)
    {
        lambdasIdx[0] = (i & 0x3); //11
        lambdasIdx[1] = (i & 0xC) >> 2; //1100
        lambdasIdx[2] = (i & 0x30) >> 4; //110000
        lambdasIdx[3] = (i & 0xC0) >> 6; //11000000
        lambdasIdx[4] = (i & 0x300) >> 8; //1100000000


        for (j = 0; j < 154; ++j)
        {
            AIdx = AIdxArray[j];
            lambdaId = lambdasIdx[AIdx];

            pCurrData = (fe_type*)(&sameLambdaData[j][lambdaId][0]);
            pCurrData[sameLambdaDataCounter[j][lambdaId]] = pNodes[j].getData(i).getElement();
            ++sameLambdaDataCounter[j][lambdaId];

            coeffPow[j][lambdaId][0] = (fe_type)1;
            coeffPow[j][lambdaId][1] = GF_W16_INLINE_MULT(LOG16, ALOG16, sigmas[j], lambdas[AIdx][lambdaId]);
            coeffPow[j][lambdaId][2] = GF_W16_INLINE_MULT(LOG16, ALOG16, sigmas_pow2[j], lambdas_pow2[AIdx][lambdaId]);

            if (j == pNodesToRecoverIdx[0])
                recNodesCoeff[i][0] = coeffPow[j][lambdaId][1];
            else if (j == pNodesToRecoverIdx[1])
                recNodesCoeff[i][1] = coeffPow[j][lambdaId][1];
            else if (j == pNodesToRecoverIdx[2])
                recNodesCoeff[i][2] = coeffPow[j][lambdaId][1];
        }
    }

    for (i = 0; i < 154; ++i)
    {
        if (i != pNodesToRecoverIdx[0] &&
            i != pNodesToRecoverIdx[1] &&
            i != pNodesToRecoverIdx[2])
        {
            AIdx = AIdxArray[i];

            availNodes[AIdx][availNodesCounter[AIdx]] = i;
            ++availNodesCounter[AIdx];
        }
    }

    unsigned int counter = 0;
    for (j = 0; j < 3; ++j)
    {
        for (AIdx = 0; AIdx < 5; ++AIdx)
        {
            for (i = 0; i < availNodesCounter[AIdx]; ++i)
            {
                nodeId = availNodes[AIdx][i];

                for (k = 0; k < 4; ++k)
                {
                    dataSrc[counter] = &sameLambdaData[nodeId][k][0];
                    resDst[counter] = &sameLambdaRes[AIdx][k][j][0];
                    currCoeff[counter] = coeffPow[nodeId][k][j];

                    tempNodeIdCounter[counter] = nodeId;
                    tempAIdxCounter[counter] = AIdx;
                    tempKCounter[counter] = k;
                    tempJCointer[counter] = j;

                    ++counter;
                }
            }
        }
    }

    auto start_time = chrono::high_resolution_clock::now();

    for (i = 0; i < 1812; ++i)
    {
        GF_multiply_region_w32(&GF, dataSrc[i], resDst[i], currCoeff[i]);
//        GF.multiply_region.w32(&GF, dataSrc[i], resDst[i], currCoeff[i], 512, 1);
    }

    for (i = 0; i < 1024; ++i)
    {
        for (j = 0; j < 3; ++j)
            currCol[j] = 0;

        lambdasIdx[0] = (i & 0x3); //11
        lambdasIdx[1] = (i & 0xC) >> 2; //1100
        lambdasIdx[2] = (i & 0x30) >> 4; //110000
        lambdasIdx[3] = (i & 0xC0) >> 6; //11000000
        lambdasIdx[4] = (i & 0x300) >> 8; //1100000000

        for (AIdx = 0; AIdx < 5; ++AIdx)
        {
            lambdaId = lambdasIdx[AIdx];

            pCurrData = (fe_type*)(&sameLambdaRes[AIdx][lambdaId][0][0]);
            currCol[0] ^= pCurrData[sameLambdaResCounter[AIdx][lambdaId]];

            pCurrData = (fe_type*)(&sameLambdaRes[AIdx][lambdaId][1][0]);
            currCol[1] ^= pCurrData[sameLambdaResCounter[AIdx][lambdaId]];

            pCurrData = (fe_type*)(&sameLambdaRes[AIdx][lambdaId][2][0]);
            currCol[2] ^= pCurrData[sameLambdaResCounter[AIdx][lambdaId]];

            ++sameLambdaResCounter[AIdx][lambdaId];
        }

        ext_vand3_inv(invMatr, recNodesCoeff[i][0], recNodesCoeff[i][1], recNodesCoeff[i][2]);

        recData[0] = GF_W16_INLINE_MULT(LOG16, ALOG16, invMatr[0][0], currCol[0]);
        recData[0] ^= GF_W16_INLINE_MULT(LOG16, ALOG16, invMatr[0][1], currCol[1]);
        recData[0] ^= GF_W16_INLINE_MULT(LOG16, ALOG16, invMatr[0][2], currCol[2]);

        recData[1] = GF_W16_INLINE_MULT(LOG16, ALOG16, invMatr[1][0], currCol[0]);
        recData[1] ^= GF_W16_INLINE_MULT(LOG16, ALOG16, invMatr[1][1], currCol[1]);
        recData[1] ^= GF_W16_INLINE_MULT(LOG16, ALOG16, invMatr[1][2], currCol[2]);

        recData[2] = GF_W16_INLINE_MULT(LOG16, ALOG16, invMatr[2][0], currCol[0]);
        recData[2] ^= GF_W16_INLINE_MULT(LOG16, ALOG16, invMatr[2][1], currCol[1]);
        recData[2] ^= GF_W16_INLINE_MULT(LOG16, ALOG16, invMatr[2][2], currCol[2]);

        pNodes[pNodesToRecoverIdx[0]].setData(i, recData[0]);
        pNodes[pNodesToRecoverIdx[1]].setData(i, recData[1]);
        pNodes[pNodesToRecoverIdx[2]].setData(i, recData[2]);

    }

    auto  end_time = chrono::high_resolution_clock::now();
    auto  elapsed = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    double elapsed_time = std::chrono::duration<double>(elapsed).count();

    return elapsed_time;
}

// 3 nodes recovery procedure
double recover_3_nodes(unsigned int* pNodesToRecoverIdx, Node* pNodes)
{
    int lambdasIdx[5] = { 0, 0, 0, 0, 0 };
    unsigned int i = 0, j = 0, k = 0;
    int AIdx, currNodeToRec = 0, currLambdaIdx = 0;

    FieldElement currSum[5], currSigmaSum[5], currSigmaPow2Sum[5], sumCol[3];
    LogFieldElement invMatr[3][3];
    FieldElement currRecData[3];

    unsigned int ADataNodesNum[5] = { 0, 0, 0, 0, 0 };

    for (i = 0; i < 1024; ++i)
    {
        lambdasIdx[0] = (i & 0x3); //11
        lambdasIdx[1] = (i & 0xC) >> 2; //1100
        lambdasIdx[2] = (i & 0x30) >> 4; //110000
        lambdasIdx[3] = (i & 0xC0) >> 6; //11000000
        lambdasIdx[4] = (i & 0x300) >> 8; //1100000000

        ADataNodesNum[0] = 0;
        ADataNodesNum[1] = 0;
        ADataNodesNum[2] = 0;
        ADataNodesNum[3] = 0;
        ADataNodesNum[4] = 0;

        for (j = 0; j < 154; ++j)
        {
            AIdx = AIdxArray[j];

            if (j != pNodesToRecoverIdx[0] &&
                j != pNodesToRecoverIdx[1] &&
                j != pNodesToRecoverIdx[2])
            {
                AData[i][AIdx][ADataNodesNum[AIdx]] = pNodes[j].getData(i);
                ++ADataNodesNum[AIdx];
            }
        }

        for (AIdx = 0; AIdx < 5; ++AIdx)
        {
            currLambdaIdx = lambdasIdx[AIdx];
            ACol[i][AIdx][0] = lambdas[AIdx][currLambdaIdx];
            ACol[i][AIdx][1] = lambdas_pow2[AIdx][currLambdaIdx];
            ACol[i][AIdx][2] = lambdas_pow3[AIdx][currLambdaIdx];
        }

        currNodeToRec = pNodesToRecoverIdx[0];
        AIdx = AIdxArray[currNodeToRec];
        currLambdaIdx = lambdasIdx[AIdx];
        nodesToRecLambdas[i][0] = (FieldElement(lambdas[AIdx][currLambdaIdx]).toLog() * sigmas[currNodeToRec]).toNormal();
        currNodeToRec = pNodesToRecoverIdx[1];
        AIdx = AIdxArray[currNodeToRec];
        currLambdaIdx = lambdasIdx[AIdx];
        nodesToRecLambdas[i][1] = (FieldElement(lambdas[AIdx][currLambdaIdx]).toLog() * sigmas[currNodeToRec]).toNormal();
        currNodeToRec = pNodesToRecoverIdx[2];
        AIdx = AIdxArray[currNodeToRec];
        currLambdaIdx = lambdasIdx[AIdx];
        nodesToRecLambdas[i][2] = (FieldElement(lambdas[AIdx][currLambdaIdx]).toLog() * sigmas[currNodeToRec]).toNormal();
    }

    ADataNodesNum[0] = 0;
    ADataNodesNum[1] = 0;
    ADataNodesNum[2] = 0;
    ADataNodesNum[3] = 0;
    ADataNodesNum[4] = 0;

    for (j = 0; j < 154; ++j)
    {
        AIdx = AIdxArray[j];

        if (j != pNodesToRecoverIdx[0] &&
            j != pNodesToRecoverIdx[1] &&
            j != pNodesToRecoverIdx[2])
        {
            ASigmas[AIdx][ADataNodesNum[AIdx]] = FieldElement(sigmas[j]);
            ASigmasPow2[AIdx][ADataNodesNum[AIdx]] = FieldElement(sigmas_pow2[j]);
            ++ADataNodesNum[AIdx];
        }
    }

    //double start_time = 0;
    //timer_start(&start_time);
    auto start_time = chrono::high_resolution_clock::now();

    for (i = 0; i < 1024; ++i)
    {
        for (j = 0; j < 5; ++j)
        {
            currSum[j] = ZERO_ELEMENT;
            currSigmaSum[j] = ZERO_ELEMENT;
            currSigmaPow2Sum[j] = ZERO_ELEMENT;

            for (k = 0; k < ADataNodesNum[j]; ++k)
            {
                currSum[j] += AData[i][j][k];
                currSigmaSum[j] += (ASigmas[j][k].toLog() * AData[i][j][k]).toNormal();
                currSigmaPow2Sum[j] += (ASigmasPow2[j][k].toLog() * AData[i][j][k]).toNormal();
            }
        }

        sumCol[0] = currSum[0] + currSum[1] + currSum[2] + currSum[3] + currSum[4];
        sumCol[1] = (currSigmaSum[0].toLog() * ACol[i][0][0]).toNormal()
            + (currSigmaSum[1].toLog() * ACol[i][1][0]).toNormal()
            + (currSigmaSum[2].toLog() * ACol[i][2][0]).toNormal()
            + (currSigmaSum[3].toLog() * ACol[i][3][0]).toNormal()
            + (currSigmaSum[4].toLog() * ACol[i][4][0]).toNormal();
        sumCol[2] = (currSigmaPow2Sum[0].toLog() * ACol[i][0][1]).toNormal()
            + (currSigmaPow2Sum[1].toLog() * ACol[i][1][1]).toNormal()
            + (currSigmaPow2Sum[2].toLog() * ACol[i][2][1]).toNormal()
            + (currSigmaPow2Sum[3].toLog() * ACol[i][3][1]).toNormal()
            + (currSigmaPow2Sum[4].toLog() * ACol[i][4][1]).toNormal();

        vand3_inv(invMatr, nodesToRecLambdas[i][0], nodesToRecLambdas[i][1], nodesToRecLambdas[i][2]);

        currRecData[0] = (invMatr[0][0] * sumCol[0]).toNormal() + (invMatr[0][1] * sumCol[1]).toNormal() + (invMatr[0][2] * sumCol[2]).toNormal();
        currRecData[1] = (invMatr[1][0] * sumCol[0]).toNormal() + (invMatr[1][1] * sumCol[1]).toNormal() + (invMatr[1][2] * sumCol[2]).toNormal();
        currRecData[2] = (invMatr[2][0] * sumCol[0]).toNormal() + (invMatr[2][1] * sumCol[1]).toNormal() + (invMatr[2][2] * sumCol[2]).toNormal();

        pNodes[pNodesToRecoverIdx[0]].setData(i, currRecData[0]);
        pNodes[pNodesToRecoverIdx[1]].setData(i, currRecData[1]);
        pNodes[pNodesToRecoverIdx[2]].setData(i, currRecData[2]);

    }

    auto  end_time = chrono::high_resolution_clock::now();
    auto  elapsed = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    double elapsed_time = std::chrono::duration<double>(elapsed).count();

    return elapsed_time;
}


// 2 nodes recovery procedure extended
double ext_recover_2_nodes(unsigned int* pNodesToRecoverIdx, Node* pNodes)
{
    int lambdasIdx[5] = { 0, 0, 0, 0, 0 };
    int i = 0, j = 0, k = 0, v = 0;
    unsigned int AIdx = 0, lambdaId = 0, nodeId = 0;
    fe_type* pCurrData = 0, * pCurrRes = 0;
    fe_type currCol[2] = { 0, 0 };
    fe_type recData[2] = { 0, 0 };
    fe_type invMatr[2][2] = { {0,0}, {0,0} };

    for (i = 0; i < 4; ++i)
        for (j = 0; j < 154; ++j)
            sameLambdaDataCounter[j][i] = 0;

    for (i = 0; i < 512; ++i)
        for (j = 0; j < 4; ++j)
            for (k = 0; k < 4; ++k)
                for (v = 0; v < 5; ++v)
                    sameLambdaRes[v][k][j][i] = (uint8_t)0;

    for (i = 0; i < 4; ++i)
        for (j = 0; j < 5; ++j)
            sameLambdaResCounter[j][i] = 0;

    for (i = 0; i < 5; ++i)
        availNodesCounter[i] = 0;

    for (i = 0; i < 1024; ++i)
    {
        lambdasIdx[0] = (i & 0x3); //11
        lambdasIdx[1] = (i & 0xC) >> 2; //1100
        lambdasIdx[2] = (i & 0x30) >> 4; //110000
        lambdasIdx[3] = (i & 0xC0) >> 6; //11000000
        lambdasIdx[4] = (i & 0x300) >> 8; //1100000000


        for (j = 0; j < 154; ++j)
        {
            AIdx = AIdxArray[j];
            lambdaId = lambdasIdx[AIdx];

            pCurrData = (fe_type*)(&sameLambdaData[j][lambdaId][0]);
            pCurrData[sameLambdaDataCounter[j][lambdaId]] = pNodes[j].getData(i).getElement();
            ++sameLambdaDataCounter[j][lambdaId];

            coeffPow[j][lambdaId][0] = (fe_type)1;
            coeffPow[j][lambdaId][1] = GF_W16_INLINE_MULT(LOG16, ALOG16, sigmas[j], lambdas[AIdx][lambdaId]);

            if (j == pNodesToRecoverIdx[0])
                recNodesCoeff[i][0] = coeffPow[j][lambdaId][1];
            else if (j == pNodesToRecoverIdx[1])
                recNodesCoeff[i][1] = coeffPow[j][lambdaId][1];
        }
    }

    for (i = 0; i < 154; ++i)
    {
        if (i != pNodesToRecoverIdx[0] &&
            i != pNodesToRecoverIdx[1])
        {
            AIdx = AIdxArray[i];

            availNodes[AIdx][availNodesCounter[AIdx]] = i;
            ++availNodesCounter[AIdx];
        }
    }

    unsigned int counter = 0;
    for (j = 0; j < 2; ++j)
    {
        for (AIdx = 0; AIdx < 5; ++AIdx)
        {
            for (i = 0; i < availNodesCounter[AIdx]; ++i)
            {
                nodeId = availNodes[AIdx][i];

                for (k = 0; k < 4; ++k)
                {
                    dataSrc[counter] = &sameLambdaData[nodeId][k][0];
                    resDst[counter] = &sameLambdaRes[AIdx][k][j][0];
                    currCoeff[counter] = coeffPow[nodeId][k][j];

                    tempNodeIdCounter[counter] = nodeId;
                    tempAIdxCounter[counter] = AIdx;
                    tempKCounter[counter] = k;
                    tempJCointer[counter] = j;

                    ++counter;
                }
            }
        }
    }

    auto start_time = chrono::high_resolution_clock::now();

    for (i = 0; i < 1216; ++i)
    {
        GF_multiply_region_w32(&GF, dataSrc[i], resDst[i], currCoeff[i]);
        //GF.multiply_region.w32(&GF, dataSrc[i], resDst[i], currCoeff[i], 512, 1);
    }

    for (i = 0; i < 1024; ++i)
    {
        for (j = 0; j < 2; ++j)
            currCol[j] = 0;

        lambdasIdx[0] = (i & 0x3); //11
        lambdasIdx[1] = (i & 0xC) >> 2; //1100
        lambdasIdx[2] = (i & 0x30) >> 4; //110000
        lambdasIdx[3] = (i & 0xC0) >> 6; //11000000
        lambdasIdx[4] = (i & 0x300) >> 8; //1100000000

        for (AIdx = 0; AIdx < 5; ++AIdx)
        {
            lambdaId = lambdasIdx[AIdx];

            pCurrData = (fe_type*)(&sameLambdaRes[AIdx][lambdaId][0][0]);
            currCol[0] ^= pCurrData[sameLambdaResCounter[AIdx][lambdaId]];

            pCurrData = (fe_type*)(&sameLambdaRes[AIdx][lambdaId][1][0]);
            currCol[1] ^= pCurrData[sameLambdaResCounter[AIdx][lambdaId]];

            ++sameLambdaResCounter[AIdx][lambdaId];
        }

        ext_vand2_inv(invMatr, recNodesCoeff[i][0], recNodesCoeff[i][1]);

        recData[0] = GF_W16_INLINE_MULT(LOG16, ALOG16, invMatr[0][0], currCol[0]);
        recData[0] ^= GF_W16_INLINE_MULT(LOG16, ALOG16, invMatr[0][1], currCol[1]);

        recData[1] = GF_W16_INLINE_MULT(LOG16, ALOG16, invMatr[1][0], currCol[0]);
        recData[1] ^= GF_W16_INLINE_MULT(LOG16, ALOG16, invMatr[1][1], currCol[1]);

        pNodes[pNodesToRecoverIdx[0]].setData(i, recData[0]);
        pNodes[pNodesToRecoverIdx[1]].setData(i, recData[1]);
    }

    auto  end_time = chrono::high_resolution_clock::now();
    auto  elapsed = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    double elapsed_time = std::chrono::duration<double>(elapsed).count();

    return elapsed_time;
}

// 2 nodes recovery procedure
double recover_2_nodes(unsigned int* pNodesToRecoverIdx, Node* pNodes)
{
    int lambdasIdx[5] = { 0, 0, 0, 0, 0 };
    unsigned int i = 0, j = 0, k = 0;
    int AIdx, currNodeToRec = 0, currLambdaIdx = 0;

    FieldElement currSum[5], currSigmaSum[5], sumCol[2];
    LogFieldElement invMatr[2][2];
    FieldElement currRecData[2];

    unsigned int ADataNodesNum[5] = { 0, 0, 0, 0, 0 };

    for (i = 0; i < 1024; ++i)
    {
        lambdasIdx[0] = (i & 0x3); //11
        lambdasIdx[1] = (i & 0xC) >> 2; //1100
        lambdasIdx[2] = (i & 0x30) >> 4; //110000
        lambdasIdx[3] = (i & 0xC0) >> 6; //11000000
        lambdasIdx[4] = (i & 0x300) >> 8; //1100000000

        ADataNodesNum[0] = 0;
        ADataNodesNum[1] = 0;
        ADataNodesNum[2] = 0;
        ADataNodesNum[3] = 0;
        ADataNodesNum[4] = 0;

        for (j = 0; j < 154; ++j)
        {
            AIdx = AIdxArray[j];

            if (j != pNodesToRecoverIdx[0] &&
                j != pNodesToRecoverIdx[1])
            {
                AData[i][AIdx][ADataNodesNum[AIdx]] = pNodes[j].getData(i);
                ++ADataNodesNum[AIdx];
            }
        }

        for (AIdx = 0; AIdx < 5; ++AIdx)
        {
            currLambdaIdx = lambdasIdx[AIdx];
            ACol[i][AIdx][0] = lambdas[AIdx][currLambdaIdx];
            ACol[i][AIdx][1] = lambdas_pow2[AIdx][currLambdaIdx];
            ACol[i][AIdx][2] = lambdas_pow3[AIdx][currLambdaIdx];
        }

        currNodeToRec = pNodesToRecoverIdx[0];
        AIdx = AIdxArray[currNodeToRec];
        currLambdaIdx = lambdasIdx[AIdx];
        nodesToRecLambdas[i][0] = (FieldElement(lambdas[AIdx][currLambdaIdx]).toLog() * sigmas[currNodeToRec]).toNormal();
        currNodeToRec = pNodesToRecoverIdx[1];
        AIdx = AIdxArray[currNodeToRec];
        currLambdaIdx = lambdasIdx[AIdx];
        nodesToRecLambdas[i][1] = (FieldElement(lambdas[AIdx][currLambdaIdx]).toLog() * sigmas[currNodeToRec]).toNormal();
    }

    ADataNodesNum[0] = 0;
    ADataNodesNum[1] = 0;
    ADataNodesNum[2] = 0;
    ADataNodesNum[3] = 0;
    ADataNodesNum[4] = 0;

    for (j = 0; j < 154; ++j)
    {
        AIdx = AIdxArray[j];

        if (j != pNodesToRecoverIdx[0] &&
            j != pNodesToRecoverIdx[1])
        {
            ASigmas[AIdx][ADataNodesNum[AIdx]] = FieldElement(sigmas[j]);
            ++ADataNodesNum[AIdx];
        }
    }

    //double start_time = 0;
    //timer_start(&start_time);
    auto start_time = chrono::high_resolution_clock::now();

    for (i = 0; i < 1024; ++i)
    {
        for (j = 0; j < 5; ++j)
        {
            currSum[j] = ZERO_ELEMENT;
            currSigmaSum[j] = ZERO_ELEMENT;

            for (k = 0; k < ADataNodesNum[j]; ++k)
            {
                currSum[j] += AData[i][j][k];
                currSigmaSum[j] += (ASigmas[j][k].toLog() * AData[i][j][k]).toNormal();
            }
        }

        sumCol[0] = currSum[0] + currSum[1] + currSum[2] + currSum[3] + currSum[4];
        sumCol[1] = (currSigmaSum[0].toLog() * ACol[i][0][0]).toNormal()
            + (currSigmaSum[1].toLog() * ACol[i][1][0]).toNormal()
            + (currSigmaSum[2].toLog() * ACol[i][2][0]).toNormal()
            + (currSigmaSum[3].toLog() * ACol[i][3][0]).toNormal()
            + (currSigmaSum[4].toLog() * ACol[i][4][0]).toNormal();

        vand2_inv(invMatr, nodesToRecLambdas[i][0], nodesToRecLambdas[i][1]);

        currRecData[0] = (invMatr[0][0] * sumCol[0]).toNormal() + (invMatr[0][1] * sumCol[1]).toNormal();
        currRecData[1] = (invMatr[1][0] * sumCol[0]).toNormal() + (invMatr[1][1] * sumCol[1]).toNormal();

        pNodes[pNodesToRecoverIdx[0]].setData(i, currRecData[0]);
        pNodes[pNodesToRecoverIdx[1]].setData(i, currRecData[1]);

    }

    auto  end_time = chrono::high_resolution_clock::now();
    auto  elapsed = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    double elapsed_time = std::chrono::duration<double>(elapsed).count();

    return elapsed_time; //(double)elapsed_time.count();
}

// 1 node recovery procedure
double recover_1_node(unsigned int* pNodesToRecoverIdx, Node* pNodes)
{
    unsigned int currPos[4], currLambdasIdx[5];
    unsigned int nodeToRecIdx = pNodesToRecoverIdx[0];
    unsigned int nodeToRecAIdx = AIdxArray[nodeToRecIdx];

    unsigned int subBlockIdx = 0, currNodeIdx = 0, currAIdx = 0, tempAIdx = 0, i = 0, j = 0;

    LogFieldElement invMatr[4][4];
    FieldElement nodeToRecSigma = FieldElement(sigmas[nodeToRecIdx]);

    vand4_inv2(invMatr, 
        (nodeToRecSigma.toLog() * lambdas[nodeToRecAIdx][0]).toNormal(),
        (nodeToRecSigma.toLog() * lambdas[nodeToRecAIdx][1]).toNormal(),
        (nodeToRecSigma.toLog() * lambdas[nodeToRecAIdx][2]).toNormal(),
        (nodeToRecSigma.toLog() * lambdas[nodeToRecAIdx][3]).toNormal());

    FieldElement diffAPartSum[4] = { ZERO_ELEMENT, ZERO_ELEMENT, ZERO_ELEMENT, ZERO_ELEMENT };
    FieldElement sameAPartSum[4] = { ZERO_ELEMENT, ZERO_ELEMENT, ZERO_ELEMENT, ZERO_ELEMENT };
    FieldElement ASum[4] = { ZERO_ELEMENT, ZERO_ELEMENT, ZERO_ELEMENT, ZERO_ELEMENT };

    unsigned int diffANodeNum[4] = { 0, 0, 0, 0 };
    unsigned int sameANodeNum = 0;

    unsigned int posArray[256][4] = { 0 };

    FieldElement recData[4];

    for (subBlockIdx = 0; subBlockIdx < 256; ++subBlockIdx)
    {
        curr_pos_get(currPos, currLambdasIdx, nodeToRecAIdx, subBlockIdx);

        posArray[subBlockIdx][0] = currPos[0];
        posArray[subBlockIdx][1] = currPos[1];
        posArray[subBlockIdx][2] = currPos[2];
        posArray[subBlockIdx][3] = currPos[3];

        sameANodeNum = 0;
        diffANodeNum[0] = 0;
        diffANodeNum[1] = 0;
        diffANodeNum[2] = 0;
        diffANodeNum[3] = 0;

        for (currNodeIdx = 0; currNodeIdx < 154; ++currNodeIdx)
        {
            if (currNodeIdx != nodeToRecIdx)
            {
                currAIdx = AIdxArray[currNodeIdx];
                if (currAIdx == nodeToRecAIdx)
                {
                    for (i = 0; i < 4; ++i)
                    {
                        sameAData[subBlockIdx][i][sameANodeNum] = pNodes[currNodeIdx].getData(currPos[i]);
                    }
                    ++sameANodeNum;
                }
                else
                {
                    tempAIdx = (currAIdx > nodeToRecAIdx) ? currAIdx - 1 : currAIdx;
                    diffAData[subBlockIdx][tempAIdx][diffANodeNum[tempAIdx]] = pNodes[currNodeIdx].getDataSum(currPos);
                    ++diffANodeNum[tempAIdx];
                }
            }
        }

        for (currAIdx = 0; currAIdx < 5; ++currAIdx)
        {
            if (currAIdx != nodeToRecAIdx)
            {
                tempAIdx = (currAIdx > nodeToRecAIdx) ? currAIdx - 1 : currAIdx;
                diffACol[subBlockIdx][tempAIdx][0] = lambdas[currAIdx][currLambdasIdx[currAIdx]];
                diffACol[subBlockIdx][tempAIdx][1] = lambdas_pow2[currAIdx][currLambdasIdx[currAIdx]];
                diffACol[subBlockIdx][tempAIdx][2] = lambdas_pow3[currAIdx][currLambdasIdx[currAIdx]];
            }
        }
    }

    sameANodeNum = 0;
    diffANodeNum[0] = 0;
    diffANodeNum[1] = 0;
    diffANodeNum[2] = 0;
    diffANodeNum[3] = 0;

    for (currNodeIdx = 0; currNodeIdx < 154; ++currNodeIdx)
    {
        if (currNodeIdx != nodeToRecIdx)
        {
            currAIdx = AIdxArray[currNodeIdx];

            if (currAIdx == nodeToRecAIdx)
            {
                sameASigma[sameANodeNum] = FieldElement(sigmas[currNodeIdx]);
                sameASigmaPow2[sameANodeNum] = FieldElement(sigmas_pow2[currNodeIdx]);
                sameASigmaPow3[sameANodeNum] = FieldElement(sigmas_pow3[currNodeIdx]);
                ++sameANodeNum;
            }
            else
            {
                tempAIdx = (currAIdx > nodeToRecAIdx) ? currAIdx - 1 : currAIdx;
                diffASigma[tempAIdx][diffANodeNum[tempAIdx]] = FieldElement(sigmas[currNodeIdx]);
                diffASigmaPow2[tempAIdx][diffANodeNum[tempAIdx]] = FieldElement(sigmas_pow2[currNodeIdx]);
                diffASigmaPow3[tempAIdx][diffANodeNum[tempAIdx]] = FieldElement(sigmas_pow3[currNodeIdx]);
                ++diffANodeNum[tempAIdx];
            }
        }
    }

    for (i = 0; i < 4; ++i)
    {
        sameACol[i][0] = FieldElement(lambdas[nodeToRecAIdx][i]);
        sameACol[i][1] = FieldElement(lambdas_pow2[nodeToRecAIdx][i]);
        sameACol[i][2] = FieldElement(lambdas_pow3[nodeToRecAIdx][i]);
    }

    //double start_time = 0;
    //timer_start(&start_time);
    auto start_time = chrono::high_resolution_clock::now();

    for (subBlockIdx = 0; subBlockIdx < 256; ++subBlockIdx)
    {
        for (i = 0; i < 4; ++i)
        {
            ASum[i] = ZERO_ELEMENT;
            recData[i] = ZERO_ELEMENT;
        }

        for (i = 0; i < 4; ++i)
        {
            diffAPartSum[0] = ZERO_ELEMENT;
            diffAPartSum[1] = ZERO_ELEMENT;
            diffAPartSum[2] = ZERO_ELEMENT;
            diffAPartSum[3] = ZERO_ELEMENT;
            for (j = 0; j < diffANodeNum[i]; ++j)
            {
                diffAPartSum[0] += diffAData[subBlockIdx][i][j];
                diffAPartSum[1] += (diffASigma[i][j].toLog() * diffAData[subBlockIdx][i][j]).toNormal();
                diffAPartSum[2] += (diffASigmaPow2[i][j].toLog() * diffAData[subBlockIdx][i][j]).toNormal();
                diffAPartSum[3] += (diffASigmaPow3[i][j].toLog() * diffAData[subBlockIdx][i][j]).toNormal();
            }

            ASum[0] += diffAPartSum[0];
            ASum[1] += (diffAPartSum[1].toLog() * diffACol[subBlockIdx][i][0]).toNormal();
            ASum[2] += (diffAPartSum[2].toLog() * diffACol[subBlockIdx][i][1]).toNormal();
            ASum[3] += (diffAPartSum[3].toLog() * diffACol[subBlockIdx][i][2]).toNormal();
        }

        for (i = 0; i < sameANodeNum; ++i)
        {
            sameAPartSum[0] = ZERO_ELEMENT;
            sameAPartSum[1] = ZERO_ELEMENT;
            sameAPartSum[2] = ZERO_ELEMENT;
            sameAPartSum[3] = ZERO_ELEMENT;
            for (j = 0; j < 4; ++j)
            {
                sameAPartSum[0] += sameAData[subBlockIdx][j][i];
                sameAPartSum[1] += (sameAData[subBlockIdx][j][i].toLog() * sameACol[j][0]).toNormal();
                sameAPartSum[2] += (sameAData[subBlockIdx][j][i].toLog() * sameACol[j][1]).toNormal();
                sameAPartSum[3] += (sameAData[subBlockIdx][j][i].toLog() * sameACol[j][2]).toNormal();;
            }
            ASum[0] += sameAPartSum[0];
            ASum[1] += (sameASigma[i].toLog() * sameAPartSum[1]).toNormal();
            ASum[2] += (sameASigmaPow2[i].toLog() * sameAPartSum[2]).toNormal();
            ASum[3] += (sameASigmaPow3[i].toLog() * sameAPartSum[3]).toNormal();
        }

        for (i = 0; i < 4; ++i)
        {
            for (j = 0; j < 4; ++j)
            {
                recData[i] += (invMatr[i][j] * ASum[j]).toNormal();
            }
        }

        pNodes[nodeToRecIdx].setData(posArray[subBlockIdx][0], recData[0]);
        pNodes[nodeToRecIdx].setData(posArray[subBlockIdx][1], recData[1]);
        pNodes[nodeToRecIdx].setData(posArray[subBlockIdx][2], recData[2]);
        pNodes[nodeToRecIdx].setData(posArray[subBlockIdx][3], recData[3]);
    }

    auto  end_time = chrono::high_resolution_clock::now();
    auto  elapsed = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    double elapsed_time = std::chrono::duration<double>(elapsed).count();

    return elapsed_time; //(double)elapsed_time.count();
}

//=======================================================================================================================================
//================================================================== Main function ======================================================
//=======================================================================================================================================


double elapsed_time[10000] = { 0 };
double inner_time1[10000] = { 0 };
double inner_time2[10000] = { 0 };

FieldElement testData[1024];
Node pNodes[154];
Node pTestNodes[4];

uint8_t src[2048];
uint8_t dst[2048];
fe_type a[1000];

int main()
{
    if (!gf_init_easy(&GF, 16))
    {
        printf("GF initialization FAILED!\n");
        exit(1);
    }

    calcHiLoTables(&GF);

    LOG16 = gf_w16_get_log_table(&GF);
    ALOG16 = gf_w16_get_mult_alog_table(&GF);
    DALOG16 = gf_w16_get_div_alog_table(&GF);

    LOG16_UI = new unsigned int[0x10000];

    for (int i = 0; i <= 0xFFFF; ++i)
    {
        LOG16_UI[i] = LOG16[i];
    }
    LOG16_UI[0] = ZeroInit;

    for (int i = 0, e = 0; i <= 0x1FFFF && e < 10; ++i)
    {
        if (LOG16[ALOG16[i]] != LogFieldElement::correctLog(i))
        {
            printf("Logspace correction error: i = %x, ALOG[i] = %x, LOG[ALOG[i]] = %x != corrected(i) = %x\n",
                i,
                ALOG16[i],
                LOG16[ALOG16[i]],
                LogFieldElement::correctLog(i));
            ++e;
        }
    }

    fe_type val = 1;

    FieldElement nzData = val;

    //data = val;
    //nzData = val;


    for (; val < 0xFFFF; ++val)
    {
        nzData = val;

        fe_type square = GF_W16_INLINE_MULT(LOG16, ALOG16, val, val);
        fe_type cube = GF_W16_INLINE_MULT(LOG16, ALOG16, square, val);

        if (cube != (nzData.toLog() * nzData.toLog() * nzData.toLog()).toNormal().getElement())
        {
            printf("Cube error for 0x%X\n", val);
            break;
        }
    }

//    printf("val = %X\n", val);
    
    FieldElement nzsq = (nzData.toLog() * nzData.toLog() * nzData.toLog()).toNormal().getElement();

//    printf("nzFE: %0X ^3  = %0X\n", nzData.getElement(), nzsq.getElement());



    /*fe_type* currData = (fe_type*)(&src[0]);

    for (unsigned int i = 0; i < 1024; ++i)
    {
        currData[i] = (fe_type)(rand() % (1 << 16));
    }

    for (unsigned int i = 0; i < 1000; ++i)
    {
        a[i] = (fe_type)(rand() % (1 << 16));
    }

    double start_time = 0;
    timer_start(&start_time);

    for (unsigned int i = 0; i < 1000; ++i)
    {
        GF.multiply_region.w32(&GF, &src[0], &dst[0], a[i], 2048, 1);
    }

    double curr_time = timer_split(&start_time);

    currData = (fe_type*)(&dst[0]);
    fe_type res = currData[0];
    for (unsigned int i = 1; i < 1024; ++i)
    {
        res ^= currData[i];
    }

    printf("Time: %f usec\n", curr_time * 1e6);
    printf("Speed: %g Mb/s\n", ((double)(16 * 1000)) / (1024*curr_time));
    printf("Res.: 0x%04x\n", res);

    exit(0);*/

    for (unsigned int i = 0; i < 154; ++i)
    {
        for (unsigned int j = 0; j < 1024; ++j)
        {
            testData[j] = FieldElement((fe_type)(rand() % (1 << 16)));
        }
        pNodes[i] = Node(FieldElement(sigmas[i]), testData);
    }

    unsigned int pNodesToRecoverIdx[4] = { 150, 151, 152, 153 };
    //recover_4_nodes(pNodesToRecoverIdx, pNodes);
    ext_recover_4_nodes(pNodesToRecoverIdx, pNodes, &inner_time1[0], &inner_time2[0]);
    printf("Encoding: Done!\n");

    printf("Code word check: ");
    int failFlag = 0;
    FieldElement S;
    int lambdasIdx[5] = { 0, 0, 0, 0, 0 };
    for (int j = 0; j < 1024; ++j)
    {
        S = ZERO_ELEMENT;
        for (int i = 0; i < 154; ++i)
        {
            S = S + pNodes[i].getData(j);
        }

        if (S.getElement() != 0)
            failFlag = 1;

        lambdasIdx[0] = (j & 0x3); //11
        lambdasIdx[1] = (j & 0xC) >> 2; //1100
        lambdasIdx[2] = (j & 0x30) >> 4; //110000
        lambdasIdx[3] = (j & 0xC0) >> 6; //11000000
        lambdasIdx[4] = (j & 0x300) >> 8; //1100000000

        for (int i = 0; i < 154; ++i)
        {
            int AIdx = AIdxArray[i];
            int currLambdaIdx = lambdasIdx[AIdx];

            S = S + (FieldElement(lambdas[AIdx][currLambdaIdx]).toLog() * FieldElement(sigmas[i]) * pNodes[i].getData(j)).toNormal();
        }

        if (S.getElement() != 0)
            failFlag = 1;

        for (int i = 0; i < 154; ++i)
        {
            int AIdx = AIdxArray[i];
            int currLambdaIdx = lambdasIdx[AIdx];

            S = S + (FieldElement(lambdas_pow2[AIdx][currLambdaIdx]).toLog() * FieldElement(sigmas_pow2[i]) * pNodes[i].getData(j)).toNormal();
        }

        if (S.getElement() != 0)
            failFlag = 1;

        for (int i = 0; i < 154; ++i)
        {
            int AIdx = AIdxArray[i];
            int currLambdaIdx = lambdasIdx[AIdx];

            S = S + (FieldElement(lambdas_pow3[AIdx][currLambdaIdx]).toLog() * FieldElement(sigmas_pow3[i]) * pNodes[i].getData(j)).toNormal();
        }

        if (S.getElement() != 0)
            failFlag = 1;
    }
    if (failFlag)
        printf("FAILED!\n");
    else
        printf("PASS!\n");

    for (int i = 0; i < 4; ++i)
    {
        pTestNodes[i] = pNodes[pNodesToRecoverIdx[i]];
        for (int j = 0; j < 1024; ++j)
            pNodes[pNodesToRecoverIdx[i]].setData(j, 0);
    }

    printf("4 nodes ext recovery check: ");
    double time1, time2;
    ext_recover_4_nodes(pNodesToRecoverIdx, pNodes, &time1, &time2);

    failFlag = 0;
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 1024; ++j)
        {
            if (pTestNodes[i].getData(j).getElement() != pNodes[pNodesToRecoverIdx[i]].getData(j).getElement())
                failFlag = 1;
        }
    }
    if (failFlag)
        printf("FAILED!\n");
    else
        printf("PASS!\n");



    for (int i = 0; i < 4; ++i)
    {
        pTestNodes[i] = pNodes[pNodesToRecoverIdx[i]];
        for (int j = 0; j < 1024; ++j)
            pNodes[pNodesToRecoverIdx[i]].setData(j, 0);
    }

    printf("4 nodes recovery check: ");
    recover_4_nodes(pNodesToRecoverIdx, pNodes);
    failFlag = 0;
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 1024; ++j)
        {
            if (pTestNodes[i].getData(j).getElement() != pNodes[pNodesToRecoverIdx[i]].getData(j).getElement())
                failFlag = 1;
        }
    }
    if (failFlag)
        printf("FAILED!\n");
    else
        printf("PASS!\n");

    for (int i = 0; i < 4; ++i)
        pTestNodes[i] = pNodes[pNodesToRecoverIdx[i]];

    for (int i = 0; i < 3; ++i)
    {
        pTestNodes[i] = pNodes[pNodesToRecoverIdx[i]];
        for (int j = 0; j < 1024; ++j)
            pNodes[pNodesToRecoverIdx[i]].setData(j, 0);
    }


    printf("3 nodes ext recovery check: ");
    ext_recover_3_nodes(pNodesToRecoverIdx, pNodes);
    failFlag = 0;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 1024; ++j)
        {
            if (pTestNodes[i].getData(j).getElement() != pNodes[pNodesToRecoverIdx[i]].getData(j).getElement())
                failFlag = 1;
        }
    }
    if (failFlag)
        printf("FAILED!\n");
    else
        printf("PASS!\n");

    for (int i = 0; i < 2; ++i)
    {
        pTestNodes[i] = pNodes[pNodesToRecoverIdx[i]];
        for (int j = 0; j < 1024; ++j)
            pNodes[pNodesToRecoverIdx[i]].setData(j, 0);
    }


    printf("2 nodes ext recovery check: ");
    ext_recover_2_nodes(pNodesToRecoverIdx, pNodes);
    failFlag = 0;
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 1024; ++j)
        {
            if (pTestNodes[i].getData(j).getElement() != pNodes[pNodesToRecoverIdx[i]].getData(j).getElement())
                failFlag = 1;
        }
    }
    if (failFlag)
        printf("FAILED!\n");
    else
        printf("PASS!\n");

    for (int i = 0; i < 1; ++i)
    {
        pTestNodes[i] = pNodes[pNodesToRecoverIdx[i]];
        for (int j = 0; j < 1024; ++j)
            pNodes[pNodesToRecoverIdx[i]].setData(j, 0);
    }


    printf("1 node recovery check: ");
    recover_1_node(pNodesToRecoverIdx, pNodes);
    failFlag = 0;
    for (int i = 0; i < 1; ++i)
    {
        for (int j = 0; j < 1024; ++j)
        {
            if (pTestNodes[i].getData(j).getElement() != pNodes[pNodesToRecoverIdx[i]].getData(j).getElement())
                failFlag = 1;
        }
    }
    if (failFlag)
        printf("FAILED!\n");
    else
        printf("PASS!\n");

    unsigned int nodesToRecoverNum = 1;

    int testsNum = 1000;

    unsigned int pTestNodesToRecoverIdx[4] = { 1, 2, 3, 4 };
    unsigned int rndId = 0;

    printf("\n==Reconstruction speed test:\n");
    for (nodesToRecoverNum = 2; nodesToRecoverNum < 5; ++nodesToRecoverNum)
    {
        printf("---- %d nodes reconstruction:\n", nodesToRecoverNum);
        for (int tests = 0; tests < testsNum; ++tests)
        {
            printf("# tests performed: %d\r", tests + 1);

            for (unsigned int i = 0; i < 154; ++i)
            {
                for (unsigned int j = 0; j < 1024; ++j)
                {
                    testData[j] = FieldElement((fe_type)(rand() % (1 << 16)));
                }
                pNodes[i] = Node(FieldElement(sigmas[i]), testData);
            }
            ext_recover_4_nodes(pNodesToRecoverIdx, pNodes, &(inner_time1[tests]), &(inner_time2[tests]));

            rndId = rand() % 154;
            for (unsigned int i = 0; i < 4; ++i)
                pTestNodesToRecoverIdx[i] = (rndId + i) % 154;

            switch (nodesToRecoverNum)
            {
            case 1:
                elapsed_time[tests] = recover_1_node(pTestNodesToRecoverIdx, pNodes);
                break;
            case 2:
                elapsed_time[tests] = ext_recover_2_nodes(pTestNodesToRecoverIdx, pNodes);
                break;
            case 3:
                elapsed_time[tests] = ext_recover_3_nodes(pTestNodesToRecoverIdx, pNodes);
                break;
            case 4:
                elapsed_time[tests] = ext_recover_4_nodes(pTestNodesToRecoverIdx, pNodes, &(inner_time1[tests]), &(inner_time2[tests]));
                break;
            default:
                printf("Invalid number of nodes to recover!");
                exit(1);
            }
        }
        printf("\n");

        double whole_elapsed_time = 0, max_elapsed_time = 0, min_elapsed_time = 1000;
        double whole_inner_time1 = 0, whole_inner_time2 = 0;

 /*       for (int k = 0; k < testsNum / 1000; ++k)
        {
            whole_elapsed_time = 0;
            for (int i = 0; i < 1000; ++i)
            {
                whole_elapsed_time += elapsed_time[k*1000 + i];
            }
            printf("%g ", ((double)(nodesToRecoverNum) * 1024 * 16 * 1000) / (whole_elapsed_time * 1e6));
        }
        printf("\n");
        whole_elapsed_time = 0; */

        for (int i = 0; i < testsNum; ++i)
        {
            whole_elapsed_time += elapsed_time[i];
            whole_inner_time1 += inner_time1[i];
            whole_inner_time2 += inner_time2[i];
            if (elapsed_time[i] < min_elapsed_time)
                min_elapsed_time = elapsed_time[i];
            if (elapsed_time[i] > max_elapsed_time)
                max_elapsed_time = elapsed_time[i];
        }
//        printf("Av. elapsed time: %f usec\n", (whole_elapsed_time * 1e6) / testsNum);
//        printf("  min. elapsed time: %f usec\n", min_elapsed_time * 1e6);
//        printf("  max. elapsed time: %f usec\n", max_elapsed_time * 1e6);
//        printf("Inner time 1: %f usec\n", (whole_inner_time1 * 1e6) / testsNum);
//        printf("Inner time 2: %f usec\n", (whole_inner_time2 * 1e6) / testsNum);
        printf("Av. speed: %g Mb/s\n", ((double)(nodesToRecoverNum) * 1024 * 16 * testsNum) / (whole_elapsed_time * 1e6));
//        printf("  max. speed: %g Mb/s\n", ((double)(nodesToRecoverNum) * 1024 * 16) / (min_elapsed_time * 1e6));
        //printf("  min. speed: %g Mb/s\n", ((double)(nodesToRecoverNum) * 1024 * 16) / (max_elapsed_time * 1e6));
    }

    printf("max log 0x%X", maxLog);
    return 0;
}