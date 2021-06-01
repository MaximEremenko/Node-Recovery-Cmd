#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <errno.h>
#include <vector>
#include <list>
#include <chrono>
#include <omp.h>

#include "Node.h"
#include "FieldDefs.h"
#include "FieldElement.h"
#include "ConstructionDefs.h"

using namespace std;

//=======================================================================================================================================
//=============================================================== Auxiliary functions ===================================================
//=======================================================================================================================================

// Inversion of Vandermond matrix of size 4x4
void vand4_inv(FieldElement B[4][4], FieldElement a, FieldElement b, FieldElement c, FieldElement d)
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
}

// Inversion of Vandermond matrix of size 3x3
void vand3_inv(FieldElement B[3][3], FieldElement a, FieldElement b, FieldElement c)
{
    FieldElement ab = a + b;
    FieldElement ac = a + c;
    FieldElement bc = b + c;

    FieldElement w0 = ~(ab * ac);
    FieldElement w1 = ~(ab * bc);
    FieldElement w2 = ~(ac * bc);

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
void vand2_inv(FieldElement B[2][2], FieldElement a, FieldElement b)
{
    FieldElement w0 = ~(a + b);

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

// 4 nodes recovery procedure
double recover_4_nodes(unsigned int* pNodesToRecoverIdx, Node* pNodes)
{
    int lambdasIdx[5] = {0, 0, 0, 0, 0};
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
        nodesToRecLambdas[i][0] = FieldElement(lambdas[AIdx][currLambdaIdx]) * sigmas[currNodeToRec];
        currNodeToRec = pNodesToRecoverIdx[1];
        AIdx = AIdxArray[currNodeToRec];
        currLambdaIdx = lambdasIdx[AIdx];
        nodesToRecLambdas[i][1] = FieldElement(lambdas[AIdx][currLambdaIdx]) * sigmas[currNodeToRec];
        currNodeToRec = pNodesToRecoverIdx[2];
        AIdx = AIdxArray[currNodeToRec];
        currLambdaIdx = lambdasIdx[AIdx];
        nodesToRecLambdas[i][2] = FieldElement(lambdas[AIdx][currLambdaIdx]) * sigmas[currNodeToRec];
        currNodeToRec = pNodesToRecoverIdx[3];
        AIdx = AIdxArray[currNodeToRec];
        currLambdaIdx = lambdasIdx[AIdx];
        nodesToRecLambdas[i][3] = FieldElement(lambdas[AIdx][currLambdaIdx]) * sigmas[currNodeToRec];
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

    auto start_time = chrono::high_resolution_clock::now();

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
                currSigmaSum[j] += ASigmas[j][k]*AData[i][j][k];
                currSigmaPow2Sum[j] += ASigmasPow2[j][k] * AData[i][j][k];
                currSigmaPow3Sum[j] += ASigmasPow3[j][k] * AData[i][j][k];
            }
        }
        sumCol[0] = currSum[0] + currSum[1] + currSum[2] + currSum[3] + currSum[4];
        sumCol[1] = currSigmaSum[0] * ACol[i][0][0] + currSigmaSum[1] * ACol[i][1][0] + currSigmaSum[2] * ACol[i][2][0] + currSigmaSum[3] * ACol[i][3][0] + currSigmaSum[4] * ACol[i][4][0];
        sumCol[2] = currSigmaPow2Sum[0] * ACol[i][0][1] + currSigmaPow2Sum[1] * ACol[i][1][1] + currSigmaPow2Sum[2] * ACol[i][2][1] + currSigmaPow2Sum[3] * ACol[i][3][1] + currSigmaPow2Sum[4] * ACol[i][4][1];
        sumCol[3] = currSigmaPow3Sum[0] * ACol[i][0][2] + currSigmaPow3Sum[1] * ACol[i][1][2] + currSigmaPow3Sum[2] * ACol[i][2][2] + currSigmaPow3Sum[3] * ACol[i][3][2] + currSigmaPow3Sum[4] * ACol[i][4][2];

        vand4_inv(invMatr, nodesToRecLambdas[i][0], nodesToRecLambdas[i][1],nodesToRecLambdas[i][2], nodesToRecLambdas[i][3]);
        
        currRecData[0] = invMatr[0][0] * sumCol[0] + invMatr[0][1] * sumCol[1] + invMatr[0][2] * sumCol[2] + invMatr[0][3] * sumCol[3];
        currRecData[1] = invMatr[1][0] * sumCol[0] + invMatr[1][1] * sumCol[1] + invMatr[1][2] * sumCol[2] + invMatr[1][3] * sumCol[3];
        currRecData[2] = invMatr[2][0] * sumCol[0] + invMatr[2][1] * sumCol[1] + invMatr[2][2] * sumCol[2] + invMatr[2][3] * sumCol[3];
        currRecData[3] = invMatr[3][0] * sumCol[0] + invMatr[3][1] * sumCol[1] + invMatr[3][2] * sumCol[2] + invMatr[3][3] * sumCol[3];

        pNodes[pNodesToRecoverIdx[0]].setData(i, currRecData[0]);
        pNodes[pNodesToRecoverIdx[1]].setData(i, currRecData[1]);
        pNodes[pNodesToRecoverIdx[2]].setData(i, currRecData[2]);
        pNodes[pNodesToRecoverIdx[3]].setData(i, currRecData[3]);

    }

    auto end_time = chrono::high_resolution_clock::now();
    auto elapsed_time = chrono::duration_cast<chrono::microseconds>(end_time - start_time);

    return (double)elapsed_time.count();
}

// 3 nodes recovery procedure
double recover_3_nodes(unsigned int* pNodesToRecoverIdx, Node* pNodes)
{
    int lambdasIdx[5] = { 0, 0, 0, 0, 0 };
    unsigned int i = 0, j = 0, k = 0;
    int AIdx, currNodeToRec = 0, currLambdaIdx = 0;

    FieldElement currSum[5], currSigmaSum[5], currSigmaPow2Sum[5], sumCol[3];
    FieldElement invMatr[3][3], currRecData[3];

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
        nodesToRecLambdas[i][0] = FieldElement(lambdas[AIdx][currLambdaIdx]) * sigmas[currNodeToRec];
        currNodeToRec = pNodesToRecoverIdx[1];
        AIdx = AIdxArray[currNodeToRec];
        currLambdaIdx = lambdasIdx[AIdx];
        nodesToRecLambdas[i][1] = FieldElement(lambdas[AIdx][currLambdaIdx]) * sigmas[currNodeToRec];
        currNodeToRec = pNodesToRecoverIdx[2];
        AIdx = AIdxArray[currNodeToRec];
        currLambdaIdx = lambdasIdx[AIdx];
        nodesToRecLambdas[i][2] = FieldElement(lambdas[AIdx][currLambdaIdx]) * sigmas[currNodeToRec];
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
                currSigmaSum[j] += ASigmas[j][k] * AData[i][j][k];
                currSigmaPow2Sum[j] += ASigmasPow2[j][k] * AData[i][j][k];
            }
        }

        sumCol[0] = currSum[0] + currSum[1] + currSum[2] + currSum[3] + currSum[4];
        sumCol[1] = currSigmaSum[0] * ACol[i][0][0] + currSigmaSum[1] * ACol[i][1][0] + currSigmaSum[2] * ACol[i][2][0] + currSigmaSum[3] * ACol[i][3][0] + currSigmaSum[4] * ACol[i][4][0];
        sumCol[2] = currSigmaPow2Sum[0] * ACol[i][0][1] + currSigmaPow2Sum[1] * ACol[i][1][1] + currSigmaPow2Sum[2] * ACol[i][2][1] + currSigmaPow2Sum[3] * ACol[i][3][1] + currSigmaPow2Sum[4] * ACol[i][4][1];

        vand3_inv(invMatr, nodesToRecLambdas[i][0], nodesToRecLambdas[i][1], nodesToRecLambdas[i][2]);

        currRecData[0] = invMatr[0][0] * sumCol[0] + invMatr[0][1] * sumCol[1] + invMatr[0][2] * sumCol[2];
        currRecData[1] = invMatr[1][0] * sumCol[0] + invMatr[1][1] * sumCol[1] + invMatr[1][2] * sumCol[2];
        currRecData[2] = invMatr[2][0] * sumCol[0] + invMatr[2][1] * sumCol[1] + invMatr[2][2] * sumCol[2];

        pNodes[pNodesToRecoverIdx[0]].setData(i, currRecData[0]);
        pNodes[pNodesToRecoverIdx[1]].setData(i, currRecData[1]);
        pNodes[pNodesToRecoverIdx[2]].setData(i, currRecData[2]);

    }

    auto end_time = chrono::high_resolution_clock::now();
    auto elapsed_time = chrono::duration_cast<chrono::microseconds>(end_time - start_time);

    return (double)elapsed_time.count();
}

// 2 nodes recovery procedure
double recover_2_nodes(unsigned int* pNodesToRecoverIdx, Node* pNodes)
{
    int lambdasIdx[5] = { 0, 0, 0, 0, 0 };
    unsigned int i = 0, j = 0, k = 0;
    int AIdx, currNodeToRec = 0, currLambdaIdx = 0;

    FieldElement currSum[5], currSigmaSum[5], sumCol[2];
    FieldElement invMatr[2][2], currRecData[2];

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
        nodesToRecLambdas[i][0] = FieldElement(lambdas[AIdx][currLambdaIdx]) * sigmas[currNodeToRec];
        currNodeToRec = pNodesToRecoverIdx[1];
        AIdx = AIdxArray[currNodeToRec];
        currLambdaIdx = lambdasIdx[AIdx];
        nodesToRecLambdas[i][1] = FieldElement(lambdas[AIdx][currLambdaIdx]) * sigmas[currNodeToRec];
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
                currSigmaSum[j] += ASigmas[j][k] * AData[i][j][k];
            }
        }

        sumCol[0] = currSum[0] + currSum[1] + currSum[2] + currSum[3] + currSum[4];
        sumCol[1] = currSigmaSum[0] * ACol[i][0][0] + currSigmaSum[1] * ACol[i][1][0] + currSigmaSum[2] * ACol[i][2][0] + currSigmaSum[3] * ACol[i][3][0] + currSigmaSum[4] * ACol[i][4][0];

        vand2_inv(invMatr, nodesToRecLambdas[i][0], nodesToRecLambdas[i][1]);

        currRecData[0] = invMatr[0][0] * sumCol[0] + invMatr[0][1] * sumCol[1];
        currRecData[1] = invMatr[1][0] * sumCol[0] + invMatr[1][1] * sumCol[1];

        pNodes[pNodesToRecoverIdx[0]].setData(i, currRecData[0]);
        pNodes[pNodesToRecoverIdx[1]].setData(i, currRecData[1]);

    }

    auto end_time = chrono::high_resolution_clock::now();
    auto elapsed_time = chrono::duration_cast<chrono::microseconds>(end_time - start_time);

    return (double)elapsed_time.count();
}

// 1 node recovery procedure
double recover_1_node(unsigned int* pNodesToRecoverIdx, Node* pNodes)
{
    unsigned int currPos[4], currLambdasIdx[5];
    unsigned int nodeToRecIdx = pNodesToRecoverIdx[0];
    unsigned int nodeToRecAIdx = AIdxArray[nodeToRecIdx];

    unsigned int subBlockIdx = 0, currNodeIdx = 0, currAIdx = 0, tempAIdx = 0, i = 0, j = 0;
    
    FieldElement invMatr[4][4];
    FieldElement nodeToRecSigma = FieldElement(sigmas[nodeToRecIdx]);

    vand4_inv(invMatr, nodeToRecSigma * lambdas[nodeToRecAIdx][0],
                       nodeToRecSigma * lambdas[nodeToRecAIdx][1],
                       nodeToRecSigma * lambdas[nodeToRecAIdx][2],
                       nodeToRecSigma * lambdas[nodeToRecAIdx][3]);

    FieldElement diffAPartSum[4] = { ZERO_ELEMENT, ZERO_ELEMENT, ZERO_ELEMENT, ZERO_ELEMENT };
    FieldElement sameAPartSum[4] = { ZERO_ELEMENT, ZERO_ELEMENT, ZERO_ELEMENT, ZERO_ELEMENT };
    FieldElement ASum[4] = { ZERO_ELEMENT, ZERO_ELEMENT, ZERO_ELEMENT, ZERO_ELEMENT };

    unsigned int diffANodeNum[4] = { 0, 0, 0, 0 };
    unsigned int sameANodeNum = 0;

    unsigned int posArray[256][4] = {0};

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
                diffAPartSum[1] += diffASigma[i][j] * diffAData[subBlockIdx][i][j];
                diffAPartSum[2] += diffASigmaPow2[i][j] * diffAData[subBlockIdx][i][j];
                diffAPartSum[3] += diffASigmaPow3[i][j] * diffAData[subBlockIdx][i][j];
            }

            ASum[0] += diffAPartSum[0];
            ASum[1] += diffAPartSum[1] * diffACol[subBlockIdx][i][0];
            ASum[2] += diffAPartSum[2] * diffACol[subBlockIdx][i][1];
            ASum[3] += diffAPartSum[3] * diffACol[subBlockIdx][i][2];
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
                sameAPartSum[1] += sameAData[subBlockIdx][j][i] * sameACol[j][0];
                sameAPartSum[2] += sameAData[subBlockIdx][j][i] * sameACol[j][1];
                sameAPartSum[3] += sameAData[subBlockIdx][j][i] * sameACol[j][2];;
            }
            ASum[0] += sameAPartSum[0];
            ASum[1] += sameASigma[i] * sameAPartSum[1];
            ASum[2] += sameASigmaPow2[i] * sameAPartSum[2];
            ASum[3] += sameASigmaPow3[i] * sameAPartSum[3];
        }

        for (i = 0; i < 4; ++i)
        {
            for (j = 0; j < 4; ++j)
            {
                recData[i] += invMatr[i][j] * ASum[j];
            }
        }

        pNodes[nodeToRecIdx].setData(posArray[subBlockIdx][0], recData[0]);
        pNodes[nodeToRecIdx].setData(posArray[subBlockIdx][1], recData[1]);
        pNodes[nodeToRecIdx].setData(posArray[subBlockIdx][2], recData[2]);
        pNodes[nodeToRecIdx].setData(posArray[subBlockIdx][3], recData[3]);
    }

    auto end_time = chrono::high_resolution_clock::now();
    auto elapsed_time = chrono::duration_cast<chrono::microseconds>(end_time - start_time);

    return (double)elapsed_time.count();
}

//=======================================================================================================================================
//================================================================== Main function ======================================================
//=======================================================================================================================================


double elapsed_time[10000] = { 0 };
FieldElement testData[1024];
Node pNodes[154];
Node pTestNodes[4];

int main()
{

    for (unsigned int i = 0; i < 154; ++i)
    {
        for (unsigned int j = 0; j < 1024; ++j)
        {
            testData[j] = FieldElement((fe_type)(rand() % (1 << 12)));
        }
        pNodes[i] = Node(FieldElement(sigmas[i]), testData);
    }

    unsigned int pNodesToRecoverIdx[4] = { 150, 151, 152, 153 };
    recover_4_nodes(pNodesToRecoverIdx, pNodes);
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

            S = S + FieldElement(lambdas[AIdx][currLambdaIdx]) * FieldElement(sigmas[i]) * pNodes[i].getData(j);
        }

        if (S.getElement() != 0)
            failFlag = 1;

        for (int i = 0; i < 154; ++i)
        {
            int AIdx = AIdxArray[i];
            int currLambdaIdx = lambdasIdx[AIdx];

            S = S + FieldElement(lambdas_pow2[AIdx][currLambdaIdx]) * FieldElement(sigmas_pow2[i]) * pNodes[i].getData(j);
        }

        if (S.getElement() != 0)
            failFlag = 1;

        for (int i = 0; i < 154; ++i)
        {
            int AIdx = AIdxArray[i];
            int currLambdaIdx = lambdasIdx[AIdx];

            S = S + FieldElement(lambdas_pow3[AIdx][currLambdaIdx]) * FieldElement(sigmas_pow3[i]) * pNodes[i].getData(j);
        }

        if (S.getElement() != 0)
            failFlag = 1;
    }
    if (failFlag)
        printf("FAILED!\n");
    else
        printf("PASS!\n");

    for (int i = 0; i < 4; ++i)
        pTestNodes[i] = pNodes[pNodesToRecoverIdx[i]];

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

    printf("3 nodes recovery check: ");
    recover_3_nodes(pNodesToRecoverIdx, pNodes);
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

    printf("2 nodes recovery check: ");
    recover_2_nodes(pNodesToRecoverIdx, pNodes);
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

    int testsNum = 10000;
    unsigned int pTestNodesToRecoverIdx[4] = { 1, 2, 3, 4 };
    unsigned int rndId = 0;

    printf("\n==Reconstruction speed test:\n");
    for (nodesToRecoverNum = 1; nodesToRecoverNum < 5; ++nodesToRecoverNum)
    {
        printf("---- %d nodes reconstruction:\n", nodesToRecoverNum);
        for (int tests = 0; tests < testsNum; ++tests)
        {
            printf("# tests performed: %d\r", tests + 1);

            for (unsigned int i = 0; i < 154; ++i)
            {
                for (unsigned int j = 0; j < 1024; ++j)
                {
                    testData[j] = FieldElement((fe_type)(rand() % (1 << 12)));
                }
                pNodes[i] = Node(FieldElement(sigmas[i]), testData);
            }
            recover_4_nodes(pNodesToRecoverIdx, pNodes);

            rndId = rand() % 154;
            for (unsigned int i = 0; i < 4; ++i)
                pTestNodesToRecoverIdx[i] = (rndId + i) % 154;

            switch (nodesToRecoverNum)
            {
            case 1:
                elapsed_time[tests] = recover_1_node(pTestNodesToRecoverIdx, pNodes);
                break;
            case 2:
                elapsed_time[tests] = recover_2_nodes(pTestNodesToRecoverIdx, pNodes);
                break;
            case 3:
                elapsed_time[tests] = recover_3_nodes(pTestNodesToRecoverIdx, pNodes);
                break;
            case 4:
                elapsed_time[tests] = recover_4_nodes(pTestNodesToRecoverIdx, pNodes);
                break;
            default:
                printf("Invalid number of nodes to recover!");
                exit(1);
            }
        }
        printf("\n");

        double whole_elapsed_time = 0, max_elapsed_time = 0, min_elapsed_time = 1000;

        for (int i = 0; i < testsNum; ++i)
        {
            whole_elapsed_time += elapsed_time[i];
            if (elapsed_time[i] < min_elapsed_time)
                min_elapsed_time = elapsed_time[i];
            if (elapsed_time[i] > max_elapsed_time)
                max_elapsed_time = elapsed_time[i];
        }

        printf("Av. elapsed time: %f usec\n", whole_elapsed_time / testsNum);
        printf("  min. elapsed time: %f usec\n", min_elapsed_time);
        printf("  max. elapsed time: %f usec\n", max_elapsed_time);
        printf("Av. speed: %g Mb/s\n", ((double)(nodesToRecoverNum) * 1024 * 12 * testsNum) / (whole_elapsed_time));
        printf("  max. speed: %g Mb/s\n", ((double)(nodesToRecoverNum) * 1024 * 12) / (min_elapsed_time));
        printf("  min. speed: %g Mb/s\n", ((double)(nodesToRecoverNum) * 1024 * 12) / (max_elapsed_time));
    }

    return 0;
}