//
//  kernel.metal
//  NumKit
//
//  Created by Koichi Oshio on 2021/01/01.
//

#include <metal_stdlib>
using namespace metal;

kernel void add_arrays (
    device const float  *inA,
    device const float  *inB,
    device float        *result,
    uint                index [[thread_position_in_grid]])
{
    result[index] = inA[index] + inB[index];
}

kernel void global (device const int    *iparam,
                    device const float  *fparam,
                    device const float  *x,
                    device const float  *y,
                    device float        *result,
                    uint2               index [[thread_position_in_grid]])
{
    float   pd1, pd2;       // local param
    float   tc1, tc2;       // local param 
    float   mnPd1, mxPd1;         // fparam (0, 1)
    float   mnPd2, mxPd2;         // fparam (0, 1)
    float   mnTc1, mxTc1;         // fparam (0, 1)
    float   mnTc2, mxTc2;         // fparam (0, 1)
    int     dim, len;       // iparam (0, 1)
    int     mode;           // iparam (2)  
    int     i, i0, i1, i2, i3;
    float   mse, val, min_mse;

    mnPd1 = fparam[0];
    mxPd1 = fparam[1];
    mnPd2 = fparam[2];
    mxPd2 = fparam[3];
    mnTc1 = fparam[4];
    mxTc1 = fparam[5];
    mnTc2 = fparam[6];
    mxTc2 = fparam[7];

    dim = iparam[0];
    len = iparam[1];
    mode = iparam[2];

    i0 = index.y;
    i1 = index.x;

    min_mse = 1e10;
    for (i2 = 0; i2 < dim; i2++) {
        for (i3 = 0; i3 < dim; i3++) {
            switch (mode) {
            default :
            case 0 :    // x: pd2, y: pd1
                pd1 = (float)i0 / dim * (mxPd1 - mnPd1) + mnPd1;
                pd2 = (float)i1 / dim * (mxPd2 - mnPd2) + mnPd2;
                tc1 = (float)i2 / dim * (mxTc1 - mnTc1) + mnTc1;
                tc2 = (float)i3 / dim * (mxTc2 - mnTc2) + mnTc2;
                break;
            case 1 :    // x: tc2, y: tc1
                pd1 = (float)i2 / dim * (mxPd1 - mnPd1) + mnPd1;
                pd2 = (float)i3 / dim * (mxPd2 - mnPd2) + mnPd2;
                tc1 = (float)i0 / dim * (mxTc1 - mnTc1) + mnTc1;
                tc2 = (float)i1 / dim * (mxTc2 - mnTc2) + mnTc2;
                break;
            case 2 :    // x : tc2, y: pd2
                pd1 = (float)i2 / dim * (mxPd1 - mnPd1) + mnPd1;
                pd2 = (float)i0 / dim * (mxPd2 - mnPd2) + mnPd2;
                tc1 = (float)i3 / dim * (mxTc1 - mnTc1) + mnTc1;
                tc2 = (float)i1 / dim * (mxTc2 - mnTc2) + mnTc2;
                break;
            }

            // error calc (1 pixel)
            mse = 0;
            for (i = 0; i < len; i++) {
                val = pd1 * exp(-x[i] / tc1) + pd2 * exp(-x[i] / tc2);
                mse += (val - y[i]) * (val - y[i]);
            }
            mse = sqrt(mse/len);
            if (min_mse > mse) {
                min_mse = mse;
            }
        }
    }

    result[index.y * dim + index.x] = min_mse;
//    switch (mode) {
//    case 0:
//        result[index.y * dim + index.x] = index.y;
//        break;
//    case 1:
//        result[index.y * dim + index.x] = index.x;
//        break;
//    case 2:
//        result[index.y * dim + index.x] = index.y + index.x;
//        break;
//    }
}

kernel void dft (device const int       *iparam,
                    device const float  *re,
                    device const float  *im,
                    device const float  *cs,
                    device const float  *sn,
                    device float        *resRe,
                    device float        *resIm,
                    
                    uint2               index [[thread_position_in_grid]])
{
    int     i, j, ix, ij, n;
    int     direction;  // 1: forward (REC_FORWARD), -1: inverse (REC_INVERSE)
    float   sumr, sumi;

    n = iparam[0];
    sumr = sumi = 0;
    for (j = 0; j < n; j++) {
        ij = index.x * j % n;
        ix = index.y * n + j;
        sumr +=  re[ix] * cs[ij] + im[ix] * sn[ij];
        sumi += -re[ix] * sn[ij] + im[ix] * cs[ij];
    }
    ix = index.y * n + index.x;
    resRe[ix] = sumr;
    resIm[ix] = sumi;
}
