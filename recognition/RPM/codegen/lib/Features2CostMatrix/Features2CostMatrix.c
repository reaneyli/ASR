/*
 * File: Features2CostMatrix.c
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 29-Nov-2021 17:59:00
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "Features2CostMatrix.h"
#include "Features2CostMatrix_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : emxArray_real_T *F1
 *                emxArray_real_T *F2
 *                const emxArray_real_T *Points1
 *                const emxArray_real_T *Points2
 *                int maxdist
 *                emxArray_real_T *C
 * Return Type  : void
 */
void Features2CostMatrix(emxArray_real_T *F1, emxArray_real_T *F2, const
  emxArray_real_T *Points1, const emxArray_real_T *Points2, int maxdist,
  emxArray_real_T *C)
{
  emxArray_real_T *b_F1;
  int i0;
  int loop_ub;
  int b_loop_ub;
  int i1;
  unsigned int bu;
  unsigned int unnamed_idx_1;
  double P;
  double a;
  int exitg1;
  long long i2;
  emxInit_real_T(&b_F1, 2);
  i0 = b_F1->size[0] * b_F1->size[1];
  b_F1->size[0] = F1->size[1];
  b_F1->size[1] = F1->size[0];
  emxEnsureCapacity_real_T(b_F1, i0);
  loop_ub = F1->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_loop_ub = F1->size[1];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_F1->data[i1 + b_F1->size[0] * i0] = F1->data[i0 + F1->size[0] * i1];
    }
  }

  i0 = F1->size[0] * F1->size[1];
  F1->size[0] = b_F1->size[0];
  F1->size[1] = b_F1->size[1];
  emxEnsureCapacity_real_T(F1, i0);
  loop_ub = b_F1->size[0] * b_F1->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    F1->data[i0] = b_F1->data[i0];
  }

  /* 559*3375 */
  i0 = b_F1->size[0] * b_F1->size[1];
  b_F1->size[0] = F2->size[1];
  b_F1->size[1] = F2->size[0];
  emxEnsureCapacity_real_T(b_F1, i0);
  loop_ub = F2->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_loop_ub = F2->size[1];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_F1->data[i1 + b_F1->size[0] * i0] = F2->data[i0 + F2->size[0] * i1];
    }
  }

  i0 = F2->size[0] * F2->size[1];
  F2->size[0] = b_F1->size[0];
  F2->size[1] = b_F1->size[1];
  emxEnsureCapacity_real_T(F2, i0);
  loop_ub = b_F1->size[0] * b_F1->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    F2->data[i0] = b_F1->data[i0];
  }

  emxFree_real_T(&b_F1);

  /* 558*3375 */
  bu = (unsigned int)F1->size[0];
  unnamed_idx_1 = (unsigned int)F2->size[0];
  i0 = C->size[0] * C->size[1];
  C->size[0] = (int)bu;
  C->size[1] = (int)unnamed_idx_1;
  emxEnsureCapacity_real_T(C, i0);
  loop_ub = (int)bu * (int)unnamed_idx_1;
  for (i0 = 0; i0 < loop_ub; i0++) {
    C->data[i0] = 0.0;
  }

  /* 559*558 */
  i0 = F1->size[1];
  for (b_loop_ub = 0; b_loop_ub < i0; b_loop_ub++) {
    P = (F1->data[F1->size[0] * b_loop_ub] + 2.2204460492503131E-16) + F2->
      data[F2->size[0] * b_loop_ub];

    /* bsxfun(@plus,F1(:,i)+eps,F2(:,i)'); */
    a = F1->data[F1->size[0] * b_loop_ub] - F2->data[F2->size[0] * b_loop_ub];

    /*      M=bsxfun(@minus,F1(:,i),F2(:,i)'); */
    i1 = C->size[0] * C->size[1];
    loop_ub = C->size[0] * C->size[1];
    emxEnsureCapacity_real_T(C, loop_ub);
    P = a * a / P;
    loop_ub = i1 - 1;
    for (i1 = 0; i1 <= loop_ub; i1++) {
      C->data[i1] += P;
    }
  }

  /*  C=zeros([size(F1,1) size(F2,1) 3375]);%559*558 */
  /*  for i=1:size(F1,1) */
  /*      P=F1(i,:)+eps+F2; */
  /*      M=F1(i,:)-F2; */
  /*      value=((M.^2)./P); */
  /*  %     tic */
  /*  %     sum_value=sum(((M.^2)./P),2); */
  /*      C(i,:,:)=value; */
  /*  end */
  /*  C=sum(C,3); */
  a = Points1->data[0] - Points2->data[0];
  P = Points1->data[Points1->size[0]] - Points2->data[Points2->size[0]];
  b_loop_ub = maxdist;
  loop_ub = 1;
  bu = 2U;
  do {
    exitg1 = 0;
    if ((bu & 1U) != 0U) {
      i2 = (long long)b_loop_ub * loop_ub;
      if (i2 > 2147483647LL) {
        i2 = 2147483647LL;
      } else {
        if (i2 < -2147483648LL) {
          i2 = -2147483648LL;
        }
      }

      loop_ub = (int)i2;
    }

    bu >>= 1U;
    if ((int)bu == 0) {
      exitg1 = 1;
    } else {
      i2 = (long long)b_loop_ub * b_loop_ub;
      if (i2 > 2147483647LL) {
        i2 = 2147483647LL;
      } else {
        if (i2 < -2147483648LL) {
          i2 = -2147483648LL;
        }
      }

      b_loop_ub = (int)i2;
    }
  } while (exitg1 == 0);

  b_loop_ub = 0;
  if (a * a + P * P > loop_ub) {
    b_loop_ub = 1;
  }

  loop_ub = b_loop_ub - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    C->data[0] = rtInf;
  }
}

/*
 * File trailer for Features2CostMatrix.c
 *
 * [EOF]
 */
