//
// File: Features2CostMatrix_or.cpp
//
// MATLAB Coder version            : 4.1
// C/C++ source code generated on  : 29-Nov-2021 18:46:19
//

// Include Files
#include "rt_nonfinite.h"
#include "Features2CostMatrix_or.h"
#include "Features2CostMatrix_or_emxutil.h"

// Function Definitions

//
// Arguments    : emxArray_real_T *F1
//                emxArray_real_T *F2
//                const emxArray_real_T *Points1
//                const emxArray_real_T *Points2
//                double maxdist
//                emxArray_real_T *C
// Return Type  : void
//
void Features2CostMatrix_or(emxArray_real_T *F1, emxArray_real_T *F2, const
  emxArray_real_T *Points1, const emxArray_real_T *Points2, double maxdist,
  emxArray_real_T *C)
{
  emxArray_real_T *b_F1;
  int i0;
  int loop_ub;
  int b_loop_ub;
  int i1;
  unsigned int unnamed_idx_0;
  unsigned int unnamed_idx_1;
  int i;
  double P;
  double a;
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

  // 559*3375
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

  // 558*3375
  unnamed_idx_0 = (unsigned int)F1->size[0];
  unnamed_idx_1 = (unsigned int)F2->size[0];
  i0 = C->size[0] * C->size[1];
  C->size[0] = (int)unnamed_idx_0;
  C->size[1] = (int)unnamed_idx_1;
  emxEnsureCapacity_real_T(C, i0);
  loop_ub = (int)unnamed_idx_0 * (int)unnamed_idx_1;
  for (i0 = 0; i0 < loop_ub; i0++) {
    C->data[i0] = 0.0;
  }

  // 559*558
  i0 = F1->size[1];
  for (i = 0; i < i0; i++) {
    P = (F1->data[F1->size[0] * i] + 2.2204460492503131E-16) + F2->data[F2->
      size[0] * i];

    // bsxfun(@plus,F1(:,i)+eps,F2(:,i)');
    a = F1->data[F1->size[0] * i] - F2->data[F2->size[0] * i];

    //      M=bsxfun(@minus,F1(:,i),F2(:,i)');
    i1 = C->size[0] * C->size[1];
    b_loop_ub = C->size[0] * C->size[1];
    emxEnsureCapacity_real_T(C, b_loop_ub);
    P = a * a / P;
    loop_ub = i1 - 1;
    for (i1 = 0; i1 <= loop_ub; i1++) {
      C->data[i1] += P;
    }
  }

  //  C=zeros([size(F1,1) size(F2,1) 3375]);%559*558
  //  for i=1:size(F1,1)
  //      P=F1(i,:)+eps+F2;
  //      M=F1(i,:)-F2;
  //      value=((M.^2)./P);
  //  %     tic
  //  %     sum_value=sum(((M.^2)./P),2);
  //      C(i,:,:)=value;
  //  end
  //  C=sum(C,3);
  a = Points1->data[0] - Points2->data[0];
  P = Points1->data[Points1->size[0]] - Points2->data[Points2->size[0]];
  b_loop_ub = 0;
  if (a * a + P * P > maxdist * maxdist) {
    b_loop_ub = 1;
  }

  loop_ub = b_loop_ub - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    C->data[0] = rtInf;
  }
}

//
// File trailer for Features2CostMatrix_or.cpp
//
// [EOF]
//
