#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  assert(nrhs==3);
  assert(nlhs==0);

#define NAME_IN   prhs[0]
#define SIZE_IN   prhs[1]
#define SIGNAL_IN prhs[2]

  // signal name
  assert(mxIsChar(NAME_IN));
  char name[EMPIRE_API_NAME_STRING_LENGTH];
  mxGetNChars(NAME_IN, name, EMPIRE_API_NAME_STRING_LENGTH);

  // size of array
  assert(mxIsDouble(SIZE_IN));
  assert(mxGetNumberOfElements(SIZE_IN) == 1);
  int sizeOfArray = (int)mxGetPr(SIZE_IN)[0]; // cast from double to int

  // signal array
  assert(mxIsDouble(SIGNAL_IN));
  assert(mxGetNumberOfElements(SIGNAL_IN) == sizeOfArray);
  double *signal = mxGetPr(SIGNAL_IN);

  EMPIRE_API_recvSignal_double(name, sizeOfArray, signal);

#undef NAME_IN
#undef SIZE_IN
#undef SIGNAL_IN
}
