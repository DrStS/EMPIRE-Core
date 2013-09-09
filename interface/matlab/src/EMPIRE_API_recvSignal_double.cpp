#include "matrix.h"
#include "mex.h"
#include <assert.h>
#include "EMPIRE_API.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  assert(nrhs==3);
  assert(nlhs==0);

  // signal name
  assert(mxIsChar(prhs[0]));
  char name[EMPIRE_API_NAME_STRING_LENGTH];
  mxGetNChars(prhs[0], name, EMPIRE_API_NAME_STRING_LENGTH);

  // size of array
  assert(mxIsDouble(prhs[1]));
  assert(mxGetNumberOfElements(prhs[1]) == 1);
  int sizeOfArray = (int)mxGetPr(prhs[1])[0]; // cast from double to int

  // signal array
  assert(mxIsDouble(prhs[2]));
  assert(mxGetNumberOfElements(prhs[2]) == sizeOfArray);
  double *signal = mxGetPr(prhs[2]);

  EMPIRE_API_recvSignal_double(name, sizeOfArray, signal);
}
