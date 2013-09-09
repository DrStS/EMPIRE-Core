#ifndef HELPERFUNCTIONS_H_
#define HELPERFUNCTIONS_H_

int *doubleArrayToIntArray(double *arrayDouble, int size) {
    int *arrayInt = new int[size];
    for (int i=0; i<size; i++)
        arrayInt[i] = (int) arrayDouble[i];
    return arrayInt;
}



#endif /* HELPERFUNCTIONS_H_ */
