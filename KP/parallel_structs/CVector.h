#ifndef C_VECTOR
#define C_VECTOR

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Complex.h"

#include <omp.h>


typedef struct CVector{
    Complex* _data;
    unsigned _size;
}CVector;

CVector create_cvector(){
    CVector ans;
    ans._data = NULL;
    ans._size = 0;
    return ans;
}

CVector new_cvector(unsigned size){
    CVector ans;
    ans._data = (Complex*) malloc(sizeof(Complex) * size);
    if(ans._data){
        ans._size = size;
    }
    return ans;
}

void delete_cvector(CVector* vec){
    free(vec->_data);
    vec->_data = NULL;
    vec->_size = 0;
}

CVector* resize_cvector(CVector* vec, unsigned size){
    vec->_data = (Complex*)realloc(vec->_data, sizeof(Complex) * size);
    if(vec->_data){
        vec->_size = size;
    }else{
        vec->_size = 0;
    }
    return vec;
}

unsigned cvector_size(CVector* vec){
    return vec->_size;
}

CVector* fill_zeros(CVector* vec){
    for(unsigned i = 0; i < cvector_size(vec); ++i){
        vec->_data[i].re = vec->_data[i].im = 0.0;
    }
    return vec;
}


Complex* i_cvec(CVector* vec, unsigned i){
    return &vec->_data[i];
}


CVector* cvector_copy(CVector* v1, CVector* v2){
    memcpy(v1->_data, v2->_data, cvector_size(v1)*sizeof(Complex));
    return v1;
}

void print_cvector(CVector* vec){
    printf("Vector of size %u:", cvector_size(vec));
    for(unsigned i = 0; i < cvector_size(vec); ++i){
        printf("|");
        print_complex(vec->_data[i]);
    }
    printf("|\n");
}


#endif

