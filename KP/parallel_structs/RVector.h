#ifndef R_VECTOR
#define R_VECTOR

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <omp.h>


typedef struct RVector{
    double* _data;
    unsigned _size;
}RVector;

RVector create_vector(){
    RVector ans;
    ans._data = NULL;
    ans._size = 0;
    return ans;
}

RVector new_vector(unsigned size){
    RVector ans;
    ans._data = (double*) malloc(sizeof(double) * size);
    if(ans._data){
        ans._size = size;
    }
    return ans;
}

void delete_vector(RVector* vec){
    free(vec->_data);
    vec->_data = NULL;
    vec->_size = 0;
}

RVector* resize_vector(RVector* vec, unsigned size){
    vec->_data = (double*)realloc(vec->_data, sizeof(double) * size);
    if(vec->_data){
        vec->_size = size;
    }else{
        vec->_size = 0;
    }
    return vec;
}


unsigned vector_size(RVector* vec){
    return vec->_size;
}

double* i_vec(RVector* vec, unsigned i){
    return &vec->_data[i];
}

double parallel_scalar_multiply(RVector* v1, RVector* v2, unsigned procnum){
    double ans = 0.0;
    omp_set_num_threads(procnum);
    #pragma omp parallel for reduction(+:ans)
    for(unsigned i = 0; i < vector_size(v1); ++i){
        ans += v1->_data[i] * v2->_data[i];
    }
    return ans;
}

/*
RVector* vector_dev_num(RVector* vec, double num){
    for(unsigned i = 0; i < vector_size(vec); ++i){
        vec->_data[i] /= num;
    }
    return vec;
}

RVector* vector_mult_num(RVector* vec, double num){
    for(unsigned i = 0; i < vector_size(vec); ++i){
        vec->_data[i] *= num;
    }
    return vec;
}

*/
RVector* vector_copy(RVector* v1, RVector* v2){
    memcpy(v1->_data, v2->_data, vector_size(v1)*sizeof(double));
    return v1;
}

/*
RVector* vector_minus(RVector* v1, RVector* v2){
    for(unsigned i = 0; i < vector_size(v1); ++i){
        v1->_data[i] -= v2->_data[i];
    }
    return v1;
}



RVector* vector_plus(RVector* v1, RVector* v2){
    for(unsigned i = 0; i < vector_size(v1); ++i){
        v1->_data[i] += v2->_data[i];
    }
    return v1;
}
*/
RVector* parallel_fill_random(RVector* vec, unsigned procnum){
    srand(42);
    omp_set_num_threads(procnum);
    #pragma omp parallel for
    for(unsigned i = 0; i < vector_size(vec); ++i){
        vec->_data[i] = (double) rand();
    }
}


double parallel_vector_norm(RVector* vec, unsigned procnum){
    return sqrt(parallel_scalar_multiply(vec, vec, procnum));
}

/*
RVector* nomilize_vector(RVector* vec){
    double norm = vector_norm(vec);
    for(unsigned i = 0; i < vector_size(vec); ++i){
        vec->_data[i] /= norm;
    }
    return vec;
}
*/

void print_vector(RVector* vec){
    printf("Vector of size %u:", vector_size(vec));
    for(unsigned i = 0; i < vector_size(vec); ++i){
        printf("|%lg", vec->_data[i]);
    }
    printf("|\n");
}


#endif

