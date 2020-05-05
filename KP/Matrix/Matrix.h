#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <omp.h>

typedef struct QMatrix{
    double** _matrix;
    unsigned _size;
}QMatrix;

typedef struct RVector{
    double* _data;
    unsigned _size;
}RVector;


void new_RVector(RVector* v, unsigned size);
void new_matrix(QMatrix* matrix, unsigned size);

void delete_matrix(QMatrix* matrix);
void delete_RVector(RVector* v);

unsigned matrix_size(QMatrix* matrix);
unsigned RVector_size(RVector* vec);

double* matrix_elem(QMatrix* matrix, unsigned i, unsigned j);
double* RVector_elem(RVector* vec, unsigned i);

void set_random(RVector* vec);

double RVector_norm(RVector* vec);

void nomilize_RVector(RVector* vec);

void matrix_x_RVector(QMatrix* matrix, RVector* vec, RVector* ans);
double scalar_multy(RVector* v1, RVector* v2);
void RVector_dev_num(RVector* vec, double num);
void RVector_copy(RVector* v1, RVector* v2);
void RVector_minus(RVector* v1, RVector* v2);

void RVector_resize(RVector* vec, unsigned size);

void read_matrix(QMatrix* matrix, const char* path);

void print_matrix(QMatrix* matrix);

void gramma_shmidt(RVector* z, RVector* q, unsigned m);

void get_values(QMatrix* matrix, RVector* values, unsigned k);


void lancosh_method(QMatrix* matrix, RVector* a, RVector* b, unsigned k);

void QR_values(RVector* a, RVector* b, RVector* values);




#endif