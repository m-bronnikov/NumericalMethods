#include "Matrix.h"

const double epsilon = 0.00000000000001;

unsigned matrix_size(QMatrix* matrix){
    return matrix->_size;
}

unsigned RVector_size(RVector* vec){
    return vec->_size;
}

double* matrix_elem(QMatrix* matrix, unsigned i, unsigned j){
    return &matrix->_matrix[i][j];
}

double* RVector_elem(RVector* vec, unsigned i){
    return &vec->_data[i];
}

void RVector_random(RVector* vec){
    srand(42);
    for(unsigned i = 0; i < RVector_size(vec); ++i){
        vec->_data[i] = (double) rand();
    }
}

double RVector_norm(RVector* vec){
    double norm = 0.0;
    for(unsigned i = 0; i < RVector_size(vec); ++i){
        norm += vec->_data[i] * vec->_data[i];
    }
    return sqrt(norm);
}

void nomilize_RVector(RVector* vec){
    double norm = RVector_norm(vec);
    for(unsigned i = 0; i < RVector_size(vec); ++i){
        vec->_data[i] /= norm;
    }
}

void new_RVector(RVector* v, unsigned size){
    v->_size = size;
    v->_data = (double *)malloc(size * sizeof(double));
}

void new_matrix(QMatrix* matrix, unsigned size){
    matrix->_size = size;
    matrix->_matrix = malloc(size * sizeof(double*));
    for(unsigned i = 0; i < size; ++i){
        matrix->_matrix[i] = malloc(size * sizeof(double));
    }
}

void read_matrix(QMatrix* matrix, const char* path){
    FILE* f = fopen(path, "r");
    unsigned size;
    fscanf(f, "%u", &size);
    new_matrix(matrix, size);
    for(unsigned i = 0; i < size; ++i){
        for(unsigned j = 0; j < size; ++j){
            fscanf(f, "%lg", &matrix->_matrix[i][j]);
        }
    }
    fclose(f);
}

void delete_matrix(QMatrix* matrix){
    for(unsigned i = 0; i < matrix_size(matrix); ++i){
        free(matrix->_matrix[i]);
    }
    free(matrix->_matrix);
    matrix->_size = 0;
}

void print_matrix(QMatrix* matrix){
    printf("Matrix:\n");
    for(unsigned i = 0; i < matrix_size(matrix); ++i){
        for(unsigned  j = 0; j < matrix_size(matrix); ++j){
            printf("%lg ", matrix->_matrix[i][j]);
        }
        printf("\n");
    }
}

void delete_RVector(RVector* v){
    v->_size = 0;
    free(v->_data);
}

void matrix_x_RVector(QMatrix* matrix, RVector* vec, RVector* ans){
    for(unsigned i = 0; i < RVector_size(vec); ++i){
        ans->_data[i] = 0.0;
        for(unsigned j = 0; j < matrix_size(matrix); ++j){
            ans->_data[i] += matrix->_matrix[i][j] * vec->_data[j];
        }   
    }
}

double scalar_multy(RVector* v1, RVector* v2){
    double ans = 0.0;
    for(unsigned i = 0; i < RVector_size(v1); ++i){
        ans += v1->_data[i] * v2->_data[i];
    }
    return ans;
}

void RVector_dev_num(RVector* vec, double num){
    for(unsigned i = 0; i < RVector_size(vec); ++i){
        vec->_data[i] /= num;
    }
}

void RVector_copy(RVector* v1, RVector* v2){
    memcpy(v1->_data, v2->_data, RVector_size(v1)*sizeof(double));
}

void RVector_minus(RVector* v1, RVector* v2){
    for(unsigned i = 0; i < RVector_size(v1); ++i){
        v1->_data[i] -= v2->_data[i];
    }
}

void RVector_resize(RVector* vec, unsigned size){
    vec->_size = size;
    vec->_data = realloc(vec->_data, size*sizeof(double));
}

void gramma_shmidt(RVector* z, RVector* q, unsigned m){
    RVector z_one;
    new_RVector(&z_one, RVector_size(z));
    RVector_copy(&z_one, z);
    for(unsigned i = 0; i <= m; ++i){
        double koeff = scalar_multy(&z_one, &q[i]);
        for(unsigned j = 0; j < RVector_size(z); ++j){
            z->_data[j] -= koeff * q[i]._data[j];
        }
    }
    delete_RVector(&z_one);
}


void get_values(QMatrix* matrix, RVector* values, unsigned k){
    RVector a;
    RVector b;
    lancosh_method(matrix, &a, &b, k);
    QR_values(&a, &b, values);
    delete_RVector(&a);
    delete_RVector(&b);
}

void lancosh_method(QMatrix* matrix, RVector* a, RVector* b, unsigned k){
    if(!k){
        k = matrix_size(matrix);
    }

    RVector* q = malloc(sizeof(RVector) * k);
    RVector z;

    new_RVector(&z, matrix_size(matrix));
    new_RVector(a, k);
    new_RVector(b, k);

    RVector_random(&z);

    double beta = RVector_norm(&z);
    unsigned i = 0;

    for(; i < k; ++i){
        if(beta < epsilon){
            RVector_resize(b, i);
            RVector_resize(a, i);
            break;
        }
        new_RVector(&q[i], matrix_size(matrix));
        RVector_dev_num(&z, beta);
        RVector_copy(&q[i], &z);

        matrix_x_RVector(matrix, &q[i], &z);
        a->_data[i] = scalar_multy(&q[i], &z);

        gramma_shmidt(&z, q, i);
        gramma_shmidt(&z, q, i);

        beta = b->_data[i] = RVector_norm(&z);
    }

    delete_RVector(&z);
    for(unsigned j = 0; j < i; ++j){
        delete_RVector(&q[j]);
    }
    free(q);
}

void QRstep(RVector* a, RVector* b){
    return ;
}

void QR_values(RVector* a, RVector* b, RVector* values){
    
}