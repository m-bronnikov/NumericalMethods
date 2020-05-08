#ifndef R_MATRIX
#define R_MATRIX

#include "RVector.h"
#include "CVector.h"
#include <stdbool.h>

typedef struct RLine{
    RVector _vector;
    unsigned* _indexes;
}RLine;

typedef struct RMatrix{
    RLine* _lines;
    unsigned _size;
}RMatrix;

unsigned line_length(RLine* line){
    return line->_vector._size;
}

RLine* resize_line(RLine* line, unsigned size){
    if(resize_vector(&line->_vector, size)->_size){
        line->_indexes = (unsigned*)realloc(line->_indexes, size * sizeof(unsigned));
    }
    return line;
}

RLine create_line(){
    RLine line;
    line._indexes = NULL;
    line._vector = create_vector();
    return line;
}

void delete_line(RLine* line){
    delete_vector(&line->_vector);
    free(line->_indexes);
    line->_indexes = NULL;
}

RMatrix create_matrix(){
    RMatrix matrix;
    matrix._lines = NULL;
    matrix._size = 0;
    return matrix;
}

void push_line(RLine* line, unsigned i, unsigned idx, double elem){
    line->_vector._data[i] = elem;
    line->_indexes[i] = idx;
}

unsigned matrix_idx(RMatrix* matrix, unsigned i, unsigned j){
    return matrix->_lines[i]._indexes[j];
}

double matrix_elem(RMatrix* matrix, unsigned i, unsigned j){
    return  *i_vec(&matrix->_lines[i]._vector, j);
}

unsigned matrix_size(RMatrix* matrix){
    return matrix->_size;
}

RVector* matrix_x_vector(RMatrix* matrix, RVector* vec, RVector* ans){
    for(unsigned i = 0; i < matrix_size(matrix); ++i){
        *i_vec(ans, i) = 0.0;
        for(unsigned j = 0; j < line_length(&matrix->_lines[i]); ++j){
            *i_vec(ans, i) += matrix_elem(matrix, i, j) * (*i_vec(vec, matrix_idx(matrix, i, j)));
        }
    }
    return ans;
}


RMatrix new_matrix(unsigned size){
    RMatrix matrix;
    matrix._lines = (RLine*)malloc(sizeof(RLine)*size);
    matrix._size = size;
    for(unsigned i = 0; i < size; ++i){
        matrix._lines[i] = create_line();
    }
    return matrix;
}

RMatrix* delete_matrix(RMatrix* matrix){
    for(unsigned i = 0; i < matrix_size(matrix); ++i){
        delete_line(&matrix->_lines[i]);
    }
    free(matrix->_lines);
    matrix->_lines = NULL;
    matrix->_size = 0;
    return matrix;
}

RMatrix read_matrix(const char* path){
    FILE* f = fopen(path, "r");
    unsigned size;
    fscanf(f, "%u", &size);
    RMatrix matrix = new_matrix(size);
    for(unsigned i = 0; i < size; ++i){
        unsigned cur_size = 0;
        for(unsigned j = 0; j < size; ++j){
            double num;
            fscanf(f, "%lg", &num);
            if(!num){
                continue;
            }
            if(line_length(&matrix._lines[i]) == cur_size){
                if(cur_size){
                    resize_line(&matrix._lines[i], (cur_size) << 1);
                }else{
                    resize_line(&matrix._lines[i], 2);
                }
            }
            push_line(&matrix._lines[i], cur_size++, j, num);
        }
        resize_line(&matrix._lines[i], cur_size);
    }
    fclose(f);
    return matrix;
}


RMatrix read_florida(const char* path){
    FILE* f = fopen(path, "r");
    unsigned size, elems;
    fscanf(f, "%u %u", &size, &elems);
    RMatrix matrix = new_matrix(size);
    unsigned* curr_sizes = (unsigned*) malloc(size*sizeof(unsigned));
    memset(curr_sizes, 0, size*sizeof(unsigned));
    for(unsigned k = 0; k < elems; ++k){
        unsigned i, j;
        double num;
        fscanf(f, "%u %u %lg", &i, &j, &num);
        --i; --j;
        // push i j
        if(line_length(&matrix._lines[i]) == curr_sizes[i]){
            if(curr_sizes[i]){
                resize_line(&matrix._lines[i], (curr_sizes[i]) << 1);
            }else{
                resize_line(&matrix._lines[i], 2);
            }
        }
        push_line(&matrix._lines[i], curr_sizes[i]++, j, num);
        if(i == j){
            continue;
        }
        // push j i
        if(line_length(&matrix._lines[j]) == curr_sizes[j]){
            if(curr_sizes[i]){
                resize_line(&matrix._lines[j], (curr_sizes[j]) << 1);
            }else{
                resize_line(&matrix._lines[j], 2);
            }
        }
        push_line(&matrix._lines[j], curr_sizes[j]++, i, num);
    }
    
    for(unsigned i = 0; i < matrix_size(&matrix); ++i){
        resize_line(&matrix._lines[i], curr_sizes[i]);
    }
    
    free(curr_sizes);
    fclose(f);
    return matrix;
}

void print_matrix(RMatrix* matrix){
    printf("Matrix %ux%u:\n", matrix_size(matrix), matrix_size(matrix));
    for(unsigned i = 0; i < matrix_size(matrix); ++i){
        unsigned idx = 0;
        for(unsigned j = 0; j < line_length(&matrix->_lines[i]); ++j){
            while(idx++ < matrix_idx(matrix, i, j)){
                printf("0.0\t");
            }
            printf("%lg\t", matrix_elem(matrix, i, j));
        }
        while(idx++ < matrix_size(matrix)){
                printf("0.0\t");
        }
        printf("\n");
    }
}

void gramma_shmidt(RVector* z, RVector* q, unsigned m){
    RVector z_one;
    z_one = new_vector(vector_size(z));
    vector_copy(&z_one, z);
    for(unsigned i = 0; i <= m; ++i){
        double koeff = scalar_multiply(&z_one, &q[i]);
        for(unsigned j = 0; j < vector_size(z); ++j){
            z->_data[j] -= koeff * q[i]._data[j];
        }
    }
    delete_vector(&z_one);
}

void lancosh_z(RVector* z, RVector* q1, RVector* q0, double a, double b){
    for(unsigned i = 0; i < vector_size(z); ++i){
        *i_vec(z, i) -= (*i_vec(q1, i)) * a + (*i_vec(q0, i)) * b;
    }
}

void lancosh_method(RMatrix* matrix, RVector* a, RVector* b, unsigned k){
    if(!k){
        k = matrix_size(matrix);
    }

    RVector* q = (RVector*)malloc(sizeof(RVector) * k);
    RVector z = new_vector(matrix_size(matrix));

    resize_vector(a, k);
    resize_vector(b, k);

    fill_random(&z);

    double beta = vector_norm(&z);
    unsigned i = 0;

    for(; i < k; ++i){
        if(beta < EPSILON){
            resize_vector(b, i);
            resize_vector(a, i);
            break;
        }

        q[i] = new_vector(matrix_size(matrix));
        vector_dev_num(&z, beta);
        vector_copy(&q[i], &z);

        matrix_x_vector(matrix, &q[i], &z);
        *i_vec(a, i) = scalar_multiply(&q[i], &z);

        gramma_shmidt(&z, q, i);
        gramma_shmidt(&z, q, i);

        /*
        if(i){
            lancosh_z(&z, &q[i], &q[i-1], *i_vec(a, i), beta);
        }else{
            lancosh_z(&z, &q[0], &q[0], *i_vec(a, i), 0);
        }
        */

        beta = b->_data[i] = vector_norm(&z);
    }

    delete_vector(&z);
    for(unsigned j = 0; j < i; ++j){
        delete_vector(&q[j]);
    }
    free(q);
}

void def_QR_step(double* a, double* b, unsigned len){
    if(len < 2){
        return;
    }
    for(unsigned i = 0; i < len - 1; ++i){
        if(fabs(b[i]) < EPSILON){
            def_QR_step(a, b, i + 1);
            def_QR_step(&a[i+1], &b[i+1], len - i - 1);
            return;
        }
    }

    double ak = a[0];
    double bk = b[0];
    double x = 0.0;

    double shift = a[len - 1];
    double norm = sqrt((ak - shift)*(ak - shift) + bk*bk);
    double gamma = (ak - shift) / norm;
    double sigma = bk / norm;

    a[0] = gamma*(gamma*ak + sigma*bk) + sigma*(gamma*bk + sigma*a[1]);

    b[0] = gamma*(gamma*bk + sigma*a[1]) - sigma*(sigma*bk + gamma*ak);

    x = b[1] * sigma;

    a[1] = gamma*(gamma*a[1] - sigma*bk) - sigma*(gamma*bk - sigma*ak);

    b[1] *= gamma;

    // k step:

    for(unsigned k = 1; k < len - 1; ++k){
        bk = b[k-1];
        norm = sqrt(bk*bk + x*x);
        b[k-1] = norm;

        gamma = bk / norm;
        sigma = x / norm;

        ak = a[k];
        bk = b[k];

        a[k] = gamma*(gamma*ak + sigma*bk) + sigma*(gamma*bk + sigma*a[k+1]);

        b[k] = gamma*(gamma*bk + sigma*a[k+1]) - sigma*(sigma*bk + gamma*ak);

        a[k+1] = gamma*(gamma*a[k+1] - sigma*bk) - sigma*(gamma*bk - sigma*ak);

        x = b[k+1] * sigma;
        
        b[k+1] *= gamma;
    }
}

void QR_step(RVector* a, RVector *b){
    def_QR_step(a->_data, b->_data, vector_size(a));
}


double complex_check(Complex last_1, Complex last_2, Complex cur_1, Complex cur_2){
    Complex r1, r2;
    // x[j]
    r1.re = cur_1.re - last_1.re; // x
    r1.im = cur_1.im - last_1.im; // y
    // x[j+1]
    r2.re = cur_2.re - last_2.re; // x
    r2.im = cur_2.im - last_2.im; // y
    double s1 = modulo_complex(r1);
    double s2 = modulo_complex(r2);
    return s1 > s2 ? s1 : s2;
}

void QR_values(RVector* a, RVector* b, CVector* values, double eps){
    resize_cvector(values, vector_size(a));
    bool flag = true;
    fill_zeros(values);
    Complex curr_1, curr_2;

    do{
        flag = false;
        QR_step(a, b);
        *i_vec(b, vector_size(b) - 1) = 0.0;
        for(int i = 0; i < vector_size(a); ++i){
            // real check:
            if(fabs(*i_vec(b, i)) < eps){
                i_cvec(values, i)->re = *i_vec(a, i);
                i_cvec(values, i)->im = 0.0;
            }else{
                quadrat_solver_first(1.0, -(*i_vec(a, i)) - (*i_vec(a, i+1)), 
                    (*i_vec(a, i))*(*i_vec(a, i+1)) - (*i_vec(b, i))*(*i_vec(b, i)), 
                    &curr_2, &curr_1);
                if(complex_check(curr_1, curr_2, *i_cvec(values, i), *i_cvec(values, i+1)) > eps){
                    flag = true;
                }
                *i_cvec(values, i) = curr_1;
                *i_cvec(values, i+1) = curr_2; 
                ++i;
            }
            // complex check:
        }
    }while(flag);
}

#endif