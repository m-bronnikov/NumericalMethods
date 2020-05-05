#ifndef R_MATRIX
#define R_MATRIX

#include "RVector.h"
#include "Complex.h"

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
        double koeff = scalar_multy(&z_one, &q[i]);
        for(unsigned j = 0; j < vector_size(z); ++j){
            z->_data[j] -= koeff * q[i]._data[j];
        }
    }
    delete_RVector(&z_one);
}

void lancosh_method(RMatrix* matrix, RVector* a, RVector* b, unsigned k){
    if(!k){
        k = matrix_size(matrix);
    }

    RVector* q = malloc(sizeof(RVector) * k);
    RVector z = new_vector(matrix_size(matrix));

    resize_vector(a, k);
    resize_vector(b, k);

    RVector_random(&z);

    double beta = vector_norm(&z);
    unsigned i = 0;

    for(; i < k; ++i){
        if(beta < EPSILON){
            vector_resize(b, i);
            vector_resize(a, i);
            break;
        }
        q[i] = new_vector(matrix_size(matrix));
        vector_dev_num(&z, beta);
        vector_copy(&q[i], &z);

        matrix_x_vector(matrix, &q[i], &z);
        a->_data[i] = scalar_multiply(&q[i], &z);

        gramma_shmidt(&z, q, i);
        // gramma_shmidt(&z, q, i);

        beta = b->_data[i] = vector_norm(&z);
    }

    delete_vector(&z);
    for(unsigned j = 0; j < i; ++j){
        delete_vector(&q[j]);
    }
    free(q);
}

#endif