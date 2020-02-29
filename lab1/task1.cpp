// Made by Max Bronnikov

#include <vector>
#include <iostream>
#include "Matrix/matrix.h"

using namespace std;


void compute_solution(const Matrix& U, const Matrix& L, const vector<double>& b, vector<double>& x){
    x.resize(U.get_m());
    // Compute solution:
    // vector z:

    vector<double> z(L.get_m());
    for(int i = 0; i < L.get_m(); ++i){
        z[i] = b[i];
        for(int j = 0; j < i; ++j){
            z[i] -= (L[j][i] * z[j]);
        }
    }
    
    // vector x:
    for(int i = U.get_m() - 1; i >= 0; --i){
        x[i] = z[i];
        for(int j = i + 1; j < U.get_m(); ++j){
            x[i] -= (x[j] * U[i][j]);
        }
        x[i] /= U[i][i];
    }
}


// Compute LU matrix, SLAU solution for enterd matrix (if LU not exist - exception) and returns matrix determinant
double gauss_completely(const Matrix& matrix, Matrix& L, Matrix& U, Matrix& X, const vector<double>& b, vector<double>& x){
    if(!matrix.is_quadratic()){
        throw "LU not exist! Try simple gauss method or enter lines in another order!";
    }
    // init:
    Matrix B(matrix.get_n(), matrix.get_m());
    B.make_ones();
    U = matrix;
    L = Matrix(matrix.get_n(), matrix.get_m());
    X = Matrix(matrix.get_n(), matrix.get_m());
    vector<double> b_1 = b;
    double determinant = 1.0;

    // First Gauss steps (compute L and U):
    for(int j = 0; j < U.get_m(); ++j){
        if(U[j][j] == 0.0){
            throw "LU not exist! Try simple gauss method or enter lines in another order!";
        }
        L[j][j] = 1.0;
        for(int i = j + 1; i < U.get_n(); ++i){
            double l_ij = U[i][j] / U[j][j];
            L[i][j] = l_ij;
            // line[i] - line[j]
            for(int k = j; k < U.get_m(); ++k){
                U[i][k] -= U[j][k] * l_ij;
            }
            // B matrix change for Back matrix compute:
            for(int k = 0; k < B.get_m(); ++k){
                B[i][k] -= B[j][k] * l_ij;
            }
            // B[i] - B[j]
            b_1[i] -= b_1[j] * l_ij;
        }
    }

    // Compute solution:
    compute_solution(U, L, b_1, x);

    // Compute back matrix like compute sol:
    B.transpose();
    for(int i = 0; i < B.get_n(); ++i){
        compute_solution(U, L, B[i], X[i]);
    }
    X.transpose();

    // Compute determinant:
    for(int i = 0; i < U.get_m(); ++i){
        determinant *= U[i][i];
    }

    return determinant;
}

void ave(){
    cout << "============================================================" << endl;
    cout << "|                    LABORATORY WORK №1                    |" << endl;
    cout << "|                    NUMERICAL METHODS                     |" << endl;
    cout << "|                         Task №1                          |" << endl;
    cout << "|                       Variant №4                         |" << endl;
    cout << "|                                                          |" << endl;
    cout << "|                                Student: Bronnikov M.A.   |" << endl;
    cout << "|                                   Date: 23.02.2020       |" << endl;
    cout << "|                                                          |" << endl;
    cout << "|                                                          |" << endl;
    cout << "|                       Moscow, 2020                       |" << endl;
    cout << "============================================================" << endl;
}

void bye(){
    cout << "============================================================" << endl;
    cout << "|                         EXIT                             |" << endl;
    cout << "============================================================" << endl;
}

/*
void matrix_init(Matrix& A){
    A = Matrix(4, 4);
    A[0][0] = -1.0;
    A[0][1] = -7.0;
    A[0][2] = -3.0;
    A[0][3] = -2.0;
    A[1][0] = -8.0;
    A[1][1] = 1.0;
    A[1][2] = -9.0;
    A[1][3] = 0.0;
    A[2][0] = 8.0;
    A[2][1] = 2.0;
    A[2][2] = -5.0;
    A[2][3] = -3.0;
    A[3][0] = -5.0;
    A[3][1] = 3.0;
    A[3][2] = 5.0;
    A[3][3] = -9.0;
}
*/

int size_init(){
    int size;
    cin >> size;
    return size;
}

void matrix_init(Matrix& A, int size){
    A = Matrix(size, size);
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            cin >> A[i][j];
        }
    }
}

/*
void vector_init(vector<double>& b){
    b.resize(4);
    b[0] = -12.0;
    b[1] = -60.0;
    b[2] = -91.0;
    b[3] = -43.0;
}
*/

void vector_init(vector<double>& b, int size){
    b.resize(size);
    for(int i = 0; i < size; ++i){
        cin >> b[i];
    }
}

void print_vector_x(const vector<double>& x){
    for(unsigned i = 0; i < x.size(); ++i){
        cout << 'x' << i + 1 << " = " << x[i] << " ";
    }
    cout << endl;
}

void print_solution(const Matrix& U, const Matrix& L, const Matrix& B, const vector<double>& x, double determinant){
    cout << "============================================================" << endl;
    cout << "|                      ANSWER:                             |" << endl;
    cout << "============================================================" << endl;
    cout << "L matrix:" << endl;
    cout << L << endl;
    cout << "============================================================" << endl;
    cout << "U matrix:" << endl;
    cout << U << endl;
    cout << "============================================================" << endl;
    cout << "Back matrix:" << endl;
    cout << B << endl;
    cout << "============================================================" << endl;
    cout << "x solution:" << endl;
    print_vector_x(x);
    cout << "============================================================" << endl;
    cout << "Determinant: " << determinant << endl;
    cout << "=============================================================" << endl;
    return;
}

void print_statement(const Matrix& A, const vector<double>& b){
    cout << "============================================================" << endl;
    cout << "|                      EXERCSICE:                          |" << endl;
    cout << "============================================================" << endl;
    cout << "| We need to solve SLAU by LU separating matrix with Gauss |" << endl;
    cout << "| method. Find determinant, back matrix, L and U matrix.   |" << endl;
    cout << "============================================================" << endl;
    cout << "SLAU matrix:" << endl;
    cout << A << endl;
    cout << "Free vector:" << endl;
    for(unsigned i = 0; i < b.size(); ++i){
        cout.width(8);
        cout << b[i] << endl;
    }
    cout << "=============================================================" << endl;
}

int main(){
    ave();
    Matrix L, U, B, A;
    vector<double> b, x;
    int size = size_init();
    matrix_init(A, size);
    vector_init(b, size);
    print_statement(A, b);

    // Main function use:
    double determinant = gauss_completely(A, L, U, B, b, x); 

    print_solution(U, L, B, x, determinant);
    
    bye();
    return 0;
}