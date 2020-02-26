// Made by Max Bronnikov

#include <vector>
#include <iostream>
#include "Matrix/matrix.h"

using namespace std;

void ave(){
    cout << "============================================================" << endl;
    cout << "|                    LABORATORY WORK №1                    |" << endl;
    cout << "|                    NUMERICAL METHODS                     |" << endl;
    cout << "|                         Task №2                          |" << endl;
    cout << "|                       Variant №4                         |" << endl;
    cout << "|                                                          |" << endl;
    cout << "|                                Student: Bronnikov M.A.   |" << endl;
    cout << "|                                   Date: 23.02.2020       |" << endl;
    cout << "|                                                          |" << endl;
    cout << "|                                                          |" << endl;
    cout << "|                       Moscow, 2020                       |" << endl;
    cout << "============================================================" << endl;
}

void print_statement(const Matrix& A, const vector<double>& b){
    cout << "============================================================" << endl;
    cout << "|                      EXERCSICE:                          |" << endl;
    cout << "============================================================" << endl;
    cout << "|  We need to solve SLAU by race method.                   |" << endl;
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

void bye(){
    cout << "============================================================" << endl;
    cout << "|                         EXIT                             |" << endl;
    cout << "============================================================" << endl;
}

/*
void matrix_init(Matrix& A){
    A = Matrix(4, 4);
    A[0][0] = 8.0;
    A[0][1] = -2.0;
    A[0][2] = 0.0;
    A[0][3] = 0.0;

    A[1][0] = -1.0;
    A[1][1] = 6.0;
    A[1][2] = -2.0;
    A[1][3] = 0.0;

    A[2][0] = 0.0;
    A[2][1] = 2.0;
    A[2][2] = 10.0;
    A[2][3] = -4.0;

    A[3][0] = 0.0;
    A[3][1] = 0.0;
    A[3][2] = -1.0;
    A[3][3] = 6.0;
}

void vector_init(vector<double>& b){
    b.resize(4);
    b[0] = 6.0;
    b[1] = 3.0;
    b[2] = 8.0;
    b[3] = 5.0;
}
*/

void matrix_init(Matrix& A){
    A = Matrix(5, 5);
    A[0][0] = -14.0;
    A[0][1] = -6.0;
    A[0][2] = 0.0;
    A[0][3] = 0.0;
    A[0][4] = 0.0;

    A[1][0] = -9.0;
    A[1][1] = 15.0;
    A[1][2] = -1.0;
    A[1][3] = 0.0;
    A[1][4] = 0.0;

    A[2][0] = 0.0;
    A[2][1] = 1.0;
    A[2][2] = -11.0;
    A[2][3] = 1.0;
    A[2][4] = 0.0;

    A[3][0] = 0.0;
    A[3][1] = 0.0;
    A[3][2] = -7.0;
    A[3][3] = 12.0;
    A[3][4] = 3.0;

    A[4][0] = 0.0;
    A[4][1] = 0.0;
    A[4][2] = 0.0;
    A[4][3] = 6.0;
    A[4][4] = -7.0;
}

void vector_init(vector<double>& b){
    b.resize(5);
    b[0] = -78.0;
    b[1] = -73.0;
    b[2] = -38.0;
    b[3] = 77.0;
    b[4] = 91.0;
}


void print_vector_x(const vector<double>& x){
    for(unsigned i = 0; i < x.size(); ++i){
        cout << 'x' << i + 1 << " = " << x[i] << " ";
    }
    cout << endl;
}

void matrix_to_vecs(const Matrix& matrix, vector<vector<double>>& vec){
    if(!matrix.is_quadratic()){
        throw "Wrong matrix";
    }
    vec.clear();
    vec.assign(matrix.get_n(), vector<double>(3, 0.0));
    
    for(int i = 0; i < matrix.get_n(); ++i){
        vec[i][0] = i - 1 < 0 ? 0.0 : matrix[i][i - 1];
        vec[i][1] = matrix[i][i];
        vec[i][2] = i + 1 < matrix.get_m() ? matrix[i][i + 1] : 0.0;
    }
}

void race_method(vector<vector<double>>& vec, const vector<double>& b, vector<double>& x){
    x.assign(b.size(), 0.0);
    vector<double> P(b.size()), Q(b.size());
    P[0] = -vec[0][2] / vec[0][1];
    Q[0] = b[0] / vec[0][1];
    cout << "P[" << 0 << "] = " << P[0] << " Q[" << 0 << "] = " << Q[0] << endl;
    for(int i = 1; i < (int)x.size(); ++i){
        double z = (vec[i][1] + vec[i][0] * P[i-1]);
        P[i] = -vec[i][2];
        P[i] /= z;
        Q[i] = (b[i] - vec[i][0] * Q[i - 1]);
        Q[i] /= z;
        cout << "P[" << i << "] = " << P[i] << " Q[" << i << "] = " << Q[i] << endl;
    }
    x.back() = Q.back();
    for(int i = x.size() - 2; i >= 0; --i){
        x[i] = P[i]*x[i+1] + Q[i];
    }
}

void print_solution(const vector<double>& x){
    cout << "============================================================" << endl;
    cout << "|                      ANSWER:                             |" << endl;
    cout << "============================================================" << endl;
    cout << "x solution:" << endl;
    print_vector_x(x);
    cout << "============================================================" << endl;
}


int main(){
    ave();

    Matrix A;
    vector<double> x, b;
    vector<vector<double>> vec;

    matrix_init(A);
    vector_init(b);
    print_statement(A, b);

    if(A.is_three_diagonal()){
        matrix_to_vecs(A, vec);
    }else{
        cout << "ERROR! MATRIX IS NOT THREE DIAGONAL! EXIT!" << endl;
        bye();
        return 0;
    }

    race_method(vec, b, x);

    print_solution(x);
    bye();

    return 0;
}