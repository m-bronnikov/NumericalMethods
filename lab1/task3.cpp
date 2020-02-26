// Made by Max Bronnikov

#include <vector>
#include <iostream>
#include "Matrix/matrix.h"
#include <cmath>

using namespace std;

void ave(){
    cout << "============================================================" << endl;
    cout << "|                    LABORATORY WORK №1                    |" << endl;
    cout << "|                    NUMERICAL METHODS                     |" << endl;
    cout << "|                         Task №3                          |" << endl;
    cout << "|                       Variant №4                         |" << endl;
    cout << "|                                                          |" << endl;
    cout << "|                                Student: Bronnikov M.A.   |" << endl;
    cout << "|                                   Date: 23.02.2020       |" << endl;
    cout << "|                                                          |" << endl;
    cout << "|                                                          |" << endl;
    cout << "|                       Moscow, 2020                       |" << endl;
    cout << "============================================================" << endl;
}

void print_statement(const Matrix& A, const vector<double>& b, double alfa){
    cout << "============================================================" << endl;
    cout << "|                      EXERCSICE:                          |" << endl;
    cout << "============================================================" << endl;
    cout << "|  We need to solve SLAU by simple ittearation and Zeidel  |" << endl;
    cout << "|  methods. Analyse count of iiteraions with accuracy.     |" << endl;
    cout << "============================================================" << endl;
    cout << "SLAU matrix:" << endl;
    cout << A << endl;
    cout << "Free vector:" << endl;
    for(unsigned i = 0; i < b.size(); ++i){
        cout.width(8);
        cout << b[i] << endl;
    }
    cout << endl;
    cout << "Accuracy: " << alfa << endl;
    cout << "=============================================================" << endl;
}

void bye(){
    cout << "============================================================" << endl;
    cout << "|                         EXIT                             |" << endl;
    cout << "============================================================" << endl;
}

/*
void matrix_init(Matrix& A){
    A = Matrix(3, 3);
    A[0][0] = 10.0;
    A[0][1] = 1.0;
    A[0][2] = 1.0;

    A[1][0] = 2.0;
    A[1][1] = 10.0;
    A[1][2] = 1.0;

    A[2][0] = 2.0;
    A[2][1] = 2.0;
    A[2][2] = 10.0;
}

void vector_init(vector<double>& b){
    b.resize(3);
    b[0] = 12.0;
    b[1] = 13.0;
    b[2] = 14.0;
}
*/

void matrix_init(Matrix& A){
    A = Matrix(4, 4);
    A[0][0] = 26.0;
    A[0][1] = -9.0;
    A[0][2] = -8.0;
    A[0][3] = 8.0;

    A[1][0] = 9.0;
    A[1][1] = -21.0;
    A[1][2] = -2.0;
    A[1][3] = 8.0;

    A[2][0] = -3.0;
    A[2][1] = 2.0;
    A[2][2] = -18.0;
    A[2][3] = 8.0;

    A[3][0] = 1.0;
    A[3][1] = -6.0;
    A[3][2] = -1.0;
    A[3][3] = 11.0;
}

void vector_init(vector<double>& b){
    b.resize(4);
    b[0] = 20.0;
    b[1] = -164.0;
    b[2] = 140.0;
    b[3] = -81.0;
}



void print_vector_x(const vector<double>& x){
    for(unsigned i = 0; i < x.size(); ++i){
        cout << 'x' << i + 1 << " = " << x[i] << " ";
    }
    cout << endl;
}

void print_solution(const vector<double>& x, int itter){
    cout << "============================================================" << endl;
    cout << "|                      ANSWER:                             |" << endl;
    cout << "============================================================" << endl;
    cout << "x solution:" << endl;
    print_vector_x(x);
    cout << "============================================================" << endl;
    cout << "Itterations: " << itter << endl;
    cout << "============================================================" << endl;
}

vector<double> vector_minus(const vector<double>& a, const vector<double>& b){
    vector<double> minus = a;
    for(unsigned i = 0; i < minus.size(); ++i){
        minus[i] -= b[i];
    }
    return minus;
}

vector<double> vector_plus(const vector<double>& a, const vector<double>& b){
    vector<double> plus = a;
    for(unsigned i = 0; i < plus.size(); ++i){
        plus[i] += b[i];
    }
    return plus;
}

double norm_of_vector(const vector<double>& vec){
    double norm = 0.0;
    for(unsigned i = 0; i < vec.size(); ++i){
        norm += vec[i] * vec[i];
    }
    return sqrt(norm);
}

int simple_itteration(const Matrix& A, const vector<double>& b, vector<double>& x, double alfa){
    Matrix M = A;
    x.resize(b.size());
    vector<double> last(b.size(), 0.0), r = b;
    //double coeff = 0.0;

    // commpute new matrix alfa and vector beta
    if(!M.is_quadratic()){
        throw "Wrong matrix! Try again!";
    }
    for(int i = 0; i < M.get_n(); ++i){
        if(!A[i][i]){
            throw "Wrong matrix! Try again!";
        }
        for(int j = 0; j < M.get_m(); ++j){
            M[i][j] = i == j ? 0.0 : -A[i][j] / A[i][i];
        }
        r[i] /= A[i][i];
    }

    x = r;
    /*
    coeff = M.get_norm();
    coeff /= 1 - coeff;
    */

    /* MAKE ITTERATIONS HERE: */

    int itter = 0;

    for(itter = 0; norm_of_vector(vector_minus(x, last)) > alfa; ++itter){
        x.swap(last);
        x = vector_plus(r, M * last);
    }

    return itter;
}

int zeidels_method(const Matrix& A, const vector<double>& b, vector<double>& x, double alfa){
    Matrix M = A;
    x.resize(b.size());
    vector<double> last(b.size(), 0.0), r = b;
    // double coeff = 0.0;

    // commpute new matrix alfa and vector beta
    if(!M.is_quadratic()){
        throw "Wrong matrix! Try again!";
    }
    for(int i = 0; i < M.get_n(); ++i){
        if(!A[i][i]){
            throw "Wrong matrix! Try again!";
        }
        for(int j = 0; j < M.get_m(); ++j){
            M[i][j] = i == j ? 0.0 : -A[i][j] / A[i][i];
        }
        r[i] /= A[i][i];
    }

    x = r;
    /*
    coeff = M.get_norm();
    coeff /= 1 - coeff;
    */

    /* MAKE ITTERATIONS HERE: */

    int itter = 0;

    for(itter = 0; norm_of_vector(vector_minus(x, last)) > alfa; ++itter){
        x.swap(last);
        x = r;
        for(int i = 0; i < M.get_n(); ++i){
            for(int j = 0; j < i; ++j){
                x[i] += x[j] * M[i][j];
            }
            for(int j = i; j < M.get_m(); ++j){
                x[i] += last[j] * M[i][j];
            }
        }
    }

    return itter;
}


int main(){
    ave();

    Matrix A;
    vector<double> x, b;
    vector<vector<double>> vec;
    double accuracy = 0.001;

    matrix_init(A);
    vector_init(b);
    print_statement(A, b, accuracy);

    //itterations:
    cout << "Simple Itterations Method:" << endl;
    int itter = simple_itteration(A, b, x, accuracy);
    print_solution(x, itter);

    //zeidel:
    cout << "Zeidels Method of Itterations:" << endl;
    itter = zeidels_method(A, b, x, accuracy);
    print_solution(x, itter);
    bye();

    return 0;
}