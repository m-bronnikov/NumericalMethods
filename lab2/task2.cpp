#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "Matrix/matrix.h"

using namespace std;

// params:
const double a = 1.0; // 4 variant
const int n = 2; // 2 variables in equation
const double start_delta = 0.2;
const double min_delta = 0.01;
const double search_step = 0.01;

void ave(){
    cout << "============================================================" << endl;
    cout << "|                    LABORATORY WORK №2                    |" << endl;
    cout << "|                    NUMERICAL METHODS                     |" << endl;
    cout << "|                         Task №2                          |" << endl;
    cout << "|                       Variant №4                         |" << endl;
    cout << "|                                                          |" << endl;
    cout << "|                                Student: Bronnikov M.A.   |" << endl;
    cout << "|                                   Date: 28.02.2020       |" << endl;
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

void print_statement(double alfa, const vector<double>& x0){
    cout << "============================================================" << endl;
    cout << "|                      EXERCSICE:                          |" << endl;
    cout << "============================================================" << endl;
    cout << "| We need to solve equation system by simple itteraion and |" << endl;
    cout << "| Newton methods.                                          |" << endl;
    cout << "============================================================" << endl;
    cout << "Equation system:" << endl;
    cout << endl;
    cout << "x1 - cos(x2) = 1" << endl;
    cout << "x2 - lg(x1 + 1) = 1" << endl;
    cout << endl;
    cout << "Accuracy: " << alfa << endl;
    cout << "x0 = (" << x0[0];
    for(int i = 1; i < n; ++i){
        cout << ", " << x0[i];
    }
    cout << ")" << endl;
    cout << "=============================================================" << endl;
}


void print_solution(const vector<double>& x_n, int itter_n, const vector<double>& x_i, int itter_i, const double q){
    cout << "============================================================" << endl;
    cout << "|                      ANSWER:                             |" << endl;
    cout << "============================================================" << endl;
    cout << "|                   NEWTON METHOD:                         |" << endl;
    cout << "============================================================" << endl;
    cout << "x = (" << x_n[0];
    for(int i = 1; i < n; ++i){
        cout << ", " << x_n[i];
    }
    cout << ")" << endl;
    cout << "Itterations: " << itter_n << endl;
    cout << "============================================================" << endl;
    cout << "|               SIMPLE ITTERATIONS METHOD:                 |" << endl;
    cout << "============================================================" << endl;
    if(q >= 1.0){
        cout << "Sufficient condition not done!" << endl;
    }else{
        cout << "Sufficient condition done with q: " << q << endl;
    }
    cout << "x = (" << x_i[0];
    for(int i = 1; i < n; ++i){
        cout << ", " << x_i[i];
    }
    cout << ")" << endl;
    cout << "Itterations: " << itter_i << endl;
    cout << "============================================================" << endl;
}

// plus and minus:
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

// norm for methods stop
double norm_of_vector(const vector<double>& vec){
    double norm = 0.0;
    for(unsigned i = 0; i < vec.size(); ++i){
        norm += vec[i] * vec[i];
    }
    return sqrt(norm);
}

// norm for zeidel stop
double norm_of_vectors(const vector<double>& x1, const vector<double>& x2){
    double norm = 0.0;
    if(x1.size() != x2.size()){
        throw "Wrong sizes of vectors";
    }
    for(unsigned i = 0; i < x1.size(); ++i){
        norm = max(abs(x1[i] - x2[i]), norm);
    }
    return norm;
}

// For SLAU solve:
void zeidels_method(const Matrix& A, const vector<double>& b, vector<double>& x, double alfa){
    Matrix M = A;
    x.resize(b.size());
    vector<double> last(b.size(), 0.0), r = b;
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
    for(int itter = 0; norm_of_vector(vector_minus(x, last)) > alfa; ++itter){
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
}

/* This methods set my equation. redefine for your task */

double f1(const vector<double>& x){
    return x[0] - cos(x[1]) - 1.0;
}

double df1_dx1(const vector<double>& x){
    return 1.0;
}

double df1_dx2(const vector<double>& x){
    return sin(x[1]);
}


double f2(const vector<double>& x){
    return x[1] - log10(x[0] + 1.0) - a;
}


double df2_dx2(const vector<double>& x){
    return 1.0;
}

double df2_dx1(const vector<double>& x){
    return -1.0 / ((x[0] + 1) * log(10));
}

double fitta1(const vector<double>& x){
    return cos(x[1]) + 1.0;
}

double dfitta1_dx1(const vector<double>& x){
    return 0.0;
}

double dfitta1_dx2(const vector<double>& x){
    return -sin(x[1]);
}

double fitta2(const vector<double>& x){
    return log10(x[0] + 1.0) + a;
}

double dfitta2_dx1(const vector<double>& x){
    return 1.0 / ((x[0] + 1) * log(10));
}

double dfitta2_dx2(const vector<double>& x){
    return 0.0;
}


void set_matrix(Matrix& A, const vector<double>& x){
    A = Matrix(n, n);
    A[0][0] = df1_dx1(x);
    A[0][1] = df1_dx2(x);
    A[1][0] = df2_dx1(x);
    A[1][1] = df2_dx2(x);
}

void set_dfitta_matrix(Matrix& A, const vector<double>& x){
    A = Matrix(n, n);
    A[0][0] = dfitta1_dx1(x);
    A[0][1] = dfitta1_dx2(x);
    A[1][0] = dfitta2_dx1(x);
    A[1][1] = dfitta2_dx2(x);
}

void set_vector(vector<double>& fx, const vector<double>& x){
    fx.resize(n);
    fx[0] = -f1(x);
    fx[1] = -f2(x);
}

/* End of methods for redefine */


// Main methods:
vector<double> nex_step(const vector<double>& x){
    vector<double> ans(n);
    ans[0] = fitta1(x);
    ans[1] = fitta2(x);
    return ans;
}

int newton_method(const vector<double>& x0, vector<double>& x, double alfa){
    int itter = 0;
    vector<double> x_i;
    vector<double> fx, dx;
    Matrix J;
    x = x0;
    do{
        x_i.swap(x);
        set_matrix(J, x_i);
        set_vector(fx, x_i);
        zeidels_method(J, fx, dx, alfa);
        x = vector_plus(x_i, dx);
        ++itter;
    }while(norm_of_vectors(x, x_i) > alfa);
    return itter;
}

double find_q(const vector<double>& x0){
    double delta = start_delta * 2.0;
    vector<double> ans_x = x0;
    Matrix Dfitta;
    double ans;
    do{
        delta /= 2.0;
        // search x with max norm
        for(unsigned i = 0; i < x0.size(); ++i){
            double maximum = 0.0;
            vector<double> x = x0;
            for(double v = x0[i] - delta; v <= x0[i] + delta; v += search_step){
                x[i] = v;
                set_dfitta_matrix(Dfitta, x);
                double norm = Dfitta.get_norm();
                if(norm > maximum){
                    ans_x[i] = v;
                    maximum = norm;
                }
            }
        }
        set_dfitta_matrix(Dfitta, ans_x); 
        ans = Dfitta.get_norm();
    }while(ans >= 1.0 && delta >= min_delta);
    return ans;
}

int itteration_method(const vector<double>& x0, vector<double>& x, double alfa, double& q){
    int itter = 0;
    vector<double> x_i;
    vector<double> fx, dx;
    q = find_q(x0);
    if(q < 1.0){
        alfa *= (1.0 - q);
        alfa /= q;
    }
    x = x0;
    do{
        x_i.swap(x);
        x = nex_step(x_i);
        ++itter; 
    }while(norm_of_vectors(x, x_i) > alfa);
    return itter;
}

int main(){
    ave();
    vector<double> x0(n), x_n(n), x_i(n, 0);
    double alfa, q;
    int itter_n, itter_i = 0;
    // read accuracy and x0:
    cin >> alfa;
    for(int i = 0; i < n; ++i){
        cin >> x0[i];
    }
    print_statement(alfa, x0);

    itter_n = newton_method(x0, x_n, alfa);
    itter_i = itteration_method(x0, x_i, alfa, q);

    print_solution(x_n, itter_n, x_i, itter_i, q);
    bye();
    return 0;
}