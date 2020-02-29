#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

void ave(){
    cout << "============================================================" << endl;
    cout << "|                    LABORATORY WORK №2                    |" << endl;
    cout << "|                    NUMERICAL METHODS                     |" << endl;
    cout << "|                         Task №1                          |" << endl;
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

void print_statement(double alfa, double x0){
    cout << "============================================================" << endl;
    cout << "|                      EXERCSICE:                          |" << endl;
    cout << "============================================================" << endl;
    cout << "| We need to solve equation by simple itteraion and Newton |" << endl;
    cout << "| method.                                                  |" << endl;
    cout << "============================================================" << endl;
    cout << "Equation:" << endl;
    cout << "x^3 + x^2 - x + 0.5 = 0" << endl;
    cout << "Accuracy: " << alfa << endl;
    cout << "x0 = " << x0 << endl;
    cout << "=============================================================" << endl;
}


void print_solution(double x_n, int itter_n, double x_i, int itter_i){
    cout << "============================================================" << endl;
    cout << "|                      ANSWER:                             |" << endl;
    cout << "============================================================" << endl;
    cout << "|                   NEWTON METHOD:                         |" << endl;
    cout << "============================================================" << endl;
    cout << "x = " << x_n << endl;
    cout << "Itterations: " << itter_n << endl;
    cout << "============================================================" << endl;
    cout << "|               SIMPLE ITTERATIONS METHOD:                 |" << endl;
    cout << "============================================================" << endl;
    cout << "x = " << x_i << endl;
    cout << "Itterations: " << itter_i << endl;
    cout << "============================================================" << endl;
}



double f(double x){
    return x*x*x + x*x - x - 0.5;
}

double df_dx(double x){
    return 3.0*x*x + 2.0*x - 1.0;
}

double d2f_dx2(double x){
    return 6.0*x + 2.0;
}

bool newton_condition(double x0){
    return f(x0) * d2f_dx2(x0) > 0;
}

double fitta(double x){
    return pow(0.5 + x - x*x, 1.0 / 3.0);
}

double dfitta_dx(double x){
    return (1.0 - 2.0*x) / (3.0 * pow(0.5 + x - x*x, 2.0 / 3.0));
}

bool itteration_condition(double x0){
    return abs(dfitta_dx(x0)) < 1;
    //return true;
}

double newton_method(double x0, double alfa, int& itter){
    double x_j, x_k = x0;
    itter = 0;
    do{
        x_j = x_k;
        x_k -= f(x_j) / df_dx(x_j);
        ++itter;
    }while(abs(x_k - x_j) >= alfa);
    return x_k;
}

double itteration_method(double x0, double alfa, int& itter){
    double x_j, x_k = x0;
    itter = 0;
    do{
        x_j = x_k;
        x_k = fitta(x_j);
        ++itter;
    }while(abs(x_k - x_j) >= alfa);
    return x_k;
}

int main(){
    double x0, alfa, x_n, x_i;
    int itter_n, itter_i;

    ave();
    cin >> alfa;
    cin >> x0;
    print_statement(alfa, x0);

    if(!newton_condition(x0)){
        cout << "Wrong x0 for Newton method. Please try to choice another x0." << endl;
        bye();
        return 0;
    }
    if(!itteration_condition(x0)){
        cout << "Wrong x0 for simmple itteration method. Please try to choice another x0." << endl;
        bye();
        return 0;
    }

    x_n = newton_method(x0, alfa, itter_n);
    x_i = itteration_method(x0, alfa, itter_i);

    print_solution(x_n, itter_n, x_i, itter_i);
    bye();
    return 0;
}