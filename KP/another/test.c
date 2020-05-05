#include "RVector.h"
#include "RMatrix.h"

int main(){
    RMatrix matrix = read_matrix("../data/m1.txt");
    print_matrix(&matrix);
    for(unsigned i = 0; i < matrix_size(&matrix); ++i){
        printf("line %u: %u\n", i, line_length(&matrix._lines[i]));
    }
    delete_matrix(&matrix);
    printf("Program ends!");
    return 0;
}