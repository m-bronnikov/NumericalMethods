#include "parallel_structs/RVector.h"
#include "parallel_structs/CVector.h"
#include "parallel_structs/RMatrix.h"

const double eps = 0.0001;
const unsigned procnum = 2;

int main(int argc, char *argv[]){
    RMatrix matrix = read_florida("data/bcsstk16.mtx");
    RVector a = create_vector(), b = create_vector();
    CVector c = create_cvector();

    //print_matrix(&matrix);

    printf("Lancosh started!\n");

    double end, start = omp_get_wtime();

    unsigned k = 0; 
    if(argc > 1){
        k = atoi(argv[1]);
    }
    parallel_lancosh_method(&matrix, &a, &b, k, procnum);
    end = omp_get_wtime();

    printf("Lancosh ended!\n");


    printf("Lancosh time: %lg\n", end - start);

    /*
    printf("\na:\n");
    print_vector(&a);
    printf("\nb:\n");
    print_vector(&b);
    */

    QR_values(&a, &b, &c, eps);
    printf("\nValues:\n");
    print_cvector(&c);
    delete_cvector(&c);
    delete_vector(&a);
    delete_vector(&b);
    delete_matrix(&matrix);
    printf("\nProgram ends!\n");
    return 0;
}