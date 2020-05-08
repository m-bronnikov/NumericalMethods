#include "parallel_structs/RVector.h"
#include "parallel_structs/CVector.h"
#include "parallel_structs/RMatrix.h"


int main(int argc, char *argv[]){
    char filename[] = "data/filen.txt";
    FILE *fout = fopen("data/output.log", "w");

    for(char c = '3'; c < '9'; ++c){
        filename[9] = c;
        RMatrix matrix = read_matrix(filename);
        RVector a = create_vector(), b = create_vector();

        printf("%s started\n", filename);

        for(unsigned procnum = 1; procnum < 10; ++procnum){
            double end, start = omp_get_wtime();
            parallel_lancosh_method(&matrix, &a, &b, 25, procnum);
            end = omp_get_wtime();
            fprintf(fout, "%u %u %lg\n", matrix_size(&matrix), procnum, end - start);
        }

        printf("%s done!\n", filename);

        delete_vector(&a);
        delete_vector(&b);
        delete_matrix(&matrix);
    }

    fclose(fout);
    return 0;
}