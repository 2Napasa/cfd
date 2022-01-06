#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int N = 200, M = 200;
double epsilon = 1e-10;

void writeDoubleVectorToFile(const char* fileName, double *matToPrint, int rows){
    FILE *f = fopen(fileName, "w+");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        
    }
    for (int i = 0; i < rows; i++){
        fprintf(f,"%d\t",i);
        fprintf(f,"%12.5f",matToPrint[i]);
        fprintf(f,"\n");
    }
    fclose(f);
}
void writeIntegerVectorToFile(const char* fileName, int *matToPrint, int rows){
    FILE *f = fopen(fileName, "w+");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        
    }

    for (int i = 0; i < rows; i++){
        fprintf(f,"%d\t",i);
        fprintf(f,"%12d",matToPrint[i]);
        fprintf(f,"\n");
    }
    fclose(f);
}
void writeDoubleMatrixToFile(const char* fileName, double **matToPrint, int rows, int cols){
    FILE *f = fopen(fileName, "w+");
    if (f == NULL)
    {
        printf("Error opening file!\n");
    }
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols-1; j++){
                fprintf(f,"%f,",matToPrint[i][j]);
        }
        int j = cols - 1;
        fprintf(f,"%f",matToPrint[i][j]);
        fprintf(f,"\n");
    }
    fclose(f);
}
void writeIntegerMatrixToFile(const char* fileName, int **matToPrint, int rows, int cols){
    FILE *f = fopen(fileName, "w+");
    if (f == NULL)
    {
        printf("Error opening file!\n");
    }
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
                fprintf(f,"%12d,",matToPrint[i][j]);
        }
        fprintf(f,"\n");
    }
    fclose(f);
}
void showIntegerMatrix(int **a, int rows, int cols){
    printf("Showing matrix of doubles:\n");
    for (int i = 0; i < rows; ++i){
        for (int j = 0; j < cols; ++j){
            printf("%12d",a[i][j]);
            // cout << a[i][j] << " ";
        }
    // cout<<endl;
    printf("\n");
    }
}
void showDoubleMatrix(double **a, int rows, int cols){
    printf("Showing matrix of int:\n");
    for (int i = 0; i < rows; ++i){
        for (int j = 0; j < cols; ++j){
            // cout << a[i][j] << " ";
            printf("%12.2f",a[i][j]);
        }
    // cout<<endl;
    printf("\n");
    }
}
void showDoubleVector(double *a, int rows){
    printf("Showing vector of doubles:\n");
    for (int i = 0; i < rows; ++i){
        printf("%d\t%12.4f",i,a[i]);
    printf("\n");
    }
}
void showIntegerVector(int *a, int rows){
    printf("Showing vector of int:\n");
    for (int i = 0; i < rows; ++i){
        printf("%d\t%12d",i,a[i]);
    printf("\n");
    }
}
double *cjCSRdense(double *val, int *col, int *row, int valSize, int colSize, int rowSize, double *b, int ansSize){
    double *x; 
    x = (double *) malloc(ansSize * sizeof(double));
    double *r;
    r = (double *) malloc(ansSize * sizeof(double));
    double *p;
    p = (double *) malloc(ansSize * sizeof(double));
    double *Ap;
    Ap = (double *) malloc(ansSize * sizeof(double));
    double rsold, alpha, rsnew, sum;
    int colStart = 0;
    int i, j, l, k, rowStart, rowEnd, colEnd;
    for (i = 0; i < ansSize; i++){// fill x with 0
        x[i] = 0;
    }
    for (i = 0; i < ansSize; i++){
        rowStart = row[i];
        rowEnd = row[i+1];
        sum = 0;
        for (j = rowStart; j < rowEnd; j++){
            k = col[j];
            sum = sum + val[j]*x[k];
        }
        r[i] = b[i] - sum;
    }
    for (i = 0; i < ansSize; i++){
        p[i] = r[i];
    }
    rsold = 0;
    for (i = 0; i < ansSize; i++){
        rsold = rsold + r[i] * r[i];
    }
    int ctr;
    for (i = 0; i < ansSize; i++){
        for (l = 0; l < ansSize; l++){
            rowStart = row[l];
            rowEnd = row[l+1];
            sum = 0;
            for (j = rowStart; j < rowEnd; j++){
                k = col[j];
                sum = sum + val[j]*p[k];
            }
            Ap[l] = sum;
        }
        sum = 0;
        for (j = 0; j < ansSize; j++){
            sum = sum + p[j] * Ap[j];
        }
        alpha = rsold / sum;
        for (j = 0; j < ansSize; j++){
            x[j] = x[j] + alpha * p[j];
            r[j] = r[j] - alpha * Ap[j];
        }
        rsnew = 0;
        for (j = 0; j < ansSize; j++){
            rsnew = rsnew + r[j] * r[j];
        }
        if (i%5000==0 && i > 1){
            printf("%5d/%d. residual = %.4e\n", i, ansSize, sqrt(rsnew));
        }
        if (sqrt(rsnew) < epsilon){
            printf("%5d/%d. residual = %.4e\n", ctr, ansSize, sqrt(rsnew));
            return x;
        }
        for (j = 0; j < ansSize; j++){
            p[j] = r[j] + (rsnew/rsold)  * p[j];
        }
        rsold = rsnew;
        ctr = i;
    }
    free(r);
    free(p);
    free(Ap);
    printf("%5d/%d. residual = %.4e\n", ctr, ansSize, sqrt(rsnew));
    return x;
}

int main(){
    double L = 1;
    int kmax = M*N;
    double dx = L/(N-1);
    double dy = L/(M-1);
    double dx2 = dx*dx; 
    double dy2 = dy*dy;
    // memory allocation 
    double *source; 
    source = (double *) malloc(N*M * sizeof(double));
    double *q; 
    q = (double *) malloc(N*M * sizeof(double));
        // sparse CSR/CRS/Yale format
    double *values; 
    values = (double *) malloc(N*M*5 * sizeof(double));
    int *col_index;
    col_index = (int *) malloc(N*M*5 * sizeof(int));
    int *row_index;
    row_index = (int *) malloc((N*M+1) * sizeof(int));
        // my own sparse format
    int **Acoord = (int**)malloc(kmax * sizeof(int*));
    double **Aval = (double**)malloc(kmax * sizeof(double*));
    for (int i = 0; i < kmax; i++){
        Acoord[i] = (int*)malloc(5 * sizeof(int));
        Aval[i] = (double*)malloc(5 * sizeof(double));
    }
    // end memory allocation

    // fill source rhs f function
    double sum = 0;
    for (int i = 0; i< N; i++){
        for (int j = 0; j < M ; j++){
            int k = j*N+i;
            source[k] = -exp(-pow(i*dx-L/2,2)-pow(j*dy-L/2,2));
            sum = sum + source[k];
        }
    }
    sum = sum/N/M;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < M; j++){
            int k = j*N+i;
            source[k] = source[k] - sum;
        }
    }

    // fill inside and boundary  
    {
    
    for (int i = 1; i < N-1; i++){
        for (int j = 1; j < M-1; j++){
            int k = j*N+i;
            Aval[k][2] = -(2/dx2+2/dy2);
            Aval[k][1] = 1/dx2;
            Aval[k][3] = 1/dx2;
            Aval[k][0] = 1/dy2;
            Aval[k][4] = 1/dy2;

            Acoord[k][0] = k - N;
            Acoord[k][1] = k - 1;
            Acoord[k][2] = k;
            Acoord[k][3] = k + 1;
            Acoord[k][4] = k + N;
            q[k] = source[k];
        }
    }
    
    
    // fill boundary part
    int i = N-1; // RIGHT
    double JR = 0; // FLUX
    for (int j = 1; j < M-1; j++){ // (fi+1):M-1
        int k = j*N+i;

        Aval[k][0] = 1/dy2;
        Aval[k][1] = 1/dx2;
        Aval[k][2] = -1/dx2-2/dy2;
        Aval[k][3] = 0;
        Aval[k][4] = 1/dy2;
        
        Acoord[k][0] = k - N;
        Acoord[k][1] = k - 1;
        Acoord[k][2] = k;
        Acoord[k][3] = -1;
        Acoord[k][4] = k + N;

        q[k] = source[k];
    }

    i = 0; // LEFT
    double JL = 0; // FLUX
    for (int j = 1; j < M-1; j++){
        int k = j*N+i;

        Aval[k][0] = 1/dy2;
        Aval[k][1] = 0;
        Aval[k][2] = -1/dx2-2/dy2;
        Aval[k][3] = 1/dx2;
        Aval[k][4] = 1/dy2;

        Acoord[k][0] = k - N;
        Acoord[k][1] = - 1;
        Acoord[k][2] = k;
        Acoord[k][3] = k + 1;
        Acoord[k][4] = k + N;
        
        q[k] = source[k];
    }

    int j = 0; // BOT
    double JB = 0; // FLUX
    for (int i = 1; i < N-1; i++){
        int k = j*N+i;

        Aval[k][0] = 0;
        Aval[k][1] = 1/dx2;
        Aval[k][2] = -2/dx2-1/dy2;
        Aval[k][3] = 1/dx2;
        Aval[k][4] = 1/dy2;

        Acoord[k][0] = -1;
        Acoord[k][1] = k - 1;
        Acoord[k][2] = k;
        Acoord[k][3] = k + 1;
        Acoord[k][4] = k + N;

        q[k] = source[k];
    }

    j = M-1; // TOP
    double JT = 0; // FLUX
    for (int i = 1; i < N-1; i++){
        int k = j*N+i;

        Aval[k][0] = 1/dy2;
        Aval[k][1] = 1/dx2;
        Aval[k][2] = -2/dx2-1/dy2;
        Aval[k][3] = 1/dx2;
        Aval[k][4] = 0;

        Acoord[k][0] = k - N;
        Acoord[k][1] = k - 1;
        Acoord[k][2] = k;
        Acoord[k][3] = k + 1;
        Acoord[k][4] = -1;

        q[k] = source[k];
    }

    i = 0; // LEFT
    j = 0; // BOT 
    int k = j*N+i;

    Aval[k][0] = 0;
    Aval[k][1] = 0;
    Aval[k][2] = -1/dx2 - 1/dy2;
    Aval[k][3] = 1/dx2;
    Aval[k][4] = 1/dy2;

    Acoord[k][0] = -1;
    Acoord[k][1] = -1;
    Acoord[k][2] = k;
    Acoord[k][3] = k + 1;
    Acoord[k][4] = k + N;

    q[k] = source[k]; 

    i = 0; // LEFT
    j = M-1; // TOP 
    k = j*N+i;

    Aval[k][0] = 1/dy2;
    Aval[k][1] = 0;
    Aval[k][2] = -1/dx2 - 1/dy2;
    Aval[k][3] = 1/dx2;
    Aval[k][4] = 0;

    Acoord[k][0] = k - N;
    Acoord[k][1] = -1;
    Acoord[k][2] = k;
    Acoord[k][3] = k + 1;
    Acoord[k][4] = -1;

    q[k] = source[k]; 

    i = N-1; // RIGHT
    j = M-1; // TOP 
    k = j*N+i;

    Aval[k][0] = 1/dy2;
    Aval[k][1] = 1/dx2;
    Aval[k][2] = -1/dx2 - 1/dy2;
    Aval[k][3] = 0;
    Aval[k][4] = 0;

    Acoord[k][0] = k - N;
    Acoord[k][1] = k - 1;
    Acoord[k][2] = k;
    Acoord[k][3] = -1;
    Acoord[k][4] = -1;

    q[k] = source[k]; 

    i = N-1; // RIGHT
    j = 0; // BOT 
    k = (j)*M+i;

    Aval[k][0] = 0;
    Aval[k][1] = 1/dx2;
    Aval[k][2] = -1/dx2 - 1/dy2;
    Aval[k][3] = 0;
    Aval[k][4] = 1/dy2;

    Acoord[k][0] = -1;
    Acoord[k][1] = k - 1;
    Acoord[k][2] = k;
    Acoord[k][3] = -1;
    Acoord[k][4] = k + N;

    q[k] = source[k]; 
    }
    
    // convert to CSR/CRS/Yale
    int valCtr = 0;
    int rowCtr = 0;
    int colCtr = 0;
    row_index[0] = 0;
    for (int i = 0; i<kmax; i++){
        for (int j = 0; j<5; j++){
            if (Aval[i][j] != 0) {
                values[valCtr] = Aval[i][j];
                valCtr++;                
            }
            if (Acoord[i][j] != -1) {
                col_index[colCtr] = Acoord[i][j];
                colCtr++;                
            }
        }
        row_index[rowCtr+1]  = valCtr;
        rowCtr++;
    }
    row_index[rowCtr] = valCtr;
    rowCtr++;
    
    // showDoubleVector(values, valCtr);
    // showIntegerVector(col_index, colCtr);
    // showIntegerVector(row_index, rowCtr);

    // writeIntegerVectorToFile("testcjCSR_rindex.txt", row_index, rowCtr);
    // writeIntegerVectorToFile("testcjCSR_cindex.txt", col_index, colCtr);
    // writeDoubleVectorToFile("testcjCSR_v.txt", values, valCtr);


    // cjCSRdense solver
    double *answer = cjCSRdense(values, col_index, row_index, valCtr, colCtr, rowCtr, q, kmax);
    double** ans2D = (double**)malloc(N * sizeof(double*));
    for (int i = 0; i < N; i++){
        ans2D[i] = (double*)malloc(M * sizeof(double));
        for (int j = 0; j < M; j++){
            int k = j * N + i;
            ans2D[i][j] = answer[k];
        }
    }
    
    writeDoubleMatrixToFile("testcjCSR_2Doutput.csv", ans2D, N, M);
    // cjCSRdense solver ends
    
    
    free(answer);
    free(source);
    free(q);
    free(values);
    free(row_index);
    free(col_index);
    
    for (int i = 0; i < N; i++)
        free(ans2D[i]);
    free(ans2D);

    for (int i = 0; i < kmax; i++)
        free(Aval[i]);
    free(Aval);
    for (int i = 0; i < kmax; i++)
        free(Acoord[i]);
    free(Acoord);
    return 0;
}