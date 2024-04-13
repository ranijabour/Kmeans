#include "spkmeans.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "assert.h"
#include <string.h>
#include <math.h>
#include <ctype.h>



int main(int argc, char **argv){
    int i,j,n;
    char *inputFilename;
    int* numberOfPoints;
    struct point* points;
    double* eigenvalues;
    double* row;
    double **Wmatrix, **Dmatrix, **Lmatrix, **Jmatrix, **S_mat, **res_matrix;
    if (argc!=3){
        printf("Invalid Input!");
        return 0;
    }
    inputFilename=argv[2];
    numberOfPoints=malloc(sizeof(int));
    points=malloc(sizeof(struct point)*2000);
    if(!points){
        printf("An Error Has Occurred");
        exit(0);
    }
    ReadingPoints(inputFilename,&points,numberOfPoints);
    n=*numberOfPoints;
    Wmatrix= WeightAdjanecyMatrixFunc(points,n);
    Dmatrix= DiagonalDegreeMatrixFunc(Wmatrix,n);
    Lmatrix= lnormFunc(Wmatrix,Dmatrix,n);
    if(strcmp(argv[1],"spk")==0){
        return 0;
    }
    if (strcmp(argv[1], "wam") == 0){
        mat_print(Wmatrix,n,n);
    }
    if (strcmp(argv[1], "ddg") == 0){
        mat_print(Dmatrix,n,n);
    }
    if (strcmp(argv[1], "lnorm") == 0) {
        mat_print(Lmatrix,n,n);
    }
    if (strcmp(argv[1], "jacobi") == 0) {
        eigenvalues=(double*)malloc(sizeof(double)*n);
        S_mat=build_jaccobi_mat(points,n);
        if(check_if_symmetric(S_mat,n)==0){
            Jmatrix= jacobiFunc(S_mat,eigenvalues,n);
            res_matrix=(double**)malloc(sizeof(double*)*(n+1));
            if(!res_matrix){
                printf("An Error Has Occurred");
                exit(0);
            }
            row=(double*)malloc(n*sizeof(double));
            if(!row) {
                printf("An Error Has Occurred");
                exit(0);
            }
            for(j=0;j<n;j++){
                if(eigenvalues[j]>-0.0001 && eigenvalues[j]<0){
                    row[j]=fabs(eigenvalues[j]);
                }else{
                    row[j]=eigenvalues[j];
                }
            }
            res_matrix[0]=row;
            for(i=0;i<n;i++){
                row=(double*)malloc(n*sizeof(double));
                if(!row){
                    printf("An Error Has Occurred");
                    exit(0);
                }
                for(j=0;j<n;j++){
                    row[j]=Jmatrix[i][j];
                }

                res_matrix[i+1]=row;
            }
            mat_print(res_matrix,n+1,n);
            free(eigenvalues);
            freematrix(S_mat,n);
            freematrix(res_matrix,n+1);
            freematrix(Jmatrix,n);
        }else{
            printf("Invalid Input!");
            exit(0);
        }

    }
    free(numberOfPoints);
    for(i=n-1;i>=0;i--){
        free(points[i].cords);
    }
    free(points);
    freematrix(Wmatrix,n);
    freematrix(Dmatrix,n);
    freematrix(Lmatrix,n);
    return 0;
}

void ReadingPoints(char* filename,struct point** points,int* numberOfPoints){
    int initial_d=100,coord_index=0,point_index=0;
    double* curr_point;
    double* new_point;
    char curr_line[2000];
    char* splitted_line;
    FILE* inputFile = fopen(filename,"r");
    if(inputFile==NULL){
        printf("Invalid Input!");
        exit(0);
    }
    curr_point = (double*)malloc(sizeof(double)*initial_d);
    if (!curr_point){
        printf("An Error Has Occurred");
        exit(0);
    }
    while(fgets(curr_line,2000,inputFile)){
        splitted_line = strtok(curr_line,",");
        while(splitted_line!=NULL){
            if(coord_index==initial_d-1){
                curr_point= realloc(curr_point,sizeof(double)*2*coord_index);
                if(!curr_point){
                    printf("An Error Has Occurred");
                    exit(0);
                }
            }
            curr_point[coord_index]=atof(splitted_line);
            coord_index++;
            splitted_line= strtok(NULL,",");
        }
        if(point_index==0){
            new_point = (double*)realloc(curr_point,sizeof(double)*coord_index);
            if(!new_point){
                printf("An Error Has Occurred");
                exit(0);
            }
            (*points)[point_index].cords=new_point;
        }else{
            (*points)[point_index].cords=curr_point;
        }
        (*points)[point_index].d=coord_index;
        point_index++;
        curr_point= (double*)malloc(sizeof(double)*coord_index);
        if(!curr_point){
            printf("An Error Has Occurred");
            exit(0);
        }
        coord_index=0;
    }
    (*points)= realloc(*points,point_index*sizeof(struct point));
    if(!(*points)){
        printf("An Error Has Occurred");
        exit(0);
    }
    *numberOfPoints=point_index;
    fclose(inputFile);
    free(curr_point);
    free(splitted_line);
}

void freematrix(double** mat,int size){
    int i=0;
    for(i=0;i<size;i++){
        free(mat[i]);
    }
    free(mat);
}

double** Matrixsubtraction(double** matrix1,double** matrix2,int n){
    double** res;
    double* row;
    int i,j;
    res=(double**)malloc(n*sizeof(double*));
    if(!res){
        printf("An Error Has Occurred");
        exit(0);
    }
    for(i=0;i<n;i++) {
        row=(double*)malloc(n*sizeof(double));
        if(!row){
            printf("An Error Has Occurred");
            exit(0);
        }
        for (j = 0; j < n; j++) {
            row[j] = matrix1[i][j]-matrix2[i][j];
        }
        res[i]=row;
    }
    return res;
}

double** Matrixmultiplication(double** matrix1,double** matrix2,int n){
    double** res;
    double* row;
    int i,j,m;
    res=(double**)malloc(n*sizeof(double*));
    if(!res){
        printf("An Error Has Occurred");
        exit(0);
    }
    for(i=0;i<n;i++) {
        row = (double *) malloc(n * sizeof(double));
        if (!row) {
            printf("An Error Has Occurred");
            exit(0);
        }
        for (j = 0; j < n; j++) {
            row[j] = 0;
        }
        res[i]=row;
    }
    for(i=0;i<n;i++) {
        for (j = 0; j < n; j++) {
            res[i][j]=0;
            for (m = 0; m < n; m++) {
                res[i][j] += matrix1[i][m]*matrix2[m][j];
            }
        }
    }
    return res;
}
double** OneoverMatrixroot(double** matrix,int n){
    double** res;
    double* row;
    int i,j;
    res=(double**)malloc(n*sizeof(double*));
    if(!res){
        printf("An Error Has Occurred");
        exit(0);
    }
    for(i=0;i<n;i++){
        row=(double*)malloc(n*sizeof(double));
        if(!row){
            printf("An Error Has Occurred");
            exit(0);
        }
        for (j = 0; j < n; j++) {
            if (i == j) {
                row[j] = 1 / sqrt(matrix[i][i]);
            }
            else{
                row[j]=0;
            }
        }
        res[i]=row;
    }
    return res;
}

double** Imatrix(int n){
    double** imatrix;
    double* irow;
    int i,j;
    imatrix=(double**)malloc(n*sizeof(double*));
    if(!imatrix){
        printf("An Error Has Occurred");
        exit(0);
    }
    for(i=0;i<n;i++) {
        irow=(double*)malloc(n*sizeof(double));
        if(!irow){
            printf("An Error Has Occurred");
            exit(0);
        }
        for (j = 0; j < n; j++) {
            if(i==j){
                irow[j]=1;
            }
            else{
                irow[j]=0;
            }
        }
        imatrix[i]=irow;
    }
    return imatrix;
}

double** WeightAdjanecyMatrixFunc(struct point* points, int numOfPoint){
    double** WeightMatrix;
    double* ith_row;
    int i,j,d;
    double sum=0,num;
    WeightMatrix=(double**)malloc(numOfPoint*sizeof(double*));
    if(!WeightMatrix){
        printf("An Error Has Occurred");
        exit(0);
    }
    for(i=0;i<numOfPoint;i++){
        ith_row=(double*)malloc(numOfPoint*sizeof(double));
        if(!ith_row){
            printf("An Error Has Occurred");
            exit(0);
        }
        for(j=0;j<numOfPoint;j++) {
            if (j < i) {
                ith_row[j] = WeightMatrix[j][i];
            } else {
                if (j == i) {
                    ith_row[j] = 0.0;
                } else {
                    for (d = 0; d < points[i].d; d++) {
                        num = points[i].cords[d] - points[j].cords[d];
                        sum += pow(num, 2);
                    }
                    ith_row[j] = exp((-1 * (sqrt(sum))) / 2);
                    sum = 0;
                }
            }
        }
        WeightMatrix[i]=ith_row;
    }
    return WeightMatrix;
}

double** DiagonalDegreeMatrixFunc(double** WeightedMatrix, int numOfPoint){
    double** DiagonalMatrix;
    int i,j,z;
    double sum=0;
    double* ith_row;
    DiagonalMatrix=(double**)malloc(numOfPoint*sizeof(double*));
    if(!DiagonalMatrix){
        printf("An Error Has Occurred");
        exit(0);
    }
    for(i=0;i<numOfPoint;i++){
        ith_row=(double*)malloc(numOfPoint*sizeof (double));
        if(!ith_row){
            printf("An Error Has Occurred");
            exit(0);
        }
        for(j=0;j<numOfPoint;j++){
            if(i!=j){
                sum=0.0;
            }else{
                for(z=0;z<numOfPoint;z++){
                    sum+=WeightedMatrix[i][z];
                }
            }
            ith_row[j]=sum;
            sum=0.0;
        }
        DiagonalMatrix[i]=ith_row;
    }
    return DiagonalMatrix;
}

double** lnormFunc(double** Wmatrix,double** Dmatrix,int n){
    double** newDmatrix;
    double** imatrix;
    double** tmp;
    double** tmp2;
    double** res;
    imatrix = Imatrix(n);
    newDmatrix= OneoverMatrixroot(Dmatrix,n);
    tmp= Matrixmultiplication(newDmatrix,Wmatrix,n);
    tmp2= Matrixmultiplication(tmp,newDmatrix,n);
    res= Matrixsubtraction(imatrix,tmp2,n);
    freematrix(imatrix,n);
    freematrix(newDmatrix,n);
    freematrix(tmp,n);
    freematrix(tmp2,n);
    return res;
}

double** getRotationMatrix_P(double** A_matrix,double* pointer_to_c,double* pointer_to_s,int numOfPoints,int* i,int* j){
    int row=0,col=0,sign=-1,x=0,y=0;
    double max=0,A_ii=0,A_jj=0,theta=0,t=0,c=0,s=0;
    double** P_mat;
    double* mat_row;
    P_mat=(double**)malloc(numOfPoints*sizeof(double*));
    if(!P_mat){
        printf("An Error Has Occurred");
        exit(0);
    }
    /*finding A(i,j) such that it's the off-diagonal element with the largest absolute value*/
    for(row=0;row<numOfPoints;row++){
        for(col=row+1;col<numOfPoints;col++){
            if(fabs(A_matrix[row][col])> fabs(max)){
                max=A_matrix[row][col];
                *i=row;
                *j=col;
            }
        }
    }
    A_ii=A_matrix[*i][*i];
    A_jj=A_matrix[*j][*j];
    theta=(A_jj-A_ii)/(2*max);
    if(theta>=0){
        sign=1;
    }
    /*calculating c and s values*/
    t=(sign)/(fabs(theta)+sqrt((pow(theta,2)+1)));
    *pointer_to_c=1/sqrt(pow(t,2)+1);
    c=*pointer_to_c;
    s=(t*c);
    *pointer_to_s=s;
    for(x=0;x<numOfPoints;x++){
        mat_row=(double*)malloc(numOfPoints*sizeof(double));
        if(!mat_row){
            printf("An Error Has Occurred");
            exit(0);
        }
        for(y=0;y<numOfPoints;y++){
            if(x==y){
                mat_row[y]=1;
            }else{
                mat_row[y]=0;
            }
        }
        P_mat[x]=mat_row;
    }
    P_mat[*i][*i]=*pointer_to_c;
    P_mat[*j][*j]=*pointer_to_c;
    P_mat[*i][*j]=*pointer_to_s;
    P_mat[*j][*i]=-(*pointer_to_s);
    return P_mat;
}
double** getTransposeMatrix(double** Matrix, int n){
    double** transposed;
    double* row;
    int i,j;
    transposed=(double**)malloc(n*sizeof(double*));
    if(!transposed){
        printf("An Error Has Occurred");
        exit(0);
    }
    for(i=0;i<n;i++){
        row=(double*)malloc(n*sizeof(double));
        if(!row){
            printf("An Error Has Occurred");
            exit(0);
        }
        transposed[i]=row;
    }
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            transposed[j][i]=Matrix[i][j];
        }
    }
    return transposed;
}

double** jacobiFunc(double** A,double* eigenValues,int numOfPoints){
    double **P,**V,**A_tag,**V2;
    double *pointer_to_c, *pointer_to_s;
    int *i,*j;
    int iterations=0,row,col;
    pointer_to_s=(double*)malloc(sizeof(double));
    pointer_to_c=(double*)malloc(sizeof(double));
    i=(int*)malloc(sizeof(int));
    j=(int*)malloc(sizeof(int));
    if((!pointer_to_s) ||(!pointer_to_c)||(!i)||(!j)){
        printf("An Error Has Occurred");
        exit(0);
    }
    P = getRotationMatrix_P(A,pointer_to_c,pointer_to_s,numOfPoints,i,j);
    A_tag = getAtag(A,*i,*j,numOfPoints,*pointer_to_c,*pointer_to_s);
    V =getRotationMatrix_P(A,pointer_to_c,pointer_to_s,numOfPoints,i,j);
    V2 =getRotationMatrix_P(A,pointer_to_c,pointer_to_s,numOfPoints,i,j);
    freematrix(V, numOfPoints);
    while(((off_mat(A,numOfPoints)- off_mat(A_tag,numOfPoints))>(1.0*pow(10,-5))) && iterations<99){
        for(row=0;row<numOfPoints;row++){
            for(col=0;col<numOfPoints;col++){
                A[row][col]=A_tag[row][col];
            }
        }
        freematrix(P,numOfPoints);
        freematrix(A_tag,numOfPoints);
        P= getRotationMatrix_P(A,pointer_to_c,pointer_to_s,numOfPoints,i,j);
        A_tag= getAtag(A,*i,*j,numOfPoints,*pointer_to_c,*pointer_to_s);
        if(iterations%2==0){
            V= Matrixmultiplication(V2,P,numOfPoints);
            freematrix(V2,numOfPoints);
        }else{
            V2= Matrixmultiplication(V,P,numOfPoints);
            freematrix(V,numOfPoints);
        }
        iterations++;
    }
    for(row=0;row<numOfPoints;row++){
        eigenValues[row]=A_tag[row][row];
    }
    free(i);
    free(j);
    free(pointer_to_s);
    free(pointer_to_c);
    freematrix(P,numOfPoints);
    freematrix(A_tag,numOfPoints);
    if(iterations%2!=0){
        return V;
    }
    return V2;
}

double** getAtag(double** A,int i,int j ,int n,double c,double s){
    double** A_tag;
    int row=0,col=0;
    double* A_tag_row;
    A_tag=(double**)malloc(n*sizeof(double*));
    if(!A_tag){
        printf("An Error Has Occurred");
        exit(0);
    }
    /*copying the content of A to A_tag*/
    for(row=0;row<n;row++){
        A_tag_row=(double*)malloc(n*sizeof(double));
        if(!A_tag_row){
            printf("An Error Has Occurred");
            exit(0);
        }
        for(col=0;col<n;col++){
            A_tag_row[col]=A[row][col];
        }
        A_tag[row]=A_tag_row;
    }
    A_tag[i][i]=c*c*A[i][i]+s*s*A[j][j]-2*s*c*A[i][j];
    A_tag[j][j]=s*s*A[i][i]+c*c*A[j][j]+2*s*c*A[i][j];
    A_tag[i][j]=0.0;
    A_tag[j][i]=0.0;
    for(row=0;row<n;row++){
        if((row!=i) &&(row!=j)){
            A_tag[row][i]=c*A[row][i]-s*A[row][j];
            A_tag[row][j]=c*A[row][j]+s*A[row][i];
            A_tag[i][row]=A_tag[row][i];
            A_tag[j][row]=A_tag[row][j];
        }
    }
    return A_tag;
}

double off_mat(double** Matrix,int n){
    double res=0.0;
    int i,j;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i!=j){
                res+=pow(Matrix[i][j],2);
            }
        }
    }
    return res;
}


void mat_print(double** Matrix,int rows,int cols){
    int i,j;
    for(i=0;i<rows;i++) {
        for (j = 0; j < cols; j++) {
            if (j < cols - 1) {
                printf("%.4f,", Matrix[i][j]);
            } else {
                printf("%.4f\n", Matrix[i][j]);
            }
        }
    }
}

double** build_jaccobi_mat(struct point* points,int n){
    double** res_Mat;
    int i,j,d;
    double* row;
    res_Mat=(double**)malloc(n*sizeof(double*));
    if(!res_Mat){
        printf("An Error Has Occurred");
        exit(0);
    }
    d=points[0].d;
    for(i=0;i<n;i++){
        row=(double*)malloc(sizeof(double)*n);
        if(!row){
            printf("An Error Has Occurred");
            exit(0);
        }
        for(j=0;j<n;j++){
            if(j<d){
                row[j]=points[i].cords[j];
            }else{
                row[j]=0.0;
            }
        }
        res_Mat[i]=row;
    }
    return res_Mat;
}

/*the following func returns 0 if the input Matrix is symmetric, 1 otherwise */
int check_if_symmetric(double** mat,int n){
    int i,j;
    for(i=0;i<n;i++){
        for(j=i+1;j<n;j++){
            if(mat[i][j]!=mat[j][i]){
                return 1;
            }
        }
    }
    return 0;
}

int find_k(double* eigenValues,int n){
    int i,index=0;
    double max=0.0;
    for(i=0;i<floor(n/2);i++){
        if(fabs(eigenValues[i]-eigenValues[i+1])>max){
            max=fabs(eigenValues[i]-eigenValues[i+1]);
            index=i;
        }
    }
    return index+1;
}

int compare(const void *a,const void *b){
    double* x;
    double* y;
    x=(double*)a;
    y=(double*)b;
    if(*x>*y){
        return 1;
    }else if(*x<*y){
        return -1;
    }else{
        return 0;
    }
}

double** getUMatrix(double** matrix,int k, int n,double* eigen,double* sorted_eigen){
    int i,j,r,index;
    double** U_mat;
    double* row;
    U_mat=(double**)malloc(sizeof(double*)*n);
    if(!U_mat){
        printf("An Error Has Occurred");
        exit(0);
    }
    for(i=0;i<n;i++){
        row=(double*)malloc(sizeof(double)*k);
        if(!row){
            printf("An Error Has Occurred");
            exit(0);
        }
        U_mat[i]=row;
        for(j=0;j<k;j++){
            row[j]=0.0;
        }
    }
    for(i=0;i<k;i++){
        for(j=0;j<n;j++){
            if(eigen[j]==sorted_eigen[i]){
                /*finding the corresponding vector index*/
                index=j;
                break;
            }
        }
        for(r=0;r<n;r++){
            U_mat[r][i]=matrix[r][index];
        }
    }
    return U_mat;
}

double** getTMatrix(double** U_mat, int n, int k){
    double** T_mat;
    double* row;
    int i,j;
    double row_sum=0;
    T_mat=(double**)malloc(sizeof(double*)*n);
    if(!T_mat){
        printf("An Error Has Occurred");
        exit(0);
    }
    for(i=0;i<n;i++){
        row=(double*)malloc(sizeof(double)*k);
        if(!row){
            printf("An Error Has Occurred");
            exit(0);
        }
        for(j=0;j<k;j++){
            row[j]=0.0;
        }
        T_mat[i]=row;
    }
    for(i=0;i<n;i++){
        for(j=0;j<k;j++){
            row_sum+=pow(U_mat[i][j],2);
        }
        for(j=0;j<k;j++){
            if(row_sum==0){
                T_mat[i][j]=0.0;
            }else{
                T_mat[i][j]=U_mat[i][j]/ pow(row_sum,0.5);
            }
        }
        row_sum=0.0;
    }
    return T_mat;
}

void AssignToClosestCluster(struct cluster clusters[], struct point *Point, int k) {
    struct cluster* closestCluster;
    int dimension=Point->d,j,i;
    double* cordes=Point->cords;
    double* clusterCentroids=clusters[0].centroid;
    double sum=0,num;
    double min;
    for ( j = 0; j < dimension; j++) {
        num = cordes[j] - clusterCentroids[j];
        sum += pow(num, 2);
    }
    min=sum;
    closestCluster = &(clusters[0]);
    if(min==0){
        Point->cluster_index=closestCluster->index;
        closestCluster->clusterLen++;
        return;
    }
    for ( i = 1; i < k; i++) {
        cordes=Point->cords;
        clusterCentroids=clusters[i].centroid;
        sum=0;
        for ( j = 0; j < dimension; j++) {
            num = cordes[j] - clusterCentroids[j];
            sum += pow(num, 2);
        }
        if(sum<min){
            min = sum;
            closestCluster = &(clusters[i]);
        }
    }
    Point->cluster_index=closestCluster->index;
    closestCluster->clusterLen++;
}

void updateCentroid(struct cluster c[],int k,struct point points[],int numberOfLines) {
    int dimension=points[0].d;
    int i,y,j,z;
    double sum;
    double *arr;
    double* p1;
    for( i=0;i<k;i++) {
        for ( y= 0; y < dimension; y++) {
            sum = 0;
            for ( j = 0; j < numberOfLines; j++) {
                if (points[j].cluster_index == i) {
                    arr = malloc(sizeof(double) * dimension+1);
                    if(!arr){
                        printf("An Error Has Occurred");
                        exit(0);
                    }
                    p1 = points[j].cords;
                    for ( z = 0; z < dimension ; z++) {
                        arr[z] = p1[z];
                    }
                    sum += arr[y];
                    free(arr);
                }
            }
            c[i].centroid[y] = sum / c[i].clusterLen;
        }
    }
}
