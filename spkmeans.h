
#ifndef FINALPROJECT_SPKMEANS_H
#define FINALPROJECT_SPKMEANS_H
typedef enum {spk, wam, ddg, lnorm, jacobi
} Goal;

struct point
{
    double* cords;
    int cluster_index;
    int d;
};
typedef struct point point;

struct cluster
{
    double* centroid;
    int index;
    int clusterLen;
};
typedef struct cluster cluster;

int main(int argc,char **argv) ;

void ReadingPoints(char* filename,struct point** points,int* numberOfPoints);

/**
 * freematrix is a function that gets a matrix as an input and frees the memory that it holds.
 */

void freematrix(double** mat,int size);

/**
 * MatrixSubstraction is a function that takes two matrixes and returns a matrix
 * which is the result of substracting the second input matrix from the first input matrix .
 */

double** Matrixsubtraction(double** matrix1,double** matrix2,int n);

/**
 * Matrixmultiplication is a function that takes two matrixes and returns a matrix which is
 * the result of multiplying the second input matrix in the first input matrix.
 */

double** Matrixmultiplication(double** matrix1,double** matrix2,int n);

/**
 * OneoverMatrixroot is a function that gets a matrix as an input and calculates the matrix D^-0.5.
 */

double** OneoverMatrixroot(double** matrix,int n);

/**
 * Imatrix is a function that gets an integer as input and returns the I-matrix of size n*n.
 */

double** Imatrix(int n);

/**
 * WeightAdjanecyMatrixFunc is a function that gets as an input an array of points and the
 * number of points, calculates and returns the Weight matrix as specified in the instructions.
 */
double** WeightAdjanecyMatrixFunc(struct point* points, int numOfPoint);

/**
 * DiagonalDegreeMatrixFunc is a function that gets as an input the weight matrix and the number
 * of points(an int), calculates and returns the Diagonal matrix.
 */

double** DiagonalDegreeMatrixFunc(double** WeightedMatrix, int numOfPoint);

/**
 * Lnorm is a function that calculates and returns the Normalized laplacian matrix.
 */

double** lnormFunc(double** Wmatrix,double** Dmatrix,int n);

/**
 * getRotationMatrix_P is a function that calculates the value of c and s as specified in the instructions,
 * and returns the Rotation Matrix.
 */

double** getRotationMatrix_P(double** A_matrix,double* pointer_to_c,double* pointer_to_s,int numOfPoints,int* i,int* j);

double** getTransposeMatrix(double** Matrix, int n);

/**
 * jacobiFunc is a function that returns the V matrix and calculates the sigenValues and eigen vectors of a given matrix as an input
 */

double** jacobiFunc(double** A,double* eigenValues,int numOfPoints);

/**
 * getAtag is function that gets A matrix, calculates and returns the A tag matrix.
 */

double** getAtag(double** A,int i,int j ,int n,double c,double s);

/**
 * off_mat is a function that gets a matrix and calculates the sum of squares of all off-diagonal elements.
 */

double off_mat(double** Matrix,int n);

 /**
  * mat_print gets a matrix and prints it.
  */

void mat_print(double** Matrix,int rows,int cols);

/**
 *  build_jaccobi_mat gets as an input an array of points of type struct point,and builds a matrix from the points
 *  such that the ith_row contains the coordinates of the ith point of the array.
 */

double** build_jaccobi_mat(struct point* points,int n);

/**
 *the following func returns 0 if the input Matrix is symmetric, 1 otherwise.

*/
int check_if_symmetric(double** mat,int n);

/**
 * find_k is a function that gets an array of the eigen values sorted in increasing order,
 * calculates the optimal k(number of clusters) and returns it.
*/

int find_k(double* eigenValues,int n);

/**
 * compare is a function that gets a an b as an inputs, return 1 if a>b, returns -1 if b>a and 0 if a==b.
 */

int compare(const void *a,const void *b);

/**
 * getUMatrix is a function that calculates the U-matrix as specified in the instructions and returns it.
 */

double** getUMatrix(double** matrix,int k, int n,double* eigen,double* sorted_eigen);

/**
 * getTMatrix is a function that calculates the T-matrix as specified in the instructions and returns it.
 */

double** getTMatrix(double** U_mat, int n, int k);

void AssignToClosestCluster(struct cluster clusters[], struct point *Point, int k) ;

void updateCentroid(struct cluster c[],int k,struct point points[],int numberOfLines) ;

#endif
