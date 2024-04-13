#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdlib.h>
#include <stdbool.h>
#include "spkmeans.h"
#include <math.h>

/** kmeans_pp is a function that gets as an input an array of points,
 * clusters them into k clusters and returns the clusters.*/
static struct cluster* kmeans_pp(struct point points[],double** centroids,int k,int size,double epsilon,int max_iter){
    int i,j;
    struct cluster *clusters ;
    double** oldCentroids ;
    double** newCentroids;
    bool same=true;
    clusters=malloc(sizeof(struct cluster)*k);
    for ( i = 0; i < k; i++) {
        clusters[i].index = i;
        clusters[i].centroid = centroids[i];
        clusters[i].clusterLen = 0;
    }
    for ( i = 0; i < size; i++) {
        AssignToClosestCluster(clusters, &points[i], k);
    }
    oldCentroids=malloc(sizeof(double*)*k);
    for ( i = 0; i < k; i++) {
        oldCentroids[i] = clusters[i].centroid;
    }
    newCentroids=malloc(sizeof(double*)*k);
    while(max_iter!=0 && same==true) {
        same = false;
        updateCentroid(clusters, k, points, size);
        for (i = 0; i < k; i++) {
            newCentroids[i] = clusters[i].centroid;
            clusters[i].clusterLen = 0;
        }
        for (i = 0; i < size; i++) {
            AssignToClosestCluster(clusters, &points[i], k);
        }
        for (i = 0; i < k; i++) {
            for (j = 0; j < points[0].d; j++) {
                if (oldCentroids[i][j] != newCentroids[i][j]) {
                    same = true;
                }
            }
        }
        for (i = 0; i < k; i++) {
            oldCentroids[i] = newCentroids[i];
        }
        max_iter--;
    }
    free(newCentroids);
    free(oldCentroids);
    return clusters;
}

/**getPyPoints is a function that gets as an input a python array of points
 * and converts them into an array of points from type struct point. */
struct point* getPyPoints(PyObject *_pointsfrompy,int n,int dimension){
    struct point *points;
    int i,j;
    double* coords;
    PyObject *arr;
    double x;
    points=(struct point*)malloc(sizeof(struct point)*n);
    if(!points){
        return NULL;
    }
    for(i=0;i<n;i++){
        coords=(double*)malloc(sizeof(double)*dimension);
        if(!coords) {
            return NULL;
        }
        arr=PyList_GetItem(_pointsfrompy,i);
        for(j=0;j<dimension;j++){
            x=PyFloat_AsDouble(PyList_GetItem(arr, j));
            coords[j]=x;
        }
        points[i].d=dimension;
        points[i].cords=coords;
    }
    return points;
}

/** convertMatToPy is a function that gets as an input a matrix and coverts it to a Python Object
 */
static PyObject* convertMatToPy(double** mat,int rows,int cols){
    PyObject *py_final_cen,*py_cen;
    int i,j;
    py_final_cen=PyList_New(rows);
    if(py_final_cen==NULL){
        return NULL;
    }
    for(i=0;i<rows;i++){
        py_cen=PyList_New(cols);
        if(py_cen==NULL){
            return NULL;
        }
        for(j=0;j<cols;j++){
            PyList_SetItem(py_cen,j,PyFloat_FromDouble(mat[i][j]));
        }
        PyList_SetItem(py_final_cen,i,py_cen);
    }
    return py_final_cen;
}


/** fitoperation is a function that get the points,k,n,d and the goal read from the python interface
 * and returns the corresponding matrix to the specified goal*/
static PyObject* fitoperation(PyObject* self,PyObject* args) {
    PyObject *_Pypoints;
    struct point* points;
    int k, n=0, d, goal,i,j;
    PyObject  *Res_matrix;
    double* eigenvalues;
    double **Wmatrix, **Dmatrix, **Lmatrix, **Jmatrix, **S_mat, **res_matrix;
    double* row;
    if (!PyArg_ParseTuple(args, "Oiiii", &_Pypoints, &k, &n, &d, &goal)) {
        exit(0);
    }
    Res_matrix=PyList_New(n);
    points = getPyPoints(_Pypoints, n, d);
    row=(double*)malloc(n*sizeof(double));
    Wmatrix= WeightAdjanecyMatrixFunc(points,n);
    Dmatrix= DiagonalDegreeMatrixFunc(Wmatrix,n);
    Lmatrix= lnormFunc(Wmatrix,Dmatrix,n);
    Jmatrix=(double**)malloc(sizeof(double*)*n);
    if (goal == 0) {
        Res_matrix = convertMatToPy(Wmatrix,n,n);
    }
    if(goal==1){
        Res_matrix= convertMatToPy(Dmatrix,n,n);
    }
    if(goal==2){
        Res_matrix= convertMatToPy(Lmatrix,n,n);
    }
    if(goal==3){
        eigenvalues=(double*)malloc(sizeof(double)*n);
        S_mat=build_jaccobi_mat(points,n);
        if(check_if_symmetric(S_mat,n)==0) {
            Jmatrix = jacobiFunc(S_mat, eigenvalues, n);
            res_matrix=(double**)malloc(sizeof(double*)*(n+1));
            if(!res_matrix){
                printf("An Error Has Occured");
                exit(0);
            }
            row=(double*)malloc(n*sizeof(double));
            if(!row) {
                printf("An Error Has Occured");
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
                    printf("An Error Has Occured");
                    exit(0);
                }
                for(j=0;j<n;j++){
                    row[j]=Jmatrix[i][j];
                }
                res_matrix[i+1]=row;
            }
            Res_matrix= convertMatToPy(res_matrix,n+1,n);
            freematrix(res_matrix,n+1);
            freematrix(Jmatrix,n);
            free(eigenvalues);
            freematrix(S_mat,n);
        }else{
            printf("Invalid Input!");
            exit(0);
        }
    }
    for(i=n-1;i>=0;i--){
        free(points[i].cords);
    }
    free(points);
    freematrix(Wmatrix,n);
    freematrix(Dmatrix,n);
    freematrix(Lmatrix,n);
    return Res_matrix;
}

/**fit_getTmat is a function that returns the T-matrix as a python Object.
 */
static PyObject* fit_getTmat(PyObject *self,PyObject *args) {
    double **Lmatrix,**matrix, **Umatrix, **Tmatrix;
    double **Wmatrix,**Dmatrix;
    double *arr, *eigenarr;
    int i,j,n,d,k;
    PyObject *py_final_cen,*py_cen, *_Pypoints;
    struct point *points;
    if (!PyArg_ParseTuple(args,"Oiii",&_Pypoints,&k,&n,&d)) {
        return NULL;
    }
    points= getPyPoints(_Pypoints,n,d);
    eigenarr=(double*)malloc(n*sizeof(double));
    Wmatrix= WeightAdjanecyMatrixFunc(points,n);
    Dmatrix= DiagonalDegreeMatrixFunc(Wmatrix,n);
    Lmatrix= lnormFunc(Wmatrix,Dmatrix,n);
    matrix= jacobiFunc(Lmatrix,eigenarr,n);
    arr=(double*)malloc(n*sizeof(double));
    for(i=0;i<n;i++){
        arr[i]=eigenarr[i];
    }
    qsort(arr,n,sizeof(double),compare);
    if(k==0){
        k= find_k(arr,n);
        if(k==1){
            printf("An Error Has Occurred");
            exit(0);
        }
    }
    Umatrix= getUMatrix(matrix,k,n,eigenarr,arr);
    Tmatrix= getTMatrix(Umatrix,n,k);
    py_final_cen=PyList_New(n);
    if(py_final_cen==NULL){
        return NULL;
    }
    for(i=0;i<n;i++){
        py_cen=PyList_New(k);
        if(py_cen==NULL){
            return NULL;
        }
        for(j=0;j<k;j++){
            PyList_SetItem(py_cen,j,PyFloat_FromDouble(Tmatrix[i][j]));
        }
        PyList_SetItem(py_final_cen,i,py_cen);
    }
    freematrix(Umatrix,n);
    freematrix(Wmatrix,n);
    freematrix(Dmatrix,n);
    freematrix(Lmatrix,n);
    freematrix(Tmatrix,n);
    freematrix(matrix,n);
    free(arr);
    free(eigenarr);
    for(i=n-1;i>=0;i--){
        free(points[i].cords);
    }
    free(points);
    return py_final_cen;

}


static PyObject* fit_spk(PyObject *self,PyObject *args){
    struct point* points;
    struct cluster* clusters;
    double** centroids;
    int k,n,d,i,j;
    PyObject *py_final_cen,*py_cen;
    PyObject *_Pypoints,*Pypoint;
    PyObject *_Pycentroid;
    double *point_coords;
    double x;
    if (!PyArg_ParseTuple(args,"OOiii",&_Pypoints,&_Pycentroid,&k,&n,&d)) {
        return NULL;
    }
    points= getPyPoints(_Pypoints,n,d);
    centroids=malloc(sizeof(double*)*k);
    for(i=0;i<k;i++){
        point_coords=malloc(sizeof(double)*d);
        Pypoint=PyList_GetItem(_Pycentroid,i);
        for(j=0;j<d;j++){
            x=PyFloat_AsDouble(PyList_GetItem(Pypoint, j));
            point_coords[j]=x;
        }
        centroids[i]=point_coords;
    }
    clusters= kmeans_pp(points,centroids,k,n,0,300);
    py_final_cen=PyList_New(k);
    if(py_final_cen==NULL){
        return NULL;
    }
    for(i=0;i<k;i++){
        py_cen=PyList_New(d);
        if(py_cen==NULL){
            return NULL;
        }
        for(j=0;j<points[0].d;j++){
            PyList_SetItem(py_cen,j,Py_BuildValue("d",clusters[i].centroid[j]));
        }
        PyList_SetItem(py_final_cen,i,Py_BuildValue("O",py_cen));
    }
    for(i=n-1;i>=0;i--){
        free(points[i].cords);
    }
    free(points);
    for(i=k-1;i>=0;i--){
        free(clusters[i].centroid);
    }
    free(clusters);
    return py_final_cen;
}
static PyMethodDef capiMethods[] = {
        { "fitoperation",
                (PyCFunction) fitoperation,
                METH_VARARGS,
                     PyDoc_STR("centroids for k clusters")},
        { "fit_getTmat",
                (PyCFunction)fit_getTmat,
                METH_VARARGS,
                     PyDoc_STR("centroids for k clusters")},
        { "fit_spk",
                (PyCFunction)fit_spk,
                METH_VARARGS,
                     PyDoc_STR("centroids for k clusters")},
        {NULL,NULL,0,NULL}
};

static struct PyModuleDef _moduledef= {
        PyModuleDef_HEAD_INIT,
        "spkmeansmod",
        NULL,
        -1,
        capiMethods
};

PyMODINIT_FUNC
PyInit_spkmeansmod(void)
{
    PyObject *m;
    m=PyModule_Create(&_moduledef);
    if(!m){
        return NULL;
    }
    return m;
}