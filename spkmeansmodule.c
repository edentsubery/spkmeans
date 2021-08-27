
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "spkmeans.c"
#include "spkmeans.h"
#include "../../../Users/edent/AppData/Local/Programs/Python/Python36/include/object.h"
#include "../../../Users/edent/AppData/Local/Programs/Python/Python36/include/Python.h"

static Point *general(PyObject *args);

static Cluster *createClusterArrayWithCentroids(PyObject* pyCentroids) {
    int i,j;
    Cluster *arr;
    void *pointerForCentroid;
    void *pointerForSums;
    PyObject* pyList;
    arr = malloc(k * sizeof(Cluster));
    assert(arr!=NULL);
    for (i = 0; i < k; i++) {
        pointerForCentroid = malloc(d * sizeof(double));
        assert(pointerForCentroid!=NULL);
        pointerForSums = malloc(d * sizeof(double));
        assert(pointerForSums!=NULL);
        arr[i].centroid.coordinates = pointerForCentroid;
        pyList=PyList_GetItem(pyCentroids,i);
        for(j=0;j<d;j++){
            arr[i].centroid.coordinates[j]=PyFloat_AsDouble(PyList_GET_ITEM(pyList,j)) ;
        }
        arr[i].sum_by_coordinates.coordinates = pointerForSums;
        arr[i].size = 0;
    }
    return arr;
}

static void readAllPointsPy( Point *points,PyObject* pyPoints) {
    Point *curPoint;
    int index;
    int i;
    PyObject* pyList;
    index = 0;
    while (index < numberOfPoints) {
        curPoint = malloc(sizeof(Point));
        assert(curPoint!=NULL);
        curPoint->coordinates = malloc(d * sizeof(double));
        pyList=PyList_GetItem(pyPoints,index);
        assert(curPoint->coordinates!=NULL);
        for (i = 0; i < d; i++) {
            curPoint->coordinates[i] = PyFloat_AsDouble(PyList_GetItem(pyList,i));
        }
        addPointToArr(points, *curPoint, index++);
    }
}

static void goalFunction(PyObject *self, PyObject *args){
    Point* points=general(args);
    if(strcmp(goal,"wam") == 0){
        wam(points);
        printMatrix(WAM);
    }
    if(strcmp(goal,"ddg") == 0){
        ddm(points);
        printMatrix(DDM);
    }
    if(strcmp(goal,"lnorm") == 0){
        lnorm(points);
        printMatrix(LNORM);
    }
    if(strcmp(goal,"jacobi") == 0){
        jacobi(points);
        printEigenvalues();
        printMatrix(eigenvectors);
    }
    freePointPointer(points);
}


static void kmeanspp(PyObject *self, PyObject *args){
    Point *points1;
    Cluster *clusterArray1;
    Cluster *actualClusters;
    PyObject* pyCentroids;
    PyObject* pyPoints;
    if(!PyArg_ParseTuple(args,"ooiiii",&pyPoints, &pyCentroids,&k,&d,&MAX_ITER,&numberOfPoints)){
        printf("parsing failed\n");
    }
    points1=getPointPointer();
    actualClusters=createClusterArrayEmpty();
    readAllPointsPy(points1,pyPoints);
    clusterArray1 = createClusterArrayWithCentroids(pyCentroids);
    kmeansWithInitialCentroids(points1, clusterArray1, actualClusters);
    freeClusterArray(actualClusters);
    freeClusterArray(clusterArray1);
    freePointPointer(points1);
}

static double** getT(PyObject *self, PyObject *args){
    Point *points = general(args);
    spk(points);
    freePointPointer(points);
    printf("calcuated T");
    return T.values;
}

static Point *general(PyObject *args) {
    Point *points;
    char * filePath;
    PyArg_ParseTuple(args,"iss", &k,&goal,&filePath);
    FILE *file=fopen( filePath, "r" );
    countPoints(file);
    findD(file);
    points = getPointPointer();
    readAllPointsC(points,file);
    validateArguments();
    return points;
}


static PyMethodDef capiMethods[]={
        {"newDimension",
         (PyCFunction) getT,
         METH_VARARGS,
         PyDoc_STR("convert to a new dimension")},

        {"kmeanspp",
         (PyCFunction) kmeanspp,
         METH_VARARGS,
         PyDoc_STR("calculates k means with the initial centroids passed to it")},

        {"goal",
         (PyCFunction) goalFunction,
         METH_VARARGS,
         PyDoc_STR("calculates result according to the goal:lnorm,jacobi,ddm,wam")},

        {NULL,NULL,0,NULL}
};

static struct PyModuleDef moduledef={
    PyModuleDef_HEAD_INIT,
    "spkmeansmodule",
    NULL,
    -1,
    capiMethods
};

PyMODINIT_FUNC
PyInit_spkmeansmodule(void){
    PyObject *m;
    m=PyModule_Create(&moduledef);
    if(!m){
        return NULL;
    }
    return m;
}


