
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "spkmeans.h"
#include "../../../Users/edent/AppData/Local/Programs/Python/Python36/include/Python.h"

static Point *generalPy(PyObject *args);

static PyObject* createReturnedArray(double** matrixValues) {
    PyObject* returnedList=PyList_New(numberOfPoints);
    PyObject* returnedPoint;
    int i;
    int j;
    for (i = 0; i < numberOfPoints; i++) {
        returnedPoint=PyList_New(k);
        for (j = 0; j < k; j++) {
            PyList_SetItem(returnedPoint,j,PyFloat_FromDouble(matrixValues[i][j]));
        }
        PyList_SetItem(returnedList,i,returnedPoint);
    }
    return returnedList;
}

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

static PyObject * goalFunction(PyObject *self, PyObject *args){
    Point* points= generalPy(args);
    if(strcmp(goal,"spk")!=0){
        notSPK(points);
        Matrix T=createMatrix(numberOfPoints,k);
        return createReturnedArray(T.values);
    }
    else{
        Matrix T=spk(points);
        freePointPointer(points);
        return createReturnedArray(T.values);
    }
}


static void kmeansppC(PyObject *self, PyObject *args){
    Point *points1;
    Cluster *clusterArray1;
    PyObject* pyCentroids;
    PyObject* pyPoints;
    if(!PyArg_ParseTuple(args,"OOiiii",&pyPoints, &pyCentroids,&k,&d,&MAX_ITER,&numberOfPoints)){
        printf("An Error Has Occured\n");
        exit(0);
    }
    points1=getPointPointer();
    readAllPointsPy(points1,pyPoints);
    clusterArray1 = createClusterArrayWithCentroids(pyCentroids);
    kmeansWithInitialCentroids(points1, clusterArray1);
    freeClusterArray(clusterArray1);
    freePointPointer(points1);
}


static Point *generalPy(PyObject *args) {
    Point *points;
    char * filePath;
    PyArg_ParseTuple(args,"iss", &k,&filePath,&goal);
    FILE *file=fopen( filePath, "r" );
    countPoints(file);
    findD(file);
    points = getPointPointer();
    readAllPointsC(points,file);
    validateArguments();
    return points;
}


static PyMethodDef capiMethods[]={
        {"kmeanspp",
         (PyCFunction) kmeansppC,
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


