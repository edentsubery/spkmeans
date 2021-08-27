
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "spkmeans.h"


struct Matrix{
    double** values;
    int rows;
    int columns;
};

struct Eigenvalue{
    double value;
    int index;
};

struct Point {
    double *coordinates;
    int index;
};

struct Cluster {
    Point sum_by_coordinates;
    Point centroid;
    int size;
};

int d;
int k;
char* goal;
int numberOfPoints;
int MAX_ITER ;


static void kmeansWithInitialCentroids(Point *points1, Cluster *clusterArray1);

static int countCommas(char *str) {
    int i;
    for (i = 0; str[i]; str[i] == ',' ? i++ : *str++);
    return i;
}

static void countPoints(FILE *file) {
    char line[1000];
    int numberOfLines;
    void* result;
    numberOfLines = 0;
    do {
        result=fgets(line,999,file);
        if (result==NULL) {
            break;
        }
        numberOfLines++;
    } while(1);
    numberOfPoints = numberOfLines;
}

static void findD(FILE *file) {
    int scanresult;
    char line[1000];
    rewind(file);
    scanresult=fscanf(file, "%s", line);
    if(scanresult!=0){
        d = countCommas(line) + 1;
    }
    
}

static void readAllPointsC(Point* points, FILE *file){
    Point *curPoint;
    int index;
    int i;
    double value;
    rewind(file);
    index = 0;
    while (index < numberOfPoints) {
        curPoint = malloc(sizeof(Point));
        assert(curPoint!=NULL);
        curPoint->coordinates = malloc(d * sizeof(double));
        assert(curPoint->coordinates!=NULL);
        for (i = 0; i < d; i++) {
            if (i == 0) {
                fscanf(file, "%lf", &value);
            } else {
                fscanf(file, ",%lf", &value);
            }
            curPoint->coordinates[i] = value;
        }
        addPointToArr(points, *curPoint, index++);

        if (index != numberOfPoints) {
            fscanf(file, "\n");
        }
    }
}

static Point* getPointPointer(void) {
    Point *points = malloc(numberOfPoints * sizeof(Point));
    assert(points!=NULL);
    return points;
}

static void calculateCentroids(Point *points, Cluster *clusterArray) {
    int j, i, index, changed;
    changed = 1;
    for (j = 0; changed && j < MAX_ITER; j++) {
        emptyClusters(clusterArray);
        for (i = 0; i < numberOfPoints; i++) {
            index = closestCentroid(clusterArray, points[i]);
            addPointToCluster(clusterArray, points[i],index);
        }
        changed = calculateNewCentroid(clusterArray);
    }
}

static void freePointPointer(Point *arr) {
    free(arr);
}

static int calculateNewCentroid(Cluster *cluster) {
    int changed, m, i;
    changed = 0;
    for (m = 0; m < k; m++) {
        for (i = 0; i < k; i++) {
            double coordinateAverage;
            if (cluster[m].size != 0) {
                coordinateAverage = cluster[m].sum_by_coordinates.coordinates[i] / (cluster[m].size);
            } else {
                coordinateAverage = 0;
            }
            if (cluster[m].centroid.coordinates[i] != coordinateAverage) {
                cluster[m].centroid.coordinates[i] = coordinateAverage;
                changed = 1;
            }
        }
    }
    return changed;
}

static double distanceFromCentroid(Cluster cluster, Point point) {
    double distance;
    int i;
    distance = 0;
    for (i = 0; i < k; i++) {
        distance += (point.coordinates[i] - cluster.centroid.coordinates[i])*(point.coordinates[i] - cluster.centroid.coordinates[i]);
    }
    return distance;
}

static Cluster *createClusterArray(Point* points) {
    int i,j;
    Cluster *arr;
    void *pointerForCentroid;
    void *pointerForSums;
    arr = malloc(k * sizeof(Cluster));
    assert(arr!=NULL);
    for (i = 0; i < k; i++) {
        pointerForCentroid = malloc(k * sizeof(double));
        assert(pointerForCentroid!=NULL);
        pointerForSums = malloc(k * sizeof(double));
        assert(pointerForSums!=NULL);
        arr[i].centroid.coordinates = pointerForCentroid;
        for(j=0;j<k;j++){
            arr[i].centroid.coordinates[j]=points[i].coordinates[j] ;
        }
        arr[i].sum_by_coordinates.coordinates = pointerForSums;
        arr[i].size = 0;
    }
    return arr;
}

static void freeClusterArray(Cluster *cluster) {
    int i;
    for (i = 0; i < k; i++) {
        free(cluster[i].centroid.coordinates);
        free(cluster[i].sum_by_coordinates.coordinates);
    }
    free(cluster);
}

static int closestCentroid(Cluster *clusters, Point point) {
    int index, i;
    double minDistance, dis;
    index = 0;
    minDistance = distanceFromCentroid(clusters[index], point);
    for (i = 1; i < k; i++) {
        dis = distanceFromCentroid(clusters[i], point);
        if (dis < minDistance) {
            minDistance = dis;
            index = i;
        }
    }
    point.index=index;
    return index;
}

static void addPointToArr(Point *pointsArr, Point point, int i) {
    pointsArr[i] = point;
}

static void addPointToCluster(Cluster *cluster, Point point,int index) {
    int j;
    for (j = 0; j < k; j++) {
        cluster[index].sum_by_coordinates.coordinates[j] += point.coordinates[j];
    }
    point.index=index;
    cluster[index].size++;
}

static Matrix createMatrix(int rows,int columns){
    struct Matrix matrix;
    double** values;
    int i;
    values=calloc(rows,sizeof(double*));
    assert(values!=NULL);
    for ( i = 0; i < rows; i++){
        values[i]=calloc(columns,sizeof(double));
        assert(values[i]!=NULL);
    }
    matrix.rows=rows;
    matrix.columns=columns;
    matrix.values=values;
    return matrix;
}

static Matrix wam(Point *points){
    Matrix matrix;
    int i,j;
    matrix=createMatrix(numberOfPoints,numberOfPoints);
    for ( i = 0; i < numberOfPoints; i++){
        for ( j = 0; j < numberOfPoints; j++){
            matrix.values[i][j]=exp(-distance(points,i,j)/2);
        }
    }
    for(i=0;i<numberOfPoints;i++){
        matrix.values[i][i]=0;
    }
    return matrix;
}

static Matrix ddm(Point *points){
    Matrix matrix,WAM;
    int i;
    WAM=wam(points);
    matrix=createMatrix(numberOfPoints,numberOfPoints);
    for(i=0;i<numberOfPoints;i++){
        matrix.values[i][i]=sumForRow(WAM,i,numberOfPoints);
    }
    return matrix;
}

static Matrix sqrtDDM(Point *points){
    Matrix matrix,DDM;
    int i;
    DDM=ddm(points);
    matrix=createMatrix(numberOfPoints,numberOfPoints);
    for(i=0;i<numberOfPoints;i++){
        matrix.values[i][i]=1/sqrt(DDM.values[i][i]);
    }
    return matrix;
}

static Matrix lnorm(Point *points){
    Matrix I,matrix,dwd,sqrtddm,WAM;
    int i,j;
    I=createIdentityMatrix();
    matrix=createMatrix(numberOfPoints,numberOfPoints);
    sqrtddm=sqrtDDM(points);
    WAM=wam(points);
    dwd=multiplyMatrixes(multiplyMatrixes(sqrtddm,WAM),sqrtddm);
    for(i=0;i<numberOfPoints;i++){
        for(j=0;j<numberOfPoints;j++){
            matrix.values[i][j]=I.values[i][j]-dwd.values[i][j];
        }
    }
    return matrix;
}

static double sumForRow(Matrix matrix,int row,int length){
    double sum;
    int i;
    sum=0;
    for(i=0;i<length;i++){
        sum+=matrix.values[row][i];
    }
    return sum;

}

static Matrix multiplyMatrixes(Matrix A,Matrix B){
    struct Matrix result;
    int i,j;
    result=createMatrix(A.rows,B.columns);
    for(i=0;i<A.rows;i++){
        for(j=0;j<B.columns;j++){
            result.values[i][j]=multiplyRowAndCol(A,B,i,j);
        }
    }
    return result;
}

static double multiplyRowAndCol(Matrix AA,Matrix BB,int row,int col){
    double sum = 0;
    int i;
    for(i=0;i<AA.columns;i++){
        sum+=AA.values[row][i]*BB.values[i][col];
    }
    return sum;
}


static Matrix* findU(Matrix LNORM){
    int i = 0;
    Matrix A,P,PT,AA,U;
    Matrix* UAA=malloc(2* sizeof(Matrix));
    A=LNORM;
    P=createPivotMatrix(A);
    PT=transpose(P);
    AA=multiplyMatrixes(multiplyMatrixes(PT,A),P);
    U=P;
    while((!convergence(A,AA))&&i<100){
        A=AA;
        P=createPivotMatrix(A);
        PT=transpose(P);
        AA=multiplyMatrixes(multiplyMatrixes(PT,A),P);
        U=multiplyMatrixes(U,P);
        i++;
    }
    UAA[0]=U;
    UAA[1]=AA;
    return UAA;
}

static int* findPivot(Matrix A){
    double max;
    int *pivot;
    int i,j;
    pivot=calloc(2, sizeof(int));
    assert(pivot!=NULL);
    max=0;
    for(i=0;i<numberOfPoints;i++){
        for(j=0;j<i;j++){
            if(i!=j){
                if(fabs(A.values[i][j])>max){
                    max=fabs(A.values[i][j]);
                    pivot[0]=i;
                    pivot[1]=j;
                }
            }
        }
    }
    return pivot;
}

static int convergence(Matrix A, Matrix AA){
    if(offDiagonalSum(A)-offDiagonalSum(AA)<=0.000000000000001){
        return 1;
    }
    else{
        return 0;
    }
}

static void printMatrix(Matrix mat) {
    int i;
    int j;
    for (i = 0; i < mat.rows; i++) {
        for (j = 0; j < mat.columns; j++) {
            if (j == mat.columns - 1) {
                printf("%0.4f\n", mat.values[i][j]);
            } else {
                printf("%0.4f,", mat.values[i][j]);
            }
        }
    }
}
static Matrix transpose(Matrix p){
    Matrix pt;
    int i,j;
    pt=createMatrix(p.rows,p.columns);
    for(i=0;i<p.rows;i++){
        for(j=0;j<p.columns;j++){
            pt.values[j][i]=p.values[i][j];
        }
    }
    return pt;
}

static double offDiagonalSum(Matrix mat){
    double sum;
    int i,j;
    sum=0;
    for(i=0;i<mat.rows;i++){
        for(j=0;j<mat.rows;j++){
            if(i!=j){
                sum+=pow(mat.values[i][j],2);
            }
        }
    }
    return sum;
}
static double calculateT(int i, int j,Matrix A){
    double theta,t;
    theta=(A.values[j][j]-A.values[i][i])/(2*(A.values[i][j]));
    if(theta<0){
        t=-1/(fabs(theta)+sqrt(pow(theta,2)+1));
    }
    else{
        t=1/(fabs(theta)+sqrt(pow(theta,2)+1));
    }
    return t;

}
static Matrix createIdentityMatrix(void){
    Matrix mat=createMatrix(numberOfPoints,numberOfPoints);
    int i;
    for(i=0;i<numberOfPoints;i++){
        mat.values[i][i]=1;
    }
    return mat;
}
static Matrix createPivotMatrix(Matrix A){
    int* pivot=findPivot(A);
    int i=pivot[0];
    int j=pivot[1];
    double t=calculateT(i,j,A);
    double* cs=findCandS(t);
    double c=cs[0];
    double s=cs[1];
    Matrix mat=createIdentityMatrix();
    mat.values[i][i]=c;
    mat.values[j][j]=c;
    mat.values[i][j]=s;
    mat.values[j][i]=-s;
    return mat;
}
static double* findCandS(double t){
    double *cs;
    cs=malloc(2* sizeof(double));
    assert(cs!=NULL);
    cs[0]=1/sqrt(pow(t,2)+1);
    cs[1]=cs[0]*t;
    return cs;
}

static double distance(Point *points,int i, int j){
    double distance;
    int l;
    distance=0;
    for(l=0;l<d;l++){
        distance+=pow(points[i].coordinates[l]-points[j].coordinates[l],2);
    }
    return sqrt(distance);
}

static void emptyClusters(Cluster *clusters) {
    int i, j;
    for (i = 0; i < k; i++) {
        for (j = 0; j < k; j++) {
            clusters[i].sum_by_coordinates.coordinates[j] = 0;
        }
        clusters[i].size = 0;
    }
}

static Matrix columnsToRows(Matrix U){
    Matrix mat;
    int i,j;
    mat=createMatrix(numberOfPoints,numberOfPoints);
    for(i=0;i<numberOfPoints;i++){
        for(j=0;j<numberOfPoints;j++){
            mat.values[i][j]=U.values[j][i];
        }
    }
    return mat;
}

static void jacobi(Point *points){
    Matrix* UAA=findU(fillLnorm(points));
    Matrix eigenvectors=columnsToRows(UAA[0]);
    printEigenvalues(configureEigenvalues(UAA[1]));
    printMatrix(eigenvectors);
}

static void convertTMatrixToPoints(Point* points,Matrix T){
    Point *curPoint;
    int index;
    int i;
    index = 0;
    while (index < numberOfPoints) {
        curPoint = malloc(sizeof(Point));
        assert(curPoint!=NULL);
        curPoint->coordinates = malloc(k * sizeof(double));
        assert(curPoint->coordinates!=NULL);
        for (i = 0; i < k; i++) {
            curPoint->coordinates[i] = T.values[index][i];
        }
        addPointToArr(points, *curPoint, index++);
    }
}

static Matrix fillLnorm(Point *points){
    Matrix mat;
    int i,j;
    mat=createMatrix(numberOfPoints,numberOfPoints);
    for(i=0;i<numberOfPoints;i++){
        for(j=0;j<numberOfPoints;j++){
            mat.values[i][j]=points[i].coordinates[j];
        }
    }
    return mat;
}

static void kmeans(Matrix T){
    Point *points1;
    Cluster *clusterArray1;
    d=k;
    points1=malloc(T.rows*sizeof(struct Point));
    convertTMatrixToPoints(points1,T);
    clusterArray1 = createClusterArray(points1);
    kmeansWithInitialCentroids(points1, clusterArray1);
    freeClusterArray(clusterArray1);
    freePointPointer(points1);
}

static void kmeansWithInitialCentroids(Point *points1, Cluster *clusterArray1) {
    calculateCentroids(points1, clusterArray1);
    printCentroids(clusterArray1);
}

static void printCentroids(const Cluster *clusterArray) {
    int i;
    int j;
    for (i = 0; i < k; i++) {
        for (j = 0; j < d; j++) {
            if (j == d - 1) {
                printf("%0.4f\n", clusterArray[i].centroid.coordinates[j]);
            } else {
                printf("%0.4f,", clusterArray[i].centroid.coordinates[j]);
            }
        }
    }
}

static void printEigenvalues(Eigenvalue* eigenvalues){
    int j;
    for (j = 0; j < numberOfPoints; j++) {
        if (j == numberOfPoints- 1) {
            printf("%0.4f\n", eigenvalues[j].value);
        } else {
            printf("%0.4f,", eigenvalues[j].value);
        }
    }

}

int main(int args, char *argv[]){
    Point *points;
    FILE *file=fopen( argv[3], "r" );
    assert(args==4);
    analyzeArguments(argv);
    countPoints(file);
    findD(file);
    points = getPointPointer();
    readAllPointsC(points,file);
    validateArguments();
    if(strcmp(goal,"wam") == 0){
        printMatrix(wam(points));
    }
    else if(strcmp(goal,"ddg") == 0){
        printMatrix(ddm(points));
    }
    else if(strcmp(goal,"lnorm") == 0){
        printMatrix(lnorm(points));
    }  
    else if(strcmp(goal,"jacobi") == 0){
        jacobi(points);
    }
    else{
        kmeans(spk(points));

    }
    freePointPointer(points);
    return 0;
}

static Matrix spk(Point *points) {
    Matrix LNORM=lnorm(points);
    Matrix* UAA=findU(LNORM);
    Eigenvalue * eigenvalues=configureEigenvalues(UAA[1]);
    sortEigenvalues(eigenvalues,UAA[0]);
    if(k==0){
        printEigenvalues(eigenvalues);
        eigengapHeuristic(eigenvalues);
    }
    return createT(eigenvalues,UAA[0]);
}

static Matrix createT(Eigenvalue* eigenvalues,Matrix U){
    Matrix mat;
    int i,j,col;
    mat=createMatrix(numberOfPoints,k);
    for(i=0;i<k;i++){
        col=eigenvalues[i].index;
        for(j=0;j<numberOfPoints;j++){
            mat.values[j][i]=U.values[j][col];
        }
    }
    normalizeRows(mat);
    return mat;
}

static void normalizeRows(Matrix mat){
    int i;
    for(i=0;i<numberOfPoints;i++){
        normalizeRowForIndex(mat,i);
    }
}
static void normalizeRowForIndex(Matrix mat, int index){
    double sum,size;
    int i;
    sum=0;
    for (i=0;i<k;i++){
        sum+=pow(mat.values[index][i],2);
    }
    size=sqrt(sum);
    for (i=0;i<k;i++){
        mat.values[index][i]=mat.values[index][i]/size;
    }
}
static void sortEigenvalues (Eigenvalue* eigenvalues,Matrix U){
    qsort (eigenvalues,U.columns, sizeof(Eigenvalue), (int (*)(const void *, const void *)) compare);
}

static Eigenvalue * configureEigenvalues(Matrix AA){
    Eigenvalue * eg;
    int i;
    eg=malloc(numberOfPoints* sizeof(Eigenvalue));
    assert(eg!=NULL);
    for(i=0;i<numberOfPoints; i++) {
        eg[i].index = i;
        eg[i].value = AA.values[i][i];
    }
    return eg;
}

static int compare (const void * a, const void * b)
{
    double diff=(*(Eigenvalue *)a).value - (*(Eigenvalue *)b).value;
    if(diff==0){
        return (*(Eigenvalue *)a).index - (*(Eigenvalue *)b).index;
    }
    if(diff>0){
        return 1;
    }
    else{
        return -1;
    }
}

static void eigengapHeuristic(Eigenvalue* eigenvalues){
    k= findGap(eigenvalues);
}

static int findGap(Eigenvalue* eigenvalues){
    double max;
    int i,index;
    index=0;
    max=0;
    for(i=1;i<numberOfPoints/2;i++){
        if(fabs(eigenvalues[i].value-eigenvalues[i-1].value)>max){
            max=fabs(eigenvalues[i].value-eigenvalues[i-1].value);
            index=i;
        }
    }
    return index;
}


static void analyzeArguments(char *const *argv) {
    k = atoi(argv[1]);
    if (k < 0) {
        printf("Invalid value for number of means.\n");
    }
    goal=argv[2];
}

static void validateArguments(void){
    int valid=1;
    if (d<=0){
        printf("Invalid value for d");
        valid=0;
    }
    if(k<0){
        printf("Invalid value for k");
        valid=0;
    }
    if(numberOfPoints<k){
        printf("Invalid value for K: bigger than amount of points");
        valid=0;
    }
    if(valid==0){
        exit(0);
    }
}

