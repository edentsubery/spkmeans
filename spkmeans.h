#define PY_SSIZE_T_CLEAN


int d;
int k;
char *goal;
char *filePath;
int numberOfPoints;
int MAX_ITER;

typedef struct Matrix {
    double **values;
    int rows;
    int columns;
} Matrix;

typedef struct Eigenvalue {
    double value;
    int index;
} Eigenvalue;

typedef struct Point {
    double *coordinates;
    int index;
} Point;

typedef struct Cluster {
    Point sum_by_coordinates;
    Point centroid;
    int size;
} Cluster;

void freeMatrix(Matrix matrix);

void kmeansWithInitialCentroids(Point *points1, Cluster *clusterArray1);

int compare(const void *a, const void *b);

int countCommas(char *str);

Matrix spk(Point *points);

int EndsWithTail(char *fileName, char *tail);

void countPoints(FILE *file);

void findDTXT(FILE *file);

void readAllPointsCTXT(Point *points, FILE *file);

void readAllPointsCCSV(Point *points, FILE *file);

Point *getPointPointer(void);

Point *getPointsFromFile(void);

void freePointsArray(Point *points);

void calculateCentroids(Point *points, Cluster *clusterArray);

Point *csvCase(Point *points, FILE *file);

int calculateNewCentroid(Cluster *cluster);

Point *txtCase(Point *points, FILE *file);

double distanceFromCentroid(Cluster cluster, Point point);
void countPointsTXT(FILE *file);
void countPointsCSV(FILE *file);
Cluster *createClusterArray(Point *points);

void freeClusterArray(Cluster *cluster);

int closestCentroid(Cluster *clusters, Point point);

void addPointToArr(Point *pointsArr, Point point, int i);

void addPointToCluster(Cluster *cluster, Point point, int index);

Matrix createMatrix(int rows, int columns);

Matrix wam(Point *points);

Matrix ddm(Point *points);

Matrix sqrtDDM(Point *points);

Matrix lnorm(Point *points);

double sumForRow(Matrix matrix, int row, int length);

Matrix multiplyMatrixes(Matrix A, Matrix B);

double multiplyRowAndCol(Matrix A, Matrix B, int row, int col);

Matrix *findU(Matrix LNORM);

int *findPivot(Matrix A);

int convergence(Matrix A, Matrix AA);

Matrix fillLnorm(Point *points);

void printCentroids(const Cluster *clusterArray);

void printMatrix(Matrix mat);

Matrix transpose(Matrix p);

double offDiagonalSum(Matrix mat);

double calculateT(int i, int j, Matrix A);

Matrix createIdentityMatrix(void);

Matrix createPivotMatrix(Matrix A);

double *findCandS(double t);

double distance(Point *points, int i, int j);

void emptyClusters(Cluster *clusters);

Matrix columnsToRows(Matrix U);

void jacobi(Point *points);

void convertTMatrixToPoints(Point *points, Matrix T);

void kmeans(Matrix T);

Matrix createT(Eigenvalue *eigenvalues, Matrix U);

void normalizeRows(Matrix mat);

void normalizeRowForIndex(Matrix mat, int index);

void printEigenvalues(Eigenvalue *eigenvalues);

void eigengapHeuristic(Eigenvalue *eigenvalues);

int findGap(Eigenvalue *eigenvalues);

Eigenvalue *configureEigenvalues(Matrix U);

void sortEigenvalues(Eigenvalue *eigenvalues, Matrix U);

void validateArguments(void);

void analyzeArguments(char *const *argv);

Point *generalC(int args, char *argv[]);

void notSPK(Point *points);