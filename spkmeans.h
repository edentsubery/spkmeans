#define PY_SSIZE_T_CLEAN

typedef struct Cluster Cluster;
typedef struct Point Point;
typedef struct Matrix Matrix;
typedef struct Eigenvalue Eigenvalue;


static int  compare (const void * a, const void * b);
static int countCommas(char *str);
static void spk(Point *points);
static void countPoints(FILE *file);
static void findD(FILE *file);
static void readAllPointsC(Point* points, FILE *file);

static Point* getPointPointer(void);

static void calculateCentroids(Point *points, Cluster *clusterArray);

static void freePointPointer(Point *arr);

static int calculateNewCentroid(Cluster *cluster);

static double distanceFromCentroid(Cluster cluster, Point point);

static Cluster *createClusterArray(Point* points);

static void freeClusterArray(Cluster *cluster);

static int closestCentroid(Cluster *clusters, Point point);

static void addPointToArr(Point *pointsArr, Point point, int i);

static void addPointToCluster(Cluster *cluster, Point point,int index);

static Matrix createMatrix(int rows,int columns);
static Cluster * createClusterArrayEmpty(void);
static void wam(Point *points);

static void ddm(Point *points);

static void sqrtDDM(Point *points);
static void lnorm(Point *points);

static double sumForRow(Matrix matrix,int row,int length);

static Matrix multiplyMatrixes(Matrix A,Matrix B);

static double multiplyRowAndCol(Matrix A,Matrix B,int row,int col);

static void findU(void);

static int* findPivot(void);
static int convergence(void);
static void fillLnorm(Point *points);
static void printCentroids(const Cluster *clusterArray);
static void printMatrix(Matrix mat);
static Matrix transpose(Matrix p);
static double offDiagonalSum(Matrix mat);
static double calculateT(int i, int j);
static Matrix createIdentityMatrix(void);
static Matrix createPivotMatrix(void);
static double* findCandS(double t);

static double distance(Point *points,int i, int j);
static void emptyClusters(Cluster *clusters) ;

static void columnsToRows(void);

static void jacobi(Point *points);

static void convertTMatrixToPoints(Point* points);

static void kmeans(void);

static void evaluateClusters(Point* points1,Cluster* actualClusters);

static void createT(void);

static void normalizeRows(Matrix mat);
static void normalizeRowForIndex(Matrix mat, int index);

static void printEigenvalues(void);
static void eigengapHeuristic(void);

static int findGap(void);
static void configureEigenvalues(void);
static void sortEigenvalues (void);
static void validateArguments(void);
static void analyzeArguments( char *const *argv);