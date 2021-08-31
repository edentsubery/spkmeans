import numpy as np
import sys
import spkmeansmodule as sp


def convertToList(array):
    n = len(array)
    d = len(array[0])
    arrayList = []
    for i in range(n):
        arrayRow = []
        for j in range(d):
            arrayRow.append(array[i][j])
        arrayList.append(arrayRow)
    return arrayList


def utility(points, centroids, points_min, distribution, i):
    for index, point in enumerate(points):
        distribution[index] = distance(point, centroids[i], points_min[index])
        points_min[index] = distribution[index]
    sum = np.sum(distribution)
    distribution = distribution / sum
    return distribution


def distance(point, centroid, min):
    temp = np.sum(np.power(centroid - point, 2))
    if temp < min:
        min = temp
    return min


def kmeans_pp(result):
    n = len(result)
    k = len(result[0])
    np.random.seed(0)
    indexes = np.empty(k)
    centroids = np.empty((k, k))
    points_min = np.full(n, np.inf)
    calculate_initial_centroids(centroids, indexes, k, n, result, points_min)
    print(convertToList(centroids))
    print_indexes(indexes, k)
    sp.kmeanspp(convertToList(result), convertToList(centroids), k, k, 300, n)


def analyze_arguments():
    k = int(sys.argv[1])
    path = str(sys.argv[3])
    goal = str(sys.argv[2])
    return k, path, goal


def validate_arguments(k, d, n):
    if k >= n:
        print("Invalid Input!")
        exit(0)
    if d <= 0:
        print("Invalid Input!")
        exit(0)


def calculate_initial_centroids(centroids, indexes, k, n, points, points_min):
    index = np.random.choice(n)
    indexes[0] = index
    centroids[0] = points[index]
    distribution = np.empty(n)
    for i in range(1, k):
        distribution = utility(points, centroids, points_min, distribution, i - 1)
        index = np.random.choice(n, 1, p=distribution)
        indexes[i] = index[0]
        centroids[i] = points[index[0]]


def print_indexes(indexes, k):
    for i in range(k - 1):
        print(str(int(indexes[i])) + ",", end="")
    print(str(int(indexes[k - 1])))


def spkmeans():
    k, path, goal = analyze_arguments()
    result = sp.goal(k, path, goal)
    if goal == "spk":
        kmeans_pp(result)


spkmeans()
