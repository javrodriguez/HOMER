
// Copyright 2009, 2010, 2011, 2012 Christopher Benner <cbenner@gmail.com>
// 
// This file is part of HOMER
//
// HOMER is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// HOMER is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.




#include <stdio.h>
#include "Hashtable.h"
#include "statistics.h"
#include <limits.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>

#ifndef CLUSTERING_H
#define CLUSTERING_H

#define WHITE_SPACE 7
#define BUFFER_CLUSTER 1000000
#define CLUSTER_DISTANCE_CORRELATION 1
#define CLUSTER_DISTANCE_EUCLIDEAN 2
#define CLUSTER_SUB_INIT 7000

class TreeCluster;
class TreeNode;

class TreeCluster {
public:

	double** matrix;
	int numMatrix;
	char** names;

	double** ogMatrix;
	int ogNumMatrix;
	char** ogNames;

	double** dataMatrix;
	double** ogDataMatrix;
	int ogNumDataMatrix;
	char** ogDataMatrixNames;
	char** expNames;
	int* sub2matrix;
	int* matrix2sub;
	int numExps;
	int numSub;
	int subFlag;
	int dataCDTflag;

	TreeNode* clusters;
	int* mask;
	int rootTree;

	int fixedClusterScore;
	int fixedPositionFlag;
	int reverseFlag;
	int distanceMetric;

	TreeCluster();
	TreeCluster(double ** m, int numMatrix);
	~TreeCluster();

	double* minDistCache;
	int* minIndexCache;

	double** clusterSummary;
	int numClusters;
	
	void cluster();
	int* getOrder();
	int* getClusters(double threshold, int &numClusters);
	void init();
	void init(double** m, int numMatrix);
	void readDistMatrix(char* filename, int numAnnCols);
	void readDataMatrix(char* filename, int numAnnCols);
	void sampleDataMatrix(int sub);
	void invertDistMatrix();
	void setReverseFlag();
	void initializeMatrix();
	void calcDistMatrix(int dm);
	void printCDT(char* prefix,int matrixFlag);
	void printClusters(FILE* fp,char**names,double thresh);
	double findMinimumPair(int &left, int &right);
	double initializeCache(int &left,int &right);
	double updateCache(int &left,int &right);
	void printClusterSizes(FILE* fp,int maxClusterSizeToReport);
	void scoreClusterFiles(char** files, int numFiles, FILE* fp);
	void loadGTR(char* gtrFile);
	void expandClusters();
};

void split(char* string, char** cols, int &numCols, char delim);


class TreeNode {
public:
	int left;
	int right;
	double distance;
	double weight;
	
	int matrixIndex;

	TreeNode();
	~TreeNode();
	void print(FILE*);
};



#endif
