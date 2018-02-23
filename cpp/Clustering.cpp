
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



#include "Clustering.h"

void clusterSplit(char* string, char** cols, int &numCols, char delim);
TreeCluster::TreeCluster(double** distMatrix, int matrixLength) {
	init();
	init(distMatrix, matrixLength);
}
TreeCluster::TreeCluster() {
	init();
}
void TreeCluster::init() {
	matrix = NULL;
	numMatrix = 0;
	clusters = NULL;
	mask = NULL;
	fixedClusterScore = 0;
	fixedPositionFlag = 0;
	reverseFlag = 0;
	minDistCache = NULL;
	minIndexCache = NULL;
	names = NULL;
	distanceMetric = CLUSTER_DISTANCE_CORRELATION;
    dataMatrix= NULL;
    expNames=NULL;
	numExps=0;
	numSub = CLUSTER_SUB_INIT;
	ogNumDataMatrix = 0;
	ogDataMatrix = NULL;
	ogDataMatrixNames = NULL;
	ogMatrix = NULL;
	subFlag = 0;
	sub2matrix = NULL;
	matrix2sub = NULL;
	dataCDTflag = 0;	
	clusterSummary = NULL;

}
void TreeCluster::readDistMatrix(char* file, int numAnnCols) {
	numMatrix = 0;
	matrix = NULL;
	names = NULL;

	char* buf = new char[BUFFER_CLUSTER];
	char** line = new char*[10000];
	int numCols = 0;

	FILE* fp = fopen(file,"r");
	if (fp == NULL) {
		fprintf(stderr, "!! Could not open input file: %s\n",file);
		exit(0);
	}

	int lineCount = 0;
	while (fgets(buf, BUFFER_CLUSTER, fp) != NULL) {
		split(buf, line, numCols, '\t');
		lineCount++;
		if (lineCount == 1) {
			numMatrix = numCols-numAnnCols;
			matrix = new double*[numMatrix];
			names = new char*[numMatrix];
			for (int i=0;i<numMatrix;i++) {
				matrix[i] = new double[numMatrix];
				for (int j=0;j<numMatrix;j++) matrix[i][j]=FLT_MAX;
			}
			continue;
		}
		if (numCols < numMatrix+numAnnCols) {
			fprintf(stderr, "Something is wrong with your distance matrix!!!\n");
			exit(0);
		}
		int index = lineCount -2;
		names[index] = new char[strlen(line[0])+1];
		strcpy(names[index], line[0]);
		for (int i=numAnnCols;i<numCols;i++) {
			sscanf(line[i], "%lf", &(matrix[index][i-numAnnCols]));
		}
	}
	fclose(fp);

	delete []buf;
	delete []line;
	initializeMatrix();

}

void TreeCluster::readDataMatrix(char* file, int numAnnCols) {

	char* buf = new char[BUFFER_CLUSTER];
	char** line = new char*[10000];
	int numCols = 0;

	FILE* fp = fopen(file,"r");
	if (fp == NULL) {
		fprintf(stderr, "!! Could not open input file: %s\n",file);
		exit(0);
	}
	int lineCount = 0;
	while (fgets(buf, BUFFER_CLUSTER, fp) != NULL) {
		lineCount++;
	}
	rewind(fp);

	numMatrix = lineCount-1;
	dataMatrix = new double*[numMatrix];
	names = new char*[numMatrix];

	lineCount = 0;
	while (fgets(buf, BUFFER_CLUSTER, fp) != NULL) {
		split(buf, line, numCols, '\t');
		lineCount++;
		if (lineCount == 1) {
			numExps = numCols-numAnnCols;
			expNames = new char*[numExps];
			for (int i=numAnnCols;i<numCols;i++) {
				expNames[i-numAnnCols] = new char[strlen(line[i])+1];
				strcpy(expNames[i-numAnnCols],line[i]);
			}
			continue;
		}
		if (numCols < numExps+numAnnCols) {
			fprintf(stderr, "Something is wrong with your raw data matrix file!!!\n");
			exit(0);
		}
		int index = lineCount -2;
		names[index] = new char[strlen(line[0])+1];
		strcpy(names[index], line[0]);

		dataMatrix[index]=new double[numExps];
		for (int i=numAnnCols;i<numCols;i++) {
			sscanf(line[i], "%lf", &(dataMatrix[index][i-numAnnCols]));
		}
	}	
	fclose(fp);

	fprintf(stderr, "\tFinished reading data (%d rows by %d columns)\n", numMatrix, numExps);

	delete []buf;
	delete []line;
}
void TreeCluster::sampleDataMatrix(int sub) {


	if (sub >= (int)(numMatrix*0.9)) {
		fprintf(stderr, "!!! warning sampling level is > 90%% of total data points (skipping...)\n");
		return;
	}

	numSub = sub;
	subFlag = 1;
	
	ogNumDataMatrix = numMatrix;
	ogDataMatrix = dataMatrix;
	ogDataMatrixNames = names;



	srand(time(NULL));
	dataMatrix = new double*[numMatrix];
	names = new char*[numMatrix];
	sub2matrix=new int[numMatrix];
	matrix2sub=new int[ogNumDataMatrix];
	char* randMask = new char[ogNumDataMatrix];

	for (int i=0;i<ogNumDataMatrix;i++) {
		matrix2sub[i] = -1;
		sub2matrix[i] = -1;
		//sub2matrix[i] = 0;
		randMask[0] = 0;
	}

	int newNumMatrix=0;
	while (newNumMatrix < numSub) {
		int index = rand() % numMatrix;
		if (randMask[index]) continue; //already choosen, try again
		dataMatrix[newNumMatrix] = ogDataMatrix[index];
		names[newNumMatrix] = ogDataMatrixNames[index];
		sub2matrix[newNumMatrix] = index;
		matrix2sub[index] = newNumMatrix;
		randMask[index] = 1;
		newNumMatrix++;
	}
	numMatrix = newNumMatrix;
	delete []randMask;
}

void TreeCluster::calcDistMatrix(int distmetric) {

	distanceMetric = distmetric;
	fprintf(stderr, "\tCalculing distance matrix for %d data points\n", numMatrix);
	matrix = new double*[numMatrix];
	for (int i=0;i<numMatrix;i++) {
		matrix[i] = new double[numMatrix];
		for (int j=0;j<numMatrix;j++) {
			matrix[i][j]=FLT_MAX;
		}
	}
	if (distanceMetric == CLUSTER_DISTANCE_CORRELATION) {
		reverseFlag=1;
	}
	
	for (int i=0;i<numMatrix;i++) {
		for (int j=i;j<numMatrix;j++) {
			double score = 0.0;
			if (distanceMetric == CLUSTER_DISTANCE_CORRELATION) {
				if (i==j) {
					matrix[i][j]=1.0;
					continue;
				}
				score = correlation(dataMatrix[i], dataMatrix[j], numExps);
				matrix[i][j]=score;
				matrix[j][i]=score;
			} else if (distanceMetric == CLUSTER_DISTANCE_EUCLIDEAN) {
				if (i==j) {
					matrix[i][j]=0.0;
					continue;
				}
				score = euclidean(dataMatrix[i], dataMatrix[j], numExps);
				matrix[i][j]=score;
				matrix[j][i]=score;
			}
		}
	}
	if (reverseFlag) {
		invertDistMatrix();
	}
	initializeMatrix();
}
void TreeCluster::invertDistMatrix() {
	for (int i=0;i<numMatrix;i++) {   
		for (int j=0;j<numMatrix;j++) {
			matrix[i][j] *= -1;
		}
	}
}
void TreeCluster::init(double** distMatrix, int matrixLength) {

	init();
	fixedPositionFlag = 0;
	reverseFlag = 0;
	fixedClusterScore = 0;

	numMatrix = matrixLength;
	matrix = new double*[numMatrix];
	for (int i=0;i<numMatrix;i++) {
		matrix[i] = new double[numMatrix];
		for (int j=0;j<numMatrix;j++) {
			matrix[i][j] = distMatrix[i][j];
		}
	}
	initializeMatrix();
}

void TreeCluster::expandClusters() {


	rootTree = numMatrix-2;
	int* clusterLocation = new int[ogNumDataMatrix];
	for (int i=0;i<ogNumDataMatrix;i++) clusterLocation[i]=-1;
	TreeNode* newClusters = new TreeNode[ogNumDataMatrix-1];
	double minDistance = FLT_MAX;

	//figure out where leaves are, and update their indexes to the full matrix
	for (int i=0;i<numMatrix-1;i++) {
		newClusters[i].left = clusters[i].left;
		newClusters[i].right = clusters[i].right;
		newClusters[i].distance = clusters[i].distance;
		newClusters[i].weight = clusters[i].weight;
		newClusters[i].matrixIndex = clusters[i].matrixIndex;
		if (clusters[i].distance < minDistance) minDistance = clusters[i].distance;
	
		if (newClusters[i].left > -1) {
			int leafIndex = sub2matrix[newClusters[i].left];
			newClusters[i].left = leafIndex;
			clusterLocation[leafIndex] = i;
		}
		if (newClusters[i].right > -1) {
			int leafIndex = sub2matrix[newClusters[i].right];
			newClusters[i].right = leafIndex;
			clusterLocation[leafIndex] = i;
		}
	}
	int curClusters = numMatrix-1;

	fprintf(stderr, "\tCluster expansion progress:\n");
	fprintf(stderr, "\t|0%%                                   50%%                                  100%%|\n\t");
	double nextMark = 0.0;
	double step = 1.0/80.0;

	for (int i=0;i<ogNumDataMatrix;i++) {
		while (((double)i)/((double)ogNumDataMatrix) >= nextMark) {
			nextMark+=step;
			fprintf(stderr, "=");
		}
		if (matrix2sub[i] < 0) {
			double bestScore = FLT_MAX;
			int bestIndex = 0;
			for (int j=0;j<numMatrix;j++) {
				double score = 0.0;
				if (distanceMetric == CLUSTER_DISTANCE_CORRELATION) {
					score = correlation(ogDataMatrix[i], dataMatrix[j], numExps);
				} else if (distanceMetric == CLUSTER_DISTANCE_EUCLIDEAN) {
					score = euclidean(ogDataMatrix[i], dataMatrix[j], numExps);
				}
				if (reverseFlag) {
					score *= -1;
				}
				if (score < bestScore) {
					bestScore = score;
					bestIndex = sub2matrix[j];
				}
			}
			//Create new cluster to incorporate leaf
			//first get cluster where the leaf is located
			int nodeIndex = clusterLocation[bestIndex];
			newClusters[curClusters].left = bestIndex;
			newClusters[curClusters].right = i;
			if (newClusters[nodeIndex].left == bestIndex) {
				newClusters[nodeIndex].left = -1*(curClusters+1);
			} else if (newClusters[nodeIndex].right == bestIndex) {
				newClusters[nodeIndex].right = -1*(curClusters+1);
			} else {
				fprintf(stderr, "!!!! Something is wrong....\n");
			}
			newClusters[curClusters].distance = minDistance;
			newClusters[curClusters].weight = newClusters[nodeIndex].weight;

			//update location of leaves in the cluster
			clusterLocation[bestIndex] = curClusters;
			clusterLocation[i] = curClusters;
			
			curClusters++;
		}
	}
	fprintf(stderr, "\n");

	delete []clusters;
	clusters = newClusters;

	delete []dataMatrix;	
	dataMatrix = ogDataMatrix;
	numMatrix = ogNumDataMatrix;
	names = ogDataMatrixNames;

//fprintf(stderr, "Expand Tree:\n");
//for (int i=0;i<numMatrix-1;i++) clusters[i].print(stderr);

}

void TreeCluster::initializeMatrix() {
	fprintf(stderr, "\tDistance matrix size = %d\n", numMatrix); 

	//first make copy (it will be modified during the clustering)
	ogNumMatrix = numMatrix;
	ogMatrix = new double*[numMatrix];
	ogNames = new char*[numMatrix];
	for (int i=0;i<numMatrix;i++) {
		ogMatrix[i] = new double[numMatrix];
		ogNames[i] = new char[strlen(names[i])+1];
		strcpy(ogNames[i], names[i]);
		for (int j=0;j<numMatrix;j++) {
			ogMatrix[i][j] = matrix[i][j];
		}
	}


	clusters = new TreeNode[numMatrix-1];

	mask = new int[numMatrix*2-1];
	for (int i=0;i<numMatrix*2-1;i++) {
		mask[i] = 0;
	}
	minDistCache = new double[numMatrix];
	minIndexCache = new int[numMatrix];
	for (int i=0;i<numMatrix;i++) {
		minDistCache[i] = 1e10;
		minIndexCache[i] = 0;
	}
}
void TreeCluster::setReverseFlag() {
	reverseFlag=1;
	for (int i=0;i<numMatrix;i++) {
		for (int j=0;j<numMatrix;j++) {
			matrix[i][j] = -1*matrix[i][j];
		}
	}
}
	

TreeCluster::~TreeCluster() {
	if (clusters != NULL) delete []clusters;
	if (matrix != NULL) {
		for (int i=0;i<numMatrix;i++){ 
			delete [](matrix[i]);
		}
		delete []matrix;
	}
	if (dataMatrix != NULL) {
		for (int i=0;i<numMatrix;i++){ 
			delete [](dataMatrix[i]);
		}
		delete []dataMatrix;
	}
	if (mask != NULL) delete []mask;
	if (names != NULL) {
		for (int i=0;i<numMatrix;i++) {
			if (names[i] != NULL) delete [](names[i]);
		}
		delete []names;
	}
	if (expNames != NULL) {
		for (int i=0;i<numExps;i++) {
			if (expNames[i] != NULL) delete [](expNames[i]);
		}
		delete []expNames;
	}
	if (sub2matrix != NULL) delete []sub2matrix;
	if (matrix2sub != NULL) delete []matrix2sub;
	if (minDistCache != NULL) delete []minDistCache;
	if (minIndexCache != NULL) delete []minIndexCache;
}

void TreeCluster::loadGTR(char* gtrFile) {

	int bufSize = 1000000;
	char* buf = new char[bufSize];
	char** line = new char*[100000];
	int numCols = 0;


	char* cdtFile = new char[strlen(gtrFile)+1];
	strcpy(cdtFile, gtrFile);
	strcpy(&(cdtFile[strlen(cdtFile)-3]),"cdt");

	int lineCount = 0;
	FILE* fp = fopen(gtrFile, "r");
	while (fgets(buf,bufSize,fp) != NULL) lineCount++;
	fclose(fp);
	
	numMatrix = lineCount+1;
	fprintf(stderr, "\t\tMatrix Size: %d\n", numMatrix);
	matrix = NULL;
	names = new char*[numMatrix];

	dataCDTflag = 0;	

	Hashtable* nameHash = new Hashtable();
	Inttable* geneIndex = new Inttable();
	fp = fopen(cdtFile, "r");
	if (fp != NULL) {
		fprintf(stderr, "\tLoading cluster leaf names from cdt file\n");
		int lineCount = 0;
		while (fgets(buf,bufSize,fp) != NULL) {
			lineCount++;
			clusterSplit(buf, line, numCols, '\t');
			//if (numCols < 2) continue;
			if (lineCount == 1) {
				//numMatrix = numCols-4;
				//fprintf(stderr, "\t\tMatrix Size: %d\n", numMatrix);
				if (numCols == numMatrix+4) {
					dataCDTflag = 0;	
					matrix= new double*[numMatrix];
					for (int i=0;i<numMatrix;i++){ 
						matrix[i] = new double[numMatrix];
						for (int j=0;j<numMatrix;j++) matrix[i][j] = 0.0;
					}
				} else {
					dataCDTflag = 1;
					numExps = numCols-4;
					dataMatrix = new double*[numMatrix];
					for (int i=0;i<numMatrix;i++){ 
						dataMatrix[i] = new double[numExps];
						for (int j=0;j<numExps;j++) dataMatrix[i][j] = 0.0;
					}
				}
			}
			if (lineCount < 3) continue;
			char* n = new char[strlen(line[1])+1];
			strcpy(n,line[1]);
			nameHash->insert(n,line[0]);
			names[lineCount-3] = n;
			geneIndex->insert(lineCount-3,line[0]);
			for (int i=4;i<numCols;i++) {
				if (dataCDTflag) {
					sscanf(line[i],"%lf",&(dataMatrix[lineCount-3][i-4]));
				} else {
					sscanf(line[i],"%lf",&(matrix[lineCount-3][i-4]));
				}
			}
		}
		fclose(fp);
		//int altNumMatrix = lineCount-2;
		//fprintf(stderr, "numMatrix = %d (gtr) vs. %d (cdt)\n", numMatrix, altNumMatrix);
		fprintf(stderr, "\tNumber of data records = %d\n", numMatrix);
	} else {
		fprintf(stderr, "\tCould not find cdt file - exiting\n");
		exit(0);
	}
	delete []cdtFile;

	/*
	char** keys = nameHash->keys();
	for (int i=0;i<nameHash->total;i++) {
		char* n = (char*)nameHash->search(keys[i]);
		//fprintf(stderr, "%d\t%s\n", i,n);
		geneIndex->insert(i,keys[i]);
		delete [](keys[i]);
	}
	delete []keys;
	*/
	delete nameHash;
		

	clusters = new TreeNode[numMatrix-1];
	mask = new int[numMatrix*2-1];
	for (int i=0;i<numMatrix*2-1;i++) {
		mask[i] = 0;
	}

	fp = fopen(gtrFile, "r");
	if (fp == NULL) {
		fprintf(stderr, "!!!! Problem openning gtr file (%s)\n", gtrFile);
		exit(0);
	}
	fprintf(stderr, "\tLoading Tree from gtr file\n");
	Inttable* nodeIndex = new Inttable();
	lineCount = -1;
	double d = 0;
	while (fgets(buf,bufSize,fp) != NULL) {
		lineCount++;
		clusterSplit(buf, line, numCols, '\t');
		if (numCols < 4) continue;
		sscanf(line[3],"%lf",&d);
		clusters[lineCount].distance = -1*d;

		nodeIndex->insert(lineCount,line[0]);

		double w = 0.0;
		int left = geneIndex->search(line[1]);
		if (left == EMPTY_INT) {
			left = nodeIndex->search(line[1]);
			if (left == EMPTY_INT) {
				fprintf(stderr, "!!! Error parsing left node: %s (line:%d)\n", line[1],lineCount+1);
				exit(0);
			}
			clusters[lineCount].left = -1*(left+1);
			w += clusters[left].weight;
		} else {
			clusters[lineCount].left = left;
			w += 1.0;
		}
		int right = geneIndex->search(line[2]);
		if (right == EMPTY_INT) {
			right = nodeIndex->search(line[2]);
			if (right == EMPTY_INT) {
				fprintf(stderr, "!!! Error parsing right node: %s (line:%d)\n", line[2],lineCount+1);
				exit(0);
			}
			clusters[lineCount].right = -1*(right+1);
			w += clusters[right].weight;
		} else {
			clusters[lineCount].right = right;
			w += 1.0;
		}
	}
	fclose(fp);

	rootTree = numMatrix-2;


	delete []buf;
	delete []line;
	delete geneIndex;
	delete nodeIndex;
}


void TreeCluster::printClusters(FILE* fp, char** geneNames,double thresh) {
	//char* filename = new char[10000];
	//sprintf(filename,"%s.clusters.txt",prefix);
	//FILE* fp = fopen(filename,"w");

	if (geneNames == NULL) {
		geneNames = names;
	}

	numClusters = 0;
	int* assignment = getClusters(thresh, numClusters);

	if (dataCDTflag) {
		clusterSummary = new double*[numClusters];
		for (int i=0;i<numClusters;i++) {
			clusterSummary[i] = new double[numExps+1];
			for (int j=0;j<numExps+1;j++) {
				clusterSummary[i][j]=0.0;
			}
		}
	}

	for (int i=0;i<numMatrix;i++) {
		//int cluster = assignment[i]-1;
		fprintf(fp, "%s\t%d\n",geneNames[i],assignment[i]);
	}
	/*
		if (dataCDTflag) {
			clusterSummary[cluster][numExps]+=1.0;
			for (int j=0;j<numExps;j++) {
				clusterSummary[cluster][i]+=dataMatrix[i][j];
			}
		}
	}
	if (dataCDTflag) {
		for (int i=0;i<numClusters;i++) {
			if (clusterSummary[i][numExps] < 1) continue;
			for (int j=0;j<numExps;j++) {
				clusterSummary[i][j] /= clusterSummary[i][numExps];
			}
		}
	}*/
}

int* TreeCluster::getClusters(double thresh, int &numClusters) {

	int* assignment = new int[numMatrix];

	//int* order = new int[numMatrix];	
	int* lastTree = new int[numMatrix];	
	int* lastIndex = new int[numMatrix];	
	for (int i=0;i<numMatrix;i++) {
		//order[i] = 0;
		lastTree[i] = 0;
		lastIndex[i] = 0;
	}
	//rootTree = numMatrix-2;
	int curTree = rootTree;

	int clusterNumber = 1;
	int currentCluster = -1;

	while (!(curTree == rootTree && lastIndex[rootTree]==2)) {
		double score = clusters[curTree].distance;
		//score = (double)curTree/(double)numMatrix;
		//fprintf(stderr, "%lf\n", score);
		if (score < thresh) {
			if (currentCluster < 0) {
				currentCluster = clusterNumber++;
			}
		} else {
			currentCluster = -1;
		}
				
		if (lastIndex[curTree] == 0) {
			int i = clusters[curTree].left;
			if (i >= 0) {
				//order[curOrder++] = i;
				assignment[i] = currentCluster;
				lastIndex[curTree]++;
				continue;
			} else {
				lastIndex[curTree]++;
				int treeIndex = (-1*i)-1;
				lastTree[treeIndex] = curTree;
				curTree = treeIndex;
				continue;
			}
		} else if (lastIndex[curTree] == 1) {
			int i = clusters[curTree].right;
			if (i >= 0) {
				//order[curOrder++] = i;
				assignment[i] = currentCluster;
				lastIndex[curTree]++;
				continue;
			} else {
				lastIndex[curTree]++;
				int treeIndex = (-1*i)-1;
				lastTree[treeIndex] = curTree;
				curTree = treeIndex;
				continue;
			}
		} else {
			curTree = lastTree[curTree];
		}
	}

	fprintf(stderr, "\t%d Total clusters identified\n", clusterNumber-1);
	
	delete []lastTree;
	delete []lastIndex;

	return assignment;	
}



void TreeCluster::printCDT(char* prefix, int matrixFlag) {

	if (prefix == NULL) {
		prefix = (char*)"out";
	}
	char* filename = new char[10000];

	int* order = getOrder();
//fprintf(stderr, "Tree:\n");
//for (int i=0;i<numMatrix-1;i++) clusters[i].print(stderr);
//fprintf(stderr, "Order:\n");
//for (int i=0;i<numMatrix;i++) fprintf(stderr, " %d %d\n", i, order[i]);

	sprintf(filename, "%s.cdt",prefix);
	FILE* cdt = fopen(filename, "w");


	fprintf(cdt, "GID\tNA\tNAME\tGWEIGHT");
	if (matrixFlag) {
		for (int i=0;i<numMatrix;i++) {
			if (names != NULL) {	
				fprintf(cdt, "\t%s", names[order[i]]);
			} else {
				fprintf(cdt, "\t%d", order[i]);
			}
		}
	} else {
		for (int i=0;i<numExps;i++) {
			fprintf(cdt,"\t%s",expNames[i]);
		}
	}
	fprintf(cdt,"\n");


	fprintf(cdt, "EWEIGHT\t\t\t");
	if (matrixFlag) {
		for (int i=0;i<numMatrix;i++) {
			fprintf(cdt, "\t1.0");
		}
	} else {
		for (int i=0;i<numExps;i++) {
			fprintf(cdt, "\t1.0");
		}
	}
	fprintf(cdt,"\n");

	for (int i=0;i<numMatrix;i++) {
		fprintf(cdt, "GENE%dX", order[i]);
		if (names != NULL) {
			fprintf(cdt, "\t%s\t%s", names[order[i]],names[order[i]]);
		} else {
			fprintf(cdt, "\t%d\t%d", order[i],order[i]);
		}
		fprintf(cdt, "\t1.0");
		if (matrixFlag) {
			for (int j=0;j<numMatrix;j++) {
				if (ogMatrix != NULL) {
					fprintf(cdt, "\t%lf", ogMatrix[order[i]][order[j]]);	
				} else {
					fprintf(cdt, "\t0.0");
				}
			}
		} else {
			for (int j=0;j<numExps;j++) {
				fprintf(cdt, "\t%lf",dataMatrix[order[i]][j]);
			}
		}

		fprintf(cdt,"\n");
	}
	fclose(cdt);
		
	sprintf(filename, "%s.gtr",prefix);
	//double max = 1-clusters[0].distance;
	//fprintf(stderr, "max=%lf\n", max);
	//max = 1;
	FILE* gtr = fopen(filename, "w");
	double last = FLT_MAX;
	for (int i=numMatrix-2;i>rootTree;i--) {
		fprintf(gtr,"NODE%dX",i+1);
		if (clusters[i].left < 0) {
			fprintf(gtr,"\tNODE%dX",-1*(clusters[i].left));
		} else {
			fprintf(gtr,"\tGENE%dX",clusters[i].left);
		}
		if (clusters[i].right < 0) {
			fprintf(gtr,"\tNODE%dX",-1*(clusters[i].right));
		} else {
			fprintf(gtr,"\tGENE%dX",clusters[i].right);
		}
		//double outDist = (1-clusters[i].distance)/max;
		double outDist = -1*clusters[i].distance;
		//fprintf(stderr, "\t%d=%lf\n", i, clusters[i].distance);
		if (outDist > last) {
			outDist = last;
		}
		if (fixedClusterScore) {
			outDist = ((double)(numMatrix-1-i))/(double)numMatrix;
		}
		fprintf(gtr,"\t%lf\n",outDist);
		last = outDist;

	}
	for (int i=0;i<=rootTree;i++) {
		fprintf(gtr,"NODE%dX",i+1);
		if (clusters[i].left < 0) {
			fprintf(gtr,"\tNODE%dX",-1*(clusters[i].left));
		} else {
			fprintf(gtr,"\tGENE%dX",clusters[i].left);
		}
		if (clusters[i].right < 0) {
			fprintf(gtr,"\tNODE%dX",-1*(clusters[i].right));
		} else {
			fprintf(gtr,"\tGENE%dX",clusters[i].right);
		}
		//double outDist = (1-clusters[i].distance)/max;
		double outDist = -1*clusters[i].distance;
		//fprintf(stderr, "\t%d=%lf\n", i, clusters[i].distance);
		if (outDist > last) {
			outDist = last;
		}
		if (fixedClusterScore) {
			outDist = ((double)(numMatrix-1-i))/(double)numMatrix;
		}
		fprintf(gtr,"\t%lf\n",outDist);
		last = outDist;

	}
	fclose(gtr);

	delete []order;
	delete []filename;
}
void TreeCluster::printClusterSizes(FILE* fp,int maxSizeToReport) {

	double * sizes = new double[numMatrix+1];
	double * expectedSizes = new double[numMatrix+1];
	double * expectedSizes2 = new double[numMatrix+1];
	for (int i=0;i<numMatrix+1;i++) {
		sizes[i] = 0.0;
		expectedSizes[i] = 0.0;
		expectedSizes2[i] = 0.0;
	}
	sizes[1] = (double) numMatrix;
	expectedSizes[1] = (double) numMatrix;

	if (maxSizeToReport > numMatrix+1) {
		maxSizeToReport = numMatrix+1;
	}

	fprintf(fp , "ClusterRound\tClusterRound");
	for (int i=1;i<maxSizeToReport;i++) {
		fprintf(fp, "\t%d",i);
	}
	fprintf(fp, "\n");
	//fprintf(fp, "\tSum\tOdds\n");

	for (int i=0;i<numMatrix-1;i++) {
		int numLeft = 1;
		int numRight = 1;
		int numNow = (int)clusters[i].weight;
		if (clusters[i].left < 0) {
			numLeft = (int) clusters[-1*clusters[i].left-1].weight;
		}
		if (clusters[i].right < 0) {
			numRight = (int) clusters[-1*clusters[i].right-1].weight;
		}
		//fprintf(stderr, "%d %d %d %d %d\n", clusters[i].left, clusters[i].right,numLeft, numRight,numNow);
		sizes[numLeft]--;
		sizes[numRight]--;
		sizes[numNow]++;

		for (int j=1;j<numMatrix+1;j++) {
			expectedSizes2[j] = 0.0;
		}

		double odds = 0;
		double badOdds = 0;
		/*
		for (int j=1;j<numMatrix;j++) {
			for (int k=j;k<numMatrix;k++) {
				int m = j+k;
				double r = 0;
				if (j==k) {
					r = (expectedSizes[j]/((double)(numMatrix-i)))*(expectedSizes[k]/((double)(numMatrix-i)));
				} else {
					if (expectedSizes[k] > 0.0) {
						r = (expectedSizes[j]/((double)(numMatrix-i)))*((expectedSizes[k])/((double)(numMatrix-i)));
					}
				}
				odds += r;
				//fprintf(stderr, "\t\t%d\t%d\t%d\t%lf\n", i,j,k,r);
				if (m > numMatrix) {
					badOdds+=r;
				} else {
					expectedSizes2[m] += r;
					expectedSizes2[j] -= r;
					expectedSizes2[k] -= r;
				}
			}
		}*/
		double makeUpFraction = (odds-badOdds)/odds;
		for (int j=0;j<numMatrix+1;j++) {
			expectedSizes2[j] *= makeUpFraction;
		}


		//fprintf(stderr, "\t%d\t%lf\t%lf\n", i,odds,badOdds);
		for (int j=1;j<numMatrix+1;j++) {
			expectedSizes[j] += expectedSizes2[j];
			if (expectedSizes[j] < 0) expectedSizes[j] = 0.0;
			expectedSizes2[j] = 0.0;
		}
	
		fprintf(fp, "I%d\tI%d", i+1, i+1);
		/*
		if (odds < 0.99) {
			//break;
			fprintf(fp, "%dBad", i+1);
		} else {
			fprintf(fp, "%d", i+1);
		}*/
		//double pseudo = 1/((double)numMatrix)/5.0;
		for (int i=1;i<maxSizeToReport;i++) {
			double r = ((double)sizes[i]);
			/*
			double r2 = ((double)expectedSizes[i]);
			sum+=r2*i;
			sum2+=r*i;

			double logp = 0;
			if (r < r2) {
				logp = -1*ilogbinomialD(numMatrix-i,(int)r,r2/(double)numMatrix,numMatrix);
			} else {
				logp = logbinomialD(numMatrix-i,(int)r,r2/(double)numMatrix,numMatrix);
			}
			//r += pseudo;
			//r2 += pseudo;
			//double rr = log(r)/log(2.0) - log(r2)/log(2.0);
			*/
			fprintf(fp, "\t%lf", r);
		}
		//fprintf(fp, "\t%lf\t%lf\n",sum,odds);
		fprintf(fp, "\n");
	}

	delete []sizes;
	delete []expectedSizes;
	delete []expectedSizes2;
}

int* TreeCluster::getOrder() {

	int* order = new int[numMatrix];	
	int* lastTree = new int[numMatrix];	
	int* lastIndex = new int[numMatrix];	
	for (int i=0;i<numMatrix;i++) {
		order[i] = 0;
		lastTree[i] = 0;
		lastIndex[i] = 0;
	}
	int curTree = rootTree;
	//fprintf(stderr, "root = %d\n", curTree);
	int curOrder = 0;


	while (!(curTree == rootTree && lastIndex[rootTree]==2)) {
		if (lastIndex[curTree] == 0) {
			int i = clusters[curTree].left;
			if (i >= 0) {
				order[curOrder++] = i;
				lastIndex[curTree]++;
				continue;
			} else {
				lastIndex[curTree]++;
				int treeIndex = (-1*i)-1;
				lastTree[treeIndex] = curTree;
				curTree = treeIndex;
				continue;
			}
		} else if (lastIndex[curTree] == 1) {
			int i = clusters[curTree].right;
			if (i >= 0) {
				order[curOrder++] = i;
				lastIndex[curTree]++;
				continue;
			} else {
				lastIndex[curTree]++;
				int treeIndex = (-1*i)-1;
				lastTree[treeIndex] = curTree;
				curTree = treeIndex;
				continue;
			}
		} else {
			curTree = lastTree[curTree];
		}
	}

	
	delete []lastTree;
	delete []lastIndex;
	
	return order;

}


void TreeCluster::cluster() {

	double *newDistance = new double[numMatrix];

	fprintf(stderr, "\tClustering progress:\n");
	fprintf(stderr, "\t|0%%                                   50%%                                  100%%|\n\t");
	double nextMark = 0.0;
	double step = 1.0/80.0;

	int left=0;
	int right=0;

	for (int i=0;i<numMatrix-1;i++) {
		while (((double)i)/((double)numMatrix-1.0) >= nextMark) {
			nextMark+=step;
			fprintf(stderr, "=");
		}

		double min = FLT_MAX;
		//slow, safe way
		min = findMinimumPair(left,right);

		//new hotness
/*		if (i==0 || fixedPositionFlag) {
			min = initializeCache(left,right);
		} else {
			min = updateCache(left,right);
		}*/
	
		clusters[i].left = left;
		clusters[i].right = right;
		clusters[i].distance = min;
		double maxSubDistance = min;

		double leftWeight = 1.0;
		if (mask[left] < 0) {  //is node, not a leaf
			clusters[i].left = mask[left];
			int node = -1*mask[left] - 1;
			leftWeight = clusters[node].weight;
			if (clusters[node].distance > maxSubDistance) maxSubDistance = clusters[node].distance;
		}

		double rightWeight = 1.0;
		if (mask[right] < 0) { //is node, not a leaf
			clusters[i].right = mask[right];
			int node = -1*mask[right] - 1;
			rightWeight = clusters[node].weight;
			if (clusters[node].distance > maxSubDistance) maxSubDistance = clusters[node].distance;
		}
		//fprintf(stderr, "%d\t%lf\t%d\t%d\t%lf\t%lf\t%d\t%d\n", i, min,left,right,leftWeight,rightWeight,mask[left],mask[right] );



		clusters[i].distance = maxSubDistance;
		clusters[i].weight = leftWeight + rightWeight;
		for (int j=0;j<numMatrix;j++) {
			newDistance[j] =  (matrix[left][j]*leftWeight + matrix[right][j]*rightWeight)/clusters[i].weight;
		}

		mask[left] = -1*(i+1);
		mask[right] = 1;
		for (int j=0;j<numMatrix;j++) {
			matrix[left][j] = newDistance[j];	
			matrix[j][left] = newDistance[j];	
		}

		//fprintf(stderr,"%d----\n",i);
		//clusters[i].print(stderr);

	}
	fprintf(stderr, "\n");

	rootTree = numMatrix-2;
	delete []newDistance;

	//fprintf(stderr, "Done clustering\n");

}
double TreeCluster::updateCache(int &left, int &right) {
	int newLeft = 0;
	int newRight = 0;
	double gmin = FLT_MAX;	
	for (int i=0;i<numMatrix;i++) {
		int checkAllFlag = 0;
		if (mask[i]==1) continue;
		if (i == left) checkAllFlag = 1;
		if (i == right) fprintf(stderr, "SHOUDLNTget here!! This is bug, probably of the stupid kind...!!\n");
		if (minIndexCache[i] == left || minIndexCache[i] == right 
				|| mask[minIndexCache[i]]==1 || minIndexCache[i] == i) checkAllFlag=1;

		if (checkAllFlag) {	
			minDistCache[i] = FLT_MAX;
			for (int j=0;j<numMatrix;j++) {
				if (mask[j] == 1) continue;
				if (i==j) continue;
				if (matrix[i][j] < minDistCache[i]) {
					minDistCache[i] = matrix[i][j];
					minIndexCache[i] = j;
					if (matrix[i][j] < gmin) {
						if (i < j) {
							newLeft = i;
							newRight = j;
						} else {
							newLeft = j;
							newRight = i;
						}
						gmin = matrix[i][j];
					}
				}
			}
		} else {
			if (matrix[i][left] < minDistCache[i]) {
				minDistCache[i] = matrix[i][left];
				minIndexCache[i] = left;
			}
			if (minDistCache[i] < gmin) {
				if (i < minIndexCache[i]) {
					newLeft = i;
					newRight = minIndexCache[i];
				} else {
					newLeft = minIndexCache[i];
					newRight = i;
				}
				gmin = minDistCache[i];
			}
		}
	}
	left = newLeft;
	right = newRight;
	return gmin;
}
double TreeCluster::initializeCache(int &left, int &right) {
	double gmin = FLT_MAX;	
	for (int i=0;i<numMatrix;i++) {
		if (mask[i]==1) continue;
		minDistCache[i] = FLT_MAX;
		if (fixedPositionFlag) {
			for (int j=i-1;j>=0;j--) {
				if (mask[j] == 1) continue;
				if (matrix[i][j] < minDistCache[i]) {
					minDistCache[i] = matrix[i][j];
					minIndexCache[i] = j;
					if (matrix[i][j] < gmin) {
						if (i < j) {
							left = i;
							right = j;
						} else {
							left = j;
							right = i;
						}
						gmin = matrix[i][j];
					}
				}
				break;
			}
			for (int j=i+1;j<numMatrix;j++) {
				if (mask[j] == 1) continue;
				if (matrix[i][j] < minDistCache[i]) {
					minDistCache[i] = matrix[i][j];
					minIndexCache[i] = j;
					if (matrix[i][j] < gmin) {
						if (i < j) {
							left = i;
							right = j;
						} else {
							left = j;
							right = i;
						}
						gmin = matrix[i][j];
					}
				}
				break;
			}
		} else {
			for (int j=0;j<numMatrix;j++) {
				if (mask[j] == 1) continue;
				if (i==j) continue;
				if (matrix[i][j] < minDistCache[i]) {
					minDistCache[i] = matrix[i][j];
					minIndexCache[i] = j;
					if (matrix[i][j] < gmin) {
						if (i < j) {
							left = i;
							right = j;
						} else {
							left = j;
							right = i;
						}
						gmin = matrix[i][j];
					}
				}
			}
		}
	}
	return gmin;
}

double TreeCluster::findMinimumPair(int &left, int &right) {

	double min = FLT_MAX;	
	left = 0;
	right = 0;
	for (int i=0;i<numMatrix;i++) {
		if (mask[i]==1) continue;
		for (int j=i+1;j<numMatrix;j++) {
			if (mask[j] == 1) continue;
			if (matrix[i][j] < min) {
				left = i;
				right = j;
				min = matrix[i][j];
			}
			if (fixedPositionFlag) break;
		}
	}
	//fprintf(stderr, "Minimum: %d %d %lf\n", left, right, min);
	return min;
}

void TreeCluster::scoreClusterFiles(char** files, int numFiles, FILE* fp) {

	Inttable* regionIndex = new Inttable(numMatrix);
	for (int i=0;i<numMatrix;i++) {
		regionIndex->insert(i,names[i]);
	}
	int bufSize = 1000000;
	char* buf = new char[bufSize];
	char** line = new char*[100000];
	int numCols = 0;


	int** clusterIndex = new int*[numFiles];
	int* clusterSize = new int[numFiles];

	for (int i=0;i<numFiles;i++) {
		FILE* fp = fopen(files[i],"r");
		if (fp == NULL) {
			fprintf(stderr, "!!! Could not open %s for reading!!!\n", files[i]);
			exit(0);
		}
		clusterSize[i] = 0;
		clusterIndex[i] = new int[numMatrix];
		while (fgets(buf,bufSize,fp) != NULL) {
			clusterSplit(buf, line, numCols, '\t');
			int index = regionIndex->search(line[0]);
			if (index != EMPTY_INT) {
				clusterIndex[i][clusterSize[i]++] = index;
			}
			//fprintf(stderr, "\t\t%s\t%d\n", line[0],index);
			if (clusterSize[i] >= numMatrix) {
				fprintf(stderr,"!!! Error - cluster size for file %s exceeds # of regions\n",files[i]);
				exit(1);
			}
		}
		fclose(fp);
		fprintf(stderr, "\t%s: %d regions\n", files[i], clusterSize[i]);
	}	


	fprintf(fp, "Region Name");
	for (int i=0;i<numFiles;i++) {
		fprintf(fp, "\t%s", files[i]);
	}
	fprintf(fp, "\n");
	for (int i=0;i<numMatrix;i++) {
		fprintf(fp, "%s",names[i]);
		for (int j=0;j<numFiles;j++) {
			double score = 0.0;
			double N = 0.0;
			for (int k=0;k<clusterSize[j];k++) {
				score += matrix[i][clusterIndex[j][k]];
				N += 1.0;
			}
			if (N > 0.0) {
				score /= N;
			}
			fprintf(fp, "\t%lf",score);
		}
		fprintf(fp, "\n");
	}

	for (int i=0;i<numFiles;i++) {
		delete [](clusterIndex[i]);
	}
	delete []clusterIndex;
	delete []clusterSize;
	delete regionIndex;
	delete []buf;
	delete []line;

}







TreeNode::TreeNode() {
	left = 0;
	right = 0;
	distance = 0.0;
	weight = 0.0;
	matrixIndex = 0;
}

TreeNode::~TreeNode() {
}

void TreeNode::print(FILE* fp) {
	fprintf(fp, "\tNode:\n");
	fprintf(fp,"\t\tleft:%d\n", left);
	fprintf(fp,"\t\tright:%d\n", right);
	fprintf(fp,"\t\tdistance:%lf\n", distance);
	fprintf(fp,"\t\tweight:%lf\n", weight);
	fprintf(fp,"\t\tmatrixIndex:%d\n", matrixIndex);
}
void clusterSplit(char* string, char** cols, int &numCols, char delim) {
	cols[0] = string;
	numCols=1;
	int len = strlen(string);
	for (int i=0;i<len;i++) {
		if (string[i] == delim) {
			string[i] = '\0';
			cols[numCols] = &(string[i+1]);
			numCols++;
		} else if (string[i] == '\n') {
			string[i] = '\0';
		} else if (string[i] == '\r') {
			string[i] = '\0';
		}
	}
}

