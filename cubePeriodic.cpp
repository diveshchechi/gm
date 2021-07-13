// debugTutorial.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
// 3D calculation
// No orientation dependence
#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <string.h>
#include<iostream>
#include <string>

#define DRND(x) ((double)(x)/RAND_MAX*rand())

#define NDX 32 //64
#define NDY 32 //128
#define NDZ 32 //64
#define N 10
#define GNP 10
#define INX 400						% Drawing window 1 Pixel size of side x
#define INY 800						% Drawing window 1 Pixel size of side y
#define INZ 400						% Drawing window 1 Pixel size of side z

int ndx = NDX, ndxm = NDX - 1;
int ndy = NDY, ndym = NDY - 1;
int ndz = NDZ, ndzm = NDZ - 1;
int nm = N - 1, nmm = N - 2, GN = GNP - 1;
double PI = 3.141592, time1;
double RR = 8.3145;

double ph[N][NDX][NDY][NDZ];		//% Ph [1] [i] [j] to ph [n] [i] [j]: pf in position [i] [j]
int qh[N][NDX][NDY][NDZ];	 		//% Qh [1] [i] [j] to qh [n] [i] [j]: grain number in position [i] [j]
int n00[NDX][NDY][NDZ];				//% The number of cases where ph is not 0 in position [i] [j]
int n00p[NDX][NDY][NDZ];            //% The number of cases when ph in position [i] [j] and ph around it are not 0
int is[NDX][NDY][NDZ];
double ph2[N][NDX][NDY][NDZ];
int qh2[N][NDX][NDY][NDZ];
int seeds[NDX][NDY][NDZ];
int centroid[N][3];
int seedVolume[N];
int lloyd = 1;

double aij[GNP][GNP], wij[GNP][GNP], tij[GNP][GNP], eij[GNP][GNP];

using namespace std;
//*********** Initialize *******************************
void ini000() {
	
	int i, j, k, l, it;
	int ii, jj, kk, kkk, iGN;
	int ip, im, jp, jm, kp, km;
	double sum1, t;
	printf("India\n");
	double r;
	int x1, y1, z1, r0;
	
	
	
	//%srand(time(NULL)); % Random number initialization
	
	for (i = 0; i <= ndxm; i++) {
		for (j = 0; j <= ndym; j++) {
			for (k = 0; k <= ndzm; k++) {
				ph[1][i][j][k] = 1.0;  //qh[1][i][j][k]=GN;
				n00[i][j][k] = 1;      n00p[i][j][k] = 1;
				is[i][j][k] = GN;
			}
		}
	}
	
	r0 = 5;
	//%Write a sphere with a radius 5 with different grain number (grain number written in circle)
	for (iGN = 1; iGN <= GN - 1; iGN++) {
		x1 = ndx * DRND(1); y1 = ndy * DRND(1); z1 = ndz * DRND(1);
		while (is[x1][y1][z1] != GN) { x1 = ndx * DRND(1); y1 = ndy * DRND(1); z1 = ndz * DRND(1); }
		for (i = -r0; i <= (ndxm + r0); i++) {
			ii = i; if (i < 0) { ii = ndx + i; }  if (i > ndxm) { ii = i - ndx; }
			for (j = -r0; j <= (ndym + r0); j++) {
				jj = j; if (j < 0) { jj = ndy + j; }  if (j > ndym) { jj = j - ndy; }
				for (k = -r0; k <= (ndzm + r0); k++) {
					kk = k; if (k < 0) { kk = ndz + k; }  if (k > ndzm) { kk = k - ndz; }
					r = sqrt((double(i - x1)) * (double(i - x1))
						+ (double(j - y1)) * (double(j - y1))
						+ (double(k - z1)) * (double(k - z1)));
					if (r <= r0) { is[ii][jj][kk] = iGN; }
				}
			}
		}
	}

	/*
	r0 = 5.0; x1 = ndx / 2; y1 = ndy / 2; z1 = ndz / 2;
	for (i = 0; i <= ndxm; i++) {
		for (j = 0; j <= ndym; j++) {
			for (k = 0; k <= ndzm; k++) {
				r = sqrt((double(i - x1)) * (double(i - x1))
					+ (double(j - y1)) * (double(j - y1))
					+ (double(k - z1)) * (double(k - z1)));
				//if(r<=r0){ is[i][j][k]=1; } 
			}
		}
	}
	*/

	for (i = 0; i <= ndxm; i++) {
		for (j = 0; j <= ndym; j++) {
			for (k = 0; k <= ndzm; k++) {
				qh[1][i][j][k] = is[i][j][k];
			}
		}
	}
	
	//--- ph=0 Special case  (if there are other orientations around) added -------------------------
	for (i = 0; i <= ndxm; i++) {
		for (j = 0; j <= ndym; j++) {
			for (k = 0; k <= ndzm; k++) {

				ip = i + 1; im = i - 1;  jp = j + 1;  jm = j - 1;  kp = k + 1;  km = k - 1;
				if (i == ndxm) { ip = 0; }  if (i == 0) { im = ndxm; }
				if (j == ndym) { jp = 0; }  if (j == 0) { jm = ndym; }
				if (k == ndzm) { kp = 0; }  if (k == 0) { km = ndzm; }

				for (kk = 1; kk <= n00[ip][j][k]; kk++) {
					kkk = qh[kk][ip][j][k]; //The kk-th  in position [ip] [j] [k]
					ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
					//% If direction kkk lies in the direction of position [i] [j] [k], ii = 0, if none, ii = 1
					if (ii == 1) {
						n00p[i][j][k] += 1;
						ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
					}
				}

				for (kk = 1; kk <= n00[im][j][k]; kk++) {
					kkk = qh[kk][im][j][k];
					ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
					if (ii == 1) {
						n00p[i][j][k] += 1;
						ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
					}
				}

				for (kk = 1; kk <= n00[i][jp][k]; kk++) {
					kkk = qh[kk][i][jp][k];
					ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
					if (ii == 1) {
						n00p[i][j][k] += 1;
						ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
					}
				}

				for (kk = 1; kk <= n00[i][jm][k]; kk++) {
					kkk = qh[kk][i][jm][k];
					ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
					if (ii == 1) {
						n00p[i][j][k] += 1;
						ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
					}
				}

				for (kk = 1; kk <= n00[i][j][kp]; kk++) {
					kkk = qh[kk][i][j][kp];
					ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
					if (ii == 1) {
						n00p[i][j][k] += 1;
						ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
					}
				}

				for (kk = 1; kk <= n00[i][j][km]; kk++) {
					kkk = qh[kk][i][j][km];
					ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
					if (ii == 1) {
						n00p[i][j][k] += 1;
						ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
					}
				}

			}
		}
	}
}

void findNewPosition() {
	// for number of seeds - 2
	int i, j, k, l, it;
	int x1, y1, z1, iGN, r0, r;
	int ii, jj, kk, kkk;
	int ip, im, jp, jm, kp, km;
	double sum1, t;

	for (i = 0; i <= ndxm; i++) {
		for (j = 0; j <= ndym; j++) {
			for (k = 0; k <= ndzm; k++) {
				int seed = seeds[i][j][k];
				centroid[seed][0] += i;
				centroid[seed][1] += j;
				centroid[seed][2] += k;
				seedVolume[seed] ++;
				
				// need to reset centroid to zero later seedvolume also to zero
			}
		}
	}

	for (i = 1; i < N - 1; i++) {
		centroid[i][0] /= seedVolume[i];
		centroid[i][1] /= seedVolume[i];
		centroid[i][2] /= seedVolume[i];
		printf("For seed %d of volume %d new location is %d %d %d\n", i, seedVolume[i], centroid[i][0], centroid[i][1], centroid[i][2]);
	}

	for (i = 0; i <= ndxm; i++) {
		for (j = 0; j <= ndym; j++) {
			for (k = 0; k <= ndzm; k++) {
				ph[1][i][j][k] = 1.0;  //qh[1][i][j][k]=GN;
				n00[i][j][k] = 1;      n00p[i][j][k] = 1;
				is[i][j][k] = GN;
			}
		}
	}

	r0 = 5;

	//%Write a sphere with a radius 5 with different grain number (grain number written in circle)
	for (iGN = 1; iGN <= GN - 1; iGN++) {
		x1 = centroid[iGN][0];
		y1 = centroid[iGN][1];
		z1 = centroid[iGN][2];
		for (i = -r0; i <= (ndxm + r0); i++) {
			ii = i; if (i < 0) { ii = ndx + i; }  if (i > ndxm) { ii = i - ndx; }
			for (j = -r0; j <= (ndym + r0); j++) {
				jj = j; if (j < 0) { jj = ndy + j; }  if (j > ndym) { jj = j - ndy; }
				for (k = -r0; k <= (ndzm + r0); k++) {
					kk = k; if (k < 0) { kk = ndz + k; }  if (k > ndzm) { kk = k - ndz; }
					r = sqrt((double(i - x1)) * (double(i - x1))
						+ (double(j - y1)) * (double(j - y1))
						+ (double(k - z1)) * (double(k - z1)));
					if (r <= r0) { is[ii][jj][kk] = iGN; }
				}
			}
		}
	}

	for (i = 0; i <= ndxm; i++) {
		for (j = 0; j <= ndym; j++) {
			for (k = 0; k <= ndzm; k++) {
				for (kk = 0; kk < GN; kk++) {
					qh[kk][i][j][k] = 0;
					qh2[kk][i][j][k] = 0;
					ph2[kk][i][j][k] = 0;
				}
				qh[1][i][j][k] = is[i][j][k];
			}
		}
	}

	// resetiing centroid and seed volume
	for (i = 1; i < N - 1; i++) {
		centroid[i][0] = 0;
		centroid[i][1] = 0;
		centroid[i][2] = 0;
		seedVolume[i] = 1;
		printf("For seed %d new location is %d %d %d\n", i, centroid[i][0], centroid[i][1], centroid[i][2]);
	}

	//--- ph=0 Special case  (if there are other orientations around) added -------------------------
	for (i = 0; i <= ndxm; i++) {
		for (j = 0; j <= ndym; j++) {
			for (k = 0; k <= ndzm; k++) {

				ip = i + 1; im = i - 1;  jp = j + 1;  jm = j - 1;  kp = k + 1;  km = k - 1;
				if (i == ndxm) { ip = 0; }  if (i == 0) { im = ndxm; }
				if (j == ndym) { jp = 0; }  if (j == 0) { jm = ndym; }
				if (k == ndzm) { kp = 0; }  if (k == 0) { km = ndzm; }

				for (kk = 1; kk <= n00[ip][j][k]; kk++) {
					kkk = qh[kk][ip][j][k]; //The kk-th  in position [ip] [j] [k]
					ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
					//% If direction kkk lies in the direction of position [i] [j] [k], ii = 0, if none, ii = 1
					if (ii == 1) {
						n00p[i][j][k] += 1;
						ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
					}
				}

				for (kk = 1; kk <= n00[im][j][k]; kk++) {
					kkk = qh[kk][im][j][k];
					ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
					if (ii == 1) {
						n00p[i][j][k] += 1;
						ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
					}
				}

				for (kk = 1; kk <= n00[i][jp][k]; kk++) {
					kkk = qh[kk][i][jp][k];
					ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
					if (ii == 1) {
						n00p[i][j][k] += 1;
						ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
					}
				}

				for (kk = 1; kk <= n00[i][jm][k]; kk++) {
					kkk = qh[kk][i][jm][k];
					ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
					if (ii == 1) {
						n00p[i][j][k] += 1;
						ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
					}
				}

				for (kk = 1; kk <= n00[i][j][kp]; kk++) {
					kkk = qh[kk][i][j][kp];
					ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
					if (ii == 1) {
						n00p[i][j][k] += 1;
						ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
					}
				}

				for (kk = 1; kk <= n00[i][j][km]; kk++) {
					kkk = qh[kk][i][j][km];
					ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
					if (ii == 1) {
						n00p[i][j][k] += 1;
						ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
					}
				}

			}
		}
	}

}


//*********** Initial Data in **************************
void datin()
{
	FILE* datin0;
	int 		i, j, k, kk;

	datin0 = fopen("test.ini", "r");

	for (i = 0; i <= ndxm; i++) {
		for (j = 0; j <= ndym; j++) {
			for (k = 0; k <= ndzm; k++) {
				for (kk = 1; kk <= nm; kk++) {
					fscanf(datin0, "%d  %lf  ", &qh[kk][i][j][k], &ph[kk][i][j][k]);
				}
			}
		}
	}
	fclose(datin0);
}


//************ Data Save *******************************
void datsave(int fNum, int iter )
{
	FILE* stream;
	int 		i, j, k, kk;
	double 	col;
	//string fName = "C:\\Users\\dk18\\Work\\hiwi\\files\\29June\\cuboidRepositionDense\\test" + to_string(iter) + "num" +  to_string(fNum)+".txt";
	string fName = "C:\\Users\\dk18\\Work\\hiwi\\files\\6July\\cubeIter8seed\\test" + to_string(fNum) + ".txt";
	char* char_arr;
	char_arr = &fName[0];
	stream = fopen(char_arr, "a");
	//fprintf(stream, "%e  \n", time1);
	fprintf(stream, "xcoord,ycoord,zcoord,field,seed");
	for (i = 0; i <= ndxm; i++) {
		for (j = 0; j <= ndym; j++) {
			for (k = 0; k <= ndzm; k++) {
				/*
				fprintf(stream, "%d  \n", n00[i][j][k]);
				for (kk = 1; kk <= n00[i][j][k]; kk++) {
					fprintf(stream, "%d , %e  ", qh[kk][i][j][k], ph[kk][i][j][k]);
				}*/
				//fprintf(stream, "\n");
				col=0.; for(kk=1;kk<=n00[i][j][k];kk++){ col+=ph[kk][i][j][k]*ph[kk][i][j][k]; }
				//fprintf(stream, "\n%d,%d,%d,%e,%d,%d", i,j,k,col,is[i][j][k],qh[i][j][k]);
				double max = -1;
				int seed = 0;
				for (kk = 1; kk <= n00[i][j][k]; kk++) {
					if (max < ph[kk][i][j][k]) {
						max = ph[kk][i][j][k];
						seed = qh[kk][i][j][k];
						seeds[i][j][k] = seed;
					}
					//fprintf(stream, "%d  %e\n ", qh[kk][i][j][k], ph[kk][i][j][k]);
				}
				fprintf(stream, "\n%d,%d,%d,%e,%dseed", i, j, k, col, seed);
				//fprintf(stream, "\n%d,%d,%d,%d", i, j, k, is[i][j][k]);
			}
		}
	}
	//fprintf(stream, "---------------------------------------------------------------------------\n");
	fclose(stream);
}


void dataVtkVolMesh(int fNum)
{
	FILE* stream;
	int 		i, j, k, kk;
	double 	col;
	string fName = "C:\\Users\\dk18\\Work\\hiwi\\files\\vtkmushroom\\data" + to_string(fNum) + ".vtk";
	char* char_arr;
	char_arr = &fName[0];
	stream = fopen(char_arr, "a");
	//fprintf(stream, "%e  \n", time1);
	//int pts = ndxm * ndym * ndzm;
	int pts = ndxm * ndym * ndzm;
	fprintf(stream, "# vtk DataFile Version 2.0\n");
	fprintf(stream, "Mesh Field\n");
	fprintf(stream, "ASCII\n");
	fprintf(stream, "DATASET STRUCTURED_GRID\n");
	fprintf(stream, "DIMENSIONS %d %d %d\n", ndxm, ndym, ndzm);
	fprintf(stream, "POINTS %d int\n", pts);

	for (i = 0; i < ndxm; i++) {
		for (j = 0; j < ndym; j++) {
			for (k = 0; k < ndzm; k++) {
				fprintf(stream, "%d %d %d\n", i / 4, j / 4, k / 4);
				/*
				fprintf(stream, "%d  \n", n00[i][j][k]);
				for (kk = 1; kk <= n00[i][j][k]; kk++) {
					fprintf(stream, "%d , %e  ", qh[kk][i][j][k], ph[kk][i][j][k]);
				}*/
				//fprintf(stream, "\n");



			}
		}
	}
	fprintf(stream, "POINT_DATA %d\n", pts);
	fprintf(stream, "FIELD Field_Data 1\n");
	fprintf(stream, "field 1 %d double\n", pts);
	for (i = 0; i < ndxm; i++) {
		for (j = 0; j < ndym; j++) {
			for (k = 0; k < ndzm; k++) {
				/*
				fprintf(stream, "%d  \n", n00[i][j][k]);
				for (kk = 1; kk <= n00[i][j][k]; kk++) {
					fprintf(stream, "%d , %e  ", qh[kk][i][j][k], ph[kk][i][j][k]);
				}*/
				//fprintf(stream, "\n");
				//float v = 0.23 * k;
				col = 0.; for (kk = 1; kk <= n00[i][j][k]; kk++) { col += ph[kk][i][j][k] * ph[kk][i][j][k]; }


				double max = -1;
				int seed = 0;
				for (kk = 1; kk <= n00[i][j][k]; kk++) {
					if (max < ph[kk][i][j][k]) {
						max = ph[kk][i][j][k];
						seed = qh[kk][i][j][k];
					}
					//fprintf(stream, "%d  %e\n ", qh[kk][i][j][k], ph[kk][i][j][k]);
				}
				fprintf(stream, "%d\n", seed);
				
			}
		}
	}
	//fprintf(stream, "---------------------------------------------------------------------------\n");
	fclose(stream);
}



//************ Data Save in VTK Format *****************
void dataVtk(int fNum)
{
	FILE* stream;
	int 		i, j, k, kk;
	double 	col;
	string fName = "C:\\Users\\dk18\\Work\\hiwi\\files\\cubeVtk\\Cube" + to_string(fNum) + ".vtk";
	char* char_arr;
	char_arr = &fName[0];
	stream = fopen(char_arr, "a");
	//fprintf(stream, "%e  \n", time1);
	//int pts = ndxm * ndym * ndzm;
	int pts = ndxm * ndym * ndzm;
	fprintf(stream, "# vtk DataFile Version 2.0\n");
	fprintf(stream, "Mesh Field\n");
	fprintf(stream, "ASCII\n");
	fprintf(stream, "DATASET STRUCTURED_GRID\n");
	fprintf(stream, "DIMENSIONS %d %d %d\n", ndxm, ndym, ndzm);
	fprintf(stream, "POINTS %d int\n", pts);

	for (i = 0; i < ndxm; i++) {
		for (j = 0; j < ndym; j++) {
			for (k = 0; k < ndzm; k++) {
				fprintf(stream, "%d %d %d\n", i, j, k);
				/*
				fprintf(stream, "%d  \n", n00[i][j][k]);
				for (kk = 1; kk <= n00[i][j][k]; kk++) {
					fprintf(stream, "%d , %e  ", qh[kk][i][j][k], ph[kk][i][j][k]);
				}*/
				//fprintf(stream, "\n");



			}
		}
	}
	fprintf(stream, "POINT_DATA %d\n", pts);
	fprintf(stream, "FIELD Field_Data 1\n");
	fprintf(stream, "field 1 %d double\n", pts);
	for (i = 0; i < ndxm; i++) {
		for (j = 0; j < ndym; j++) {
			for (k = 0; k < ndzm; k++) {
				/*
				fprintf(stream, "%d  \n", n00[i][j][k]);
				for (kk = 1; kk <= n00[i][j][k]; kk++) {
					fprintf(stream, "%d , %e  ", qh[kk][i][j][k], ph[kk][i][j][k]);
				}*/
				//fprintf(stream, "\n");
				//float v = 0.23 * k;
				if (1) {
					col = 0.; for (kk = 1; kk <= n00[i][j][k]; kk++) { col += ph[kk][i][j][k] * ph[kk][i][j][k]; }
					fprintf(stream, "%f\n", col);
				}
				else {
					fprintf(stream, "0\n");
				}
			}
		}
	}
	//fprintf(stream, "---------------------------------------------------------------------------\n");
	fclose(stream);
}




int main() {
	int i, j, k, l, ii, jj, kk, ll, it, kkk;
	int ip, im, jp, jm, kp, km;
	int n1, n2, n3;
	int n0p;
	double time1max;
	double delt, L, b1, t, s;
	double dG, M1, W1, K1, E1;
	double mu_chem, mu_grad;
	double temp;
	double sum1, sum2, sum3, sxs;
	double pddtt;

	double gamma, delta, amobi;
	double aaa, vm0;

	//****** reg data ****************************************

	//printf("delt(2.0)=  "); scanf(" %lf",&delt);
	delt = 2.0;

	temp = 1000.0;
	L = 500.0;
	b1 = L / (double)NDX * 1.0e-9;

	vm0 = 7.0e-6;
	gamma = 0.5 * vm0 / RR / temp / b1;
	delta = 3.0;

	aaa = 2.0 / PI * sqrt(2.0 * delta * gamma);
	K1 = aaa * aaa;
	W1 = 4.0 * gamma / delta;
	amobi = 1.;
	M1 = amobi * PI * PI / (8.0 * delta);
	E1 = 500.0 / RR / temp;
	//E1=0.0;

	time1 = 0.;  time1max = 1.0e+20 + 1.0;

	for (i = 1; i <= GN; i++) {
		for (j = 1; j <= GN; j++) {
			wij[i][j] = W1;
			aij[i][j] = K1;
			tij[i][j] = M1;
			eij[i][j] = 0.0;	
			if ((i == GN) || (j == GN)) { eij[i][j] = E1; }
			if (i > j) { eij[i][j] = -eij[i][j]; }
			if (i == j) { wij[i][j] = 0.0;  aij[i][j] = 0.0;  tij[i][j] = 0.0; eij[i][j] = 0.0; }
		}
	}

	ini000();
	//datsave();
	//datin();
	int counting = 1000;
	datsave(counting,lloyd);
	counting++;
start:
	//if ((((int)(time1) % 100) == 0)) { datsave(); }
	//datin();
	//datsave();

	for (i = 0; i <= ndxm; i++) {
		for (j = 0; j <= ndym; j++) {
			for (k = 0; k <= ndzm; k++) {
				ip = i + 1; im = i - 1;  jp = j + 1;  jm = j - 1;  kp = k + 1;  km = k - 1;
				if (i == ndxm) { ip = 0; }  if (i == 0) { im = ndxm; }
				if (j == ndym) { jp = 0; }  if (j == 0) { jm = ndym; }
				if (k == ndzm) { kp = 0; }  if (k == 0) { km = ndzm; }

				n0p = n00p[i][j][k];
				if (n0p >= 2) {
					for (n1 = 1; n1 <= n0p; n1++) {
						ii = qh[n1][i][j][k];  pddtt = 0.0;
						for (n2 = 1; n2 <= n0p; n2++) {
							jj = qh[n2][i][j][k];  sum1 = 0.0;
							for (n3 = 1; n3 <= n0p; n3++) {
								kk = qh[n3][i][j][k]; sum2 = 0.0;
								for (l = 1; l <= n00[ip][j][k]; l++) { if (qh[l][ip][j][k] == kk) { sum2 += ph[l][ip][j][k]; } }
								for (l = 1; l <= n00[im][j][k]; l++) { if (qh[l][im][j][k] == kk) { sum2 += ph[l][im][j][k]; } }
								for (l = 1; l <= n00[i][jp][k]; l++) { if (qh[l][i][jp][k] == kk) { sum2 += ph[l][i][jp][k]; } }
								for (l = 1; l <= n00[i][jm][k]; l++) { if (qh[l][i][jm][k] == kk) { sum2 += ph[l][i][jm][k]; } }
								for (l = 1; l <= n00[i][j][kp]; l++) { if (qh[l][i][j][kp] == kk) { sum2 += ph[l][i][j][kp]; } }
								for (l = 1; l <= n00[i][j][km]; l++) { if (qh[l][i][j][km] == kk) { sum2 += ph[l][i][j][km]; } }

								sum1 += 0.5 * (aij[ii][kk] - aij[jj][kk]) * (sum2 - 6.0 * ph[n3][i][j][k])
									+ (wij[ii][kk] - wij[jj][kk]) * ph[n3][i][j][k];
							}
							pddtt += -2.0 * tij[ii][jj] / double(n00p[i][j][k])
								* (sum1 - 8.0 / PI * eij[ii][jj] * sqrt(ph[n1][i][j][k] * ph[n2][i][j][k]));
						}

						ph2[n1][i][j][k] = ph[n1][i][j][k] + pddtt * delt;
						qh2[n1][i][j][k] = qh[n1][i][j][k];
						if (ph2[n1][i][j][k] >= 1.0) { ph2[n1][i][j][k] = 1.0;  }
						if (ph2[n1][i][j][k] <= 0.0) { ph2[n1][i][j][k] = 0.0; }
						//if(ph2[n1][i][j][k]<=1.0e-4){ph2[n1][i][j][k]=0.0;}
					}
				}
				else { ph2[1][i][j][k] = ph[1][i][j][k]; qh2[1][i][j][k] = qh[1][i][j][k]; }//if
			}//k
		}//j
	}//i
	double maxSeed;
	//%--- ph2?ph, qh2?qh and ph2=0,qh2=0 -------------------------
	for (i = 0; i <= ndxm; i++) {
		for (j = 0; j <= ndym; j++) {
			for (k = 0; k <= ndzm; k++) {
				maxSeed = -100;
				for (kk = 1; kk <= nm; kk++) {
					//------------------experiment
					if (maxSeed < ph2[kk][i][j][k]) {
						//maxSeed = ph2[kk][i][j][k];
						//f << kk << endl;
						//is[i][j][k] = kk;
					}
					ph[kk][i][j][k] = ph2[kk][i][j][k];  qh[kk][i][j][k] = qh2[kk][i][j][k];
					ph2[kk][i][j][k] = 0.0;
				}
			}
		}
	}

	//%--- Sort (descending) -------------------------
	for (i = 0; i <= ndxm; i++) {
		for (j = 0; j <= ndym; j++) {
			for (k = 0; k <= ndzm; k++) {

				for (kk = 1; kk <= nm - 2; kk++) { //Sort from kk = 1 to kk = nm in descending order (large to small) (kk = 0 is ignored)
					for (l = nm - 1; l > kk; l--) {
						if (ph[l][i][j][k] > ph[l - 1][i][j][k]) {
							t = ph[l][i][j][k];  ph[l][i][j][k] = ph[l - 1][i][j][k]; ph[l - 1][i][j][k] = t;
							it = qh[l][i][j][k]; qh[l][i][j][k] = qh[l - 1][i][j][k]; qh[l - 1][i][j][k] = it;
						}
					}
				}

				//%--- Standardization-------------------------
				sum1 = 0.0; ii = 0;
				for (kk = 1; kk <= nm; kk++) { if (ph[kk][i][j][k] > 0.0) { ii++; sum1 += ph[kk][i][j][k]; } }
				n00[i][j][k] = ii; n00p[i][j][k] = ii;
				for (kk = 1; kk <= n00[i][j][k]; kk++) { ph[kk][i][j][k] = ph[kk][i][j][k] / sum1; }
				//for(kk=1;kk<=n00[i][j][k];kk++){ if(sum1>=1.0){ph[kk][i][j][k]=ph[kk][i][j][k]/sum1;} }

			}
		}
	}


	//%--- Added special case of ph = 0 (if there are other orientations around)-------------------------
	for (i = 0; i <= ndxm; i++) {
		for (j = 0; j <= ndym; j++) {
			for (k = 0; k <= ndzm; k++) {

				ip = i + 1; im = i - 1;  jp = j + 1;  jm = j - 1;  kp = k + 1;  km = k - 1;
				if (i == ndxm) { ip = 0; }  if (i == 0) { im = ndxm; }
				if (j == ndym) { jp = 0; }  if (j == 0) { jm = ndym; }
				if (k == ndzm) { kp = 0; }  if (k == 0) { km = ndzm; }

				for (kk = 1; kk <= n00[ip][j][k]; kk++) {
					kkk = qh[kk][ip][j][k]; //%The kth grain number in position [ip] [j] [k]
					ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
					//%If direction kkk lies in the direction of position [i] [j] [k], ii = 0, if none, ii = 1
					if (ii == 1) {
						n00p[i][j][k] += 1;
						ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
					}
				}

				for (kk = 1; kk <= n00[im][j][k]; kk++) {
					kkk = qh[kk][im][j][k];
					ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
					if (ii == 1) {
						n00p[i][j][k] += 1;
						ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
					}
				}

				for (kk = 1; kk <= n00[i][jp][k]; kk++) {
					kkk = qh[kk][i][jp][k];
					ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
					if (ii == 1) {
						n00p[i][j][k] += 1;
						ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
					}
				}

				for (kk = 1; kk <= n00[i][jm][k]; kk++) {
					kkk = qh[kk][i][jm][k];
					ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
					if (ii == 1) {
						n00p[i][j][k] += 1;
						ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
					}
				}

				for (kk = 1; kk <= n00[i][j][kp]; kk++) {
					kkk = qh[kk][i][j][kp];
					ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
					if (ii == 1) {
						n00p[i][j][k] += 1;
						ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
					}
				}

				for (kk = 1; kk <= n00[i][j][km]; kk++) {
					kkk = qh[kk][i][j][km];
					ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
					if (ii == 1) {
						n00p[i][j][k] += 1;
						ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
					}
				}

			}
		}
	}
	//%******[Time increase]*************************************************
	if (counting < 1003) {
		datsave(counting,lloyd);
	}
	else if (counting % 5 == 0) {
		datsave(counting,lloyd);
	}
	if ((counting-1000)%250 == 0) {
		// calculate new position of seeds
		findNewPosition();
		lloyd++;
		//counting = 1000;
		counting++;
		datsave(counting, lloyd);
		if (lloyd == 20) return 0;
		goto start;
		return 0;
	}
	counting++;
	goto start;
	return 0;
}
// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
