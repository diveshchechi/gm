// These are global variables reqd
int seeds[NDX][NDY][NDZ];
int centroid[N][3];
int seedVolume[N];
int lloyd = 1;
//************ Data Save *******************************
void datsave(int fNum, int iter )
{
	FILE* stream;
	int 		i, j, k, kk;
	double 	col;
	//string fName = "C:\\Users\\dk18\\Work\\hiwi\\files\\29June\\cuboidRepositionDense\\test" + to_string(iter) + "num" +  to_string(fNum)+".txt";
	string fName = "C:\\Users\\dk18\\Work\\hiwi\\files\\6July\\cubeIter\\test" + to_string(fNum) + ".txt";
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

