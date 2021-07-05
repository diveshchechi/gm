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
		//printf("For seed %d new location is %d %d %d\n", i, centroid[i][0], centroid[i][1], centroid[i][2]);
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
