void meshin()
{
	FILE* datin0;
	int 		i, j, k, kk;

	datin0 = fopen("C:\\Users\\dk18\\Work\\hiwi\\files\\eightparam\\eighty.txt", "r");
	//datin0 = fopen("C:\\Users\\dk18\\Work\\hiwi\\files\\mushroom\\mushroom.txt", "r");
	int x, y, z;
	x = 32;
	y = 32;
	//z = 32;
	// For eightparam  z = 64
	z = 64;
	for (i = 0; i < x; i++) {
		for (j = 0; j < y; j++) {
			for (k = 0; k < z; k++) {
				fscanf(datin0, "%d  ", &inside[i][j][k]);
			}
		}
	}
	for (i = 0; i < x; i++) {
		for (j = 0; j < y; j++) {
			for (k = 0; k < z; k++) {
				fscanf(datin0, "%lf %lf %lf", &boxCoordinates[i][j][k][0], &boxCoordinates[i][j][k][1], &boxCoordinates[i][j][k][2]);
			}
		}
	}
	fclose(datin0);
}
