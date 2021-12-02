#include <stdlib.h>
#include <stdio.h>

#define LEN 256

int main(){
	
	// Dilwsi metavlitwn
		// Dilwsi range
		int max_num = 3, max_numB = 10;
		int min_num = -3, min_numB = -10;
		int max = 2*max_num, min = 0, maxB = 2*max_numB, minB = 0;
		
		// Dilwsi upoloipwn metavlitwn
		int x, y, i, j, k = 0, n;

	// Diavasma megethous tou pinaka
	printf("Vale to megethos tou pinaka: ");
	scanf("%d", &n);
	printf("\n");
		
	// Dilwsi teleutaiwn metavlitwn
	int ni = n;
	int nj = n;
		
	// Dilwsi tou pinaka
	int pinakas[n][n], b[n];
		
	// Gemisma tou pinakas[i][j] kai tou b[i][j] me tuxaies times
		for (i=0;i<ni;i++){
			for (j=0;j<nj;j++){
				if(i==j){
					pinakas[i][j] = 0;	
				}else{
					x = rand()%(max + 1 - min) + min;
					pinakas[i][j] = x-3;
				}
			}
			y = rand()%(maxB + 1 - minB) + minB;
			b[i] = y - 10;
		}
		
	// Make strongly diagwnio pinaka
		for (i=0;i<ni;i++){
			k=0;
			for (j=0;j<nj;j++){
				if(i!=j){
					k = k+abs(pinakas[i][j]);	
				}
			}
			pinakas[i][i] = k+2;
		}
	
	
	// Apothikeusi tou A
	FILE *fpA;
	fpA = fopen("C:\\Users\\akog4\\Desktop\\Spil\\jacobiA.txt", "w");
	if (fpA == NULL){
	    printf("Error opening file A!\n");
	    exit(1);
	}
	for (i=0;i<ni;i++){
		for (j=0;j<nj;j++){
			fprintf(fpA, "%d ",pinakas[i][j]);
		}
		fprintf(fpA, "\n");
	}
	fclose(fpA);
	
	// Apothikeusi tou B
	FILE *fpB;
	fpB = fopen("C:\\Users\\akog4\\Desktop\\Spil\\jacobiB.txt", "w");
	if (fpB == NULL){
	    printf("Error opening file B!\n");
	    exit(1);
	}
	for (i=0;i<ni;i++){
		fprintf(fpB, "%d ",b[i]);
		fprintf(fpB, "\n");
	}
	fclose(fpB);
}
