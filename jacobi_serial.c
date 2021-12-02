#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <time.h>

int main(){
    int n,i,j,l=0;
	printf("enter the size of equation: ");
    scanf("%d",&n);
    double a[n][n],b[n][1],x[n][1],T[n][1],e,k;
    


clock_t begin = clock();//arxizei na metra xrono
    // Diavazei ton pinaka A
    FILE *fpA;
	fpA = fopen( "inA.txt", "r" );
	for (i=0; i<n; i++) {
	    for (j=0; j<n; j++) {
		fscanf( fpA, "%lf", &a[i][j] );
	    }
	    fscanf( fpA, "\n" );
	}

    // Diavazei ton pinaka B
	FILE *fpB;
	fpB =fopen( "inB.txt", "r" );
	for (i=0; i<n; i++) {
		fscanf( fpB, "%lf", &b[i] );
		fscanf( fpB, "\n" );
	}
	
	//Diavazei tin akriveia
	printf("\nEnter the Accuracy = ");
	scanf("%lf",&e);

    for (i=0;i<n;i++)
        T[i][0]=0;
    while (l!=n)
    {
        l=0;
        for (i=0;i<n;i++)
        {
            x[i][0]=(1/a[i][i])*(b[i][0]);
            for (j=0;j<n;j++)
            {
                if (j!=i)
                x[i][0]=x[i][0]-(1/a[i][i])*(a[i][j]*T[j][0]);
            }
        }
        for(i=0;i<n;i++)
        {
            k=fabs(x[i][0]-T[i][0]);
            if (k<=e)
            {
                l=l+1;
            }
        }
    for (i=0;i<n;i++)
        T[i][0]=x[i][0];
    }
    for (i=0;i<n;i++)
       printf("x%d =%lf \n",i+1,x[i][0]);
   
  clock_t end = clock();
double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;//xronos p ekane synolika gia na tre3ei
printf("%lf \n",time_spent);
    return 0;
}

