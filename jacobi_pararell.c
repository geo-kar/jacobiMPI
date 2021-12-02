#include <stdio.h>
#include "mpi.h"
#include <math.h>
#include <time.h>

#define Swap(x,y) {float* temp; temp = x; x = y; y = temp;}

#define MAX_DIM 2048
#define file_size 2048

typedef  float MATRIX_T[MAX_DIM][MAX_DIM];


int Parallel_jacobi(
        MATRIX_T  A_local    /* in  */, 
        float     x_local[]  /* out */, 
        float     b_local[]  /* in  */, 
        int       n          /* in  */, 
        float     tol        /* in  */, 
        int       max_iter   /* in  */,
        int       p          /* in  */, 
        int       my_rank    /* in  */);

void Read_matrix(char* prompt, MATRIX_T A_local, int n,
         int my_rank, int p);
void Read_vector(char* prompt, float x_local[], int n, int my_rank,
         int p);
void Print_matrix(char* title, MATRIX_T A_local, int n, 
         int my_rank, int p);
void Print_vector(char* title, float x_local[], int n, int my_rank,
         int p);

main(int argc, char* argv[]) {
    int        p;
    int        my_rank;
    MATRIX_T   A_local;
    float      x_local[MAX_DIM];
    float      b_local[MAX_DIM];
    int        n;
    float      tol;
    int        max_iter;
    int        converged;
	
	
	
double t1, t2; 


    MPI_Init(&argc, &argv);//initialize mpi
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  
  
  //eisagoume to megethos t pikana A to sfalma p theloume kai to megisto arithmo epanalipsewn
    if (my_rank == 0) {
        printf("Enter size of equation, tolerance, and max number of iterations\n");
        scanf("%d %f %d", &n, &tol, &max_iter);
    }
//stelnoyme tous arithmous sto mpi me int
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tol, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_iter, 1, MPI_INT, 0, MPI_COMM_WORLD);
//diabazoume ton pinaka A kai B
    Read_matrix("Matrix file(inA.txt) will br inserter automatically", A_local, n, my_rank, p);
    Read_vector("Matrix file(inB.txt) will br inserter automatically", b_local, n, my_rank, p);
t1 = MPI_Wtime();//metrisi xronou me clock.h kai mpi_time
clock_t begin = clock();

//kaloume tin synarthsh p efarmozei to algoritho 
// ean mporei na lythei o pinaka kanoyme print ti lisi
    converged = Parallel_jacobi(A_local, x_local, b_local, n,
        tol, max_iter, p, my_rank);

    if (converged)
        Print_vector("The solution is", x_local, n, my_rank, p);
    else
        if (my_rank == 0)
            printf("Failed to converge in %d iterations\n", max_iter);


	
	
	
	t2 = MPI_Wtime(); 
printf( "Elapsed time is %f\n", t2 - t1 ); 
	
	
	   MPI_Finalize();
	
	 clock_t end = clock();//metrhsh xronou teliki
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;//xronos p ekane synolika gia na tre3ei
printf("O xronos pou ekane: %4.20lf \n",time_spent);
}  /* main */


/*********************************************************************/
/* Return 1 if iteration converged, 0 otherwise */
/* MATRIX_T is a 2-dimensional array            */
int Parallel_jacobi(
        MATRIX_T  A_local    /* in  */, 
        float     x_local[]  /* out */, 
        float     b_local[]  /* in  */, 
        int       n          /* in  */, 
        float     tol        /* in  */, 
        int       max_iter   /* in  */,
        int       p          /* in  */, 
        int       my_rank    /* in  */) {
    int     i_local, i_global, j;
    int     n_bar;
    int     iter_num;
    float   x_temp1[MAX_DIM];
    float   x_temp2[MAX_DIM];
    float*  x_old;
    float*  x_new;

    float Distance(float x[], float y[], int n);
//xwrizoume ton pinaka X se isa kommatia analoga me ton arithmo twn ypologistwn
    n_bar = n/p;
//pernoume ton pinaka B    
    /* Initialize x */
    MPI_Allgather(b_local, n_bar, MPI_FLOAT, x_temp1,
        n_bar, MPI_FLOAT, MPI_COMM_WORLD);
    x_new = x_temp1;
    x_old = x_temp2;
//ftiaxnoyme 2 vectors lisewn, to palio X k kainoyrgio X 
    iter_num = 0;
    do {
        iter_num++;
        
        /* allazoume tis times t x_old kai x_new, ekteloume ton algorithmo dld */
        Swap(x_old, x_new);
        for (i_local = 0; i_local < n_bar; i_local++){
            i_global = i_local + my_rank*n_bar;
            x_local[i_local] = b_local[i_local];
            for (j = 0; j < i_global; j++)
                x_local[i_local] = x_local[i_local] -  
                    A_local[i_local][j]*x_old[j];
            for (j = i_global+1; j < n; j++)
                x_local[i_local] = x_local[i_local] -   
                    A_local[i_local][j]*x_old[j];
            x_local[i_local] = x_local[i_local]/
                    A_local[i_local][i_global];
        }
//pernoyme ton pinaka X me tis kainoyrgies liseis
        MPI_Allgather(x_local, n_bar, MPI_FLOAT, x_new,
            n_bar, MPI_FLOAT, MPI_COMM_WORLD);
    } while ((iter_num < max_iter) && 
             (Distance(x_new,x_old,n) >= tol));

    if (Distance(x_new,x_old,n) < tol)
        return 1;
    else
        return 0;
} /* Jacobi */


/*********************************************************************/
float Distance(float x[], float y[], int n) {
    int i;
    float sum = 0.0;

    for (i = 0; i < n; i++) {
        sum = sum + (x[i] - y[i])*(x[i] - y[i]);
    }
    return sqrt(sum);
} /* Distance */


/*********************************************************************/
void Read_matrix(
         char*     prompt   /* in  */,
         MATRIX_T  A_local  /* out */,
         int       n        /* in  */,
         int       my_rank  /* in  */,
         int       p        /* in  */) {

    int       i, j;
    MATRIX_T  temp;
    int       n_bar;
 
    n_bar = n/p;


    if (my_rank == 0) {
		
		
		FILE *fp;
	fp =fopen( "inA.txt", "r" );
	if (!fp) MPI_Abort( MPI_COMM_WORLD, 1 );
	/* This includes the top and bottom edge */
	for (i=0; i<file_size; i++) {
	    for (j=0; j<file_size; j++) {
		fscanf( fp, "%f", &temp[i][j] );
	    }
	    fscanf( fp, "\n" );
	}
		
		
		
		
		
		
		
		
	/*	
        printf("%s\n", prompt);
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                scanf("%f",&temp[i][j]);
			*/
    }    
  //stenoume sto mpi ton pinaka A 
    MPI_Scatter(temp, n_bar*MAX_DIM, MPI_FLOAT, A_local,
        n_bar*MAX_DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);


}  /* Read_matrix */

/*********************************************************************/
void Read_vector(
         char*  prompt     /* in  */,
         float  x_local[]  /* out */,
         int    n          /* in  */,
         int    my_rank    /* in  */,
         int    p          /* in  */) {

    int   i;
    float temp[MAX_DIM];
    int   n_bar;
    
    n_bar = n/p;
	
	FILE *fp;
	fp =fopen( "inB.txt", "r" );
	if (!fp) MPI_Abort( MPI_COMM_WORLD, 1 );
	/* This includes the top and bottom edge */
	for (i=0; i<file_size; i++) {
	    
		fscanf( fp, "%f", &temp[i] );
	  fscanf( fp, "\n" );
	}
    
	
	
/*
    if (my_rank == 0) {
        printf("%s\n", prompt);
        for (i = 0; i < n; i++)
            scanf("%f", &temp[i]);
    }
*/	
//stelnoume to pinaka B
    MPI_Scatter(temp, n_bar, MPI_FLOAT, x_local, n_bar, MPI_FLOAT,
        0, MPI_COMM_WORLD);

}  /* Read_vector */


/*********************************************************************/
void Print_matrix(char* title, MATRIX_T A_local, int n, 
         int my_rank, int p);
void Print_matrix(
         char*     title      /* in */,
         MATRIX_T  A_local    /* in */,
         int       n          /* in */,
         int       my_rank    /* in */,
         int       p          /* in */) {

    int       i, j;
    MATRIX_T  temp;
    int       n_bar;

    n_bar = n/p;

    MPI_Gather(A_local, n_bar*MAX_DIM, MPI_FLOAT, temp,
         n_bar*MAX_DIM, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        printf("%s\n", title);
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++)
                printf("%4.1f ", temp[i][j]);
            printf("\n");
        }
    }
}  /* Print_matrix */


/*********************************************************************/
void Print_vector(char* title, float x_local[], int n, int my_rank,
         int p);
void Print_vector(
         char*  title      /* in */,
         float  x_local[]  /* in */,
         int    n          /* in */,
         int    my_rank    /* in */,
         int    p          /* in */) {

    int   i;
    float temp[MAX_DIM];
    int   n_bar;

    n_bar = n/p;

    MPI_Gather(x_local, n_bar, MPI_FLOAT, temp, n_bar, MPI_FLOAT,
        0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        printf("%s\n", title);
        for (i = 0; i < n; i++)
            printf("%4.10f ", temp[i]);
        printf("\n");
    }
	
	FILE *fpB;
	fpB = fopen("output.txt", "w");
	if (fpB == NULL){
	    printf("Error opening output file!\n");
	    exit(1);
	}
	for (i=0;i<n;i++){
		fprintf(fpB, "%f ",temp[i]);
		fprintf(fpB, "\n");
	}
	fclose(fpB);
}  /* Print_vector */
