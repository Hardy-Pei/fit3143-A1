//////////////////////////////////////////////////////////////////////////////////////
// Student Name: Zhihao Pei
// Student ID: 28294335
// Student Email: zpei0001@student.monash.edu
//////////////////////////////////////////////////////////////////////////////////////

// Main program
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

// Main program
int main(int argc, char* argv[])
 {
	/* screen ( integer) coordinate */
	int iX,iY;
	const int iXmax = 8000; // default
	const int iYmax = 8000; // default

	/* world ( double) coordinate = parameter plane*/
	double Cx, Cy;
	const double CxMin = -2.5;
	const double CxMax = 1.5;
	const double CyMin = -2.0;
	const double CyMax = 2.0;

	/* */
	double PixelWidth = (CxMax - CxMin)/iXmax;
	double PixelHeight = (CyMax - CyMin)/iYmax;

	/* color component ( R or G or B) is coded from 0 to 255 */
	/* it is 24 bit color RGB file */
	const int MaxColorComponentValue = 255; 
	FILE * fp;
	char *filename = "Mandelbrot.ppm";
	char *comment = "# ";	/* comment should start with # */

	// RGB color array
	static unsigned char *colors = NULL;

	/* Z = Zx + Zy*i;	Z0 = 0 */
	double Zx, Zy;
	double Zx2, Zy2; /* Zx2 = Zx*Zx;  Zy2 = Zy*Zy  */
	/*  */
	int Iteration;
	const int IterationMax = 2000; // default

	/* bail-out value , radius of circle ;  */
	const double EscapeRadius = 400;
	double ER2 = EscapeRadius * EscapeRadius;

	// rank = index of the processor, size = number of the processor
	int rank, size;

	// processor 0 becomes the master processor
	int root = 0;

	// how many rows assigned to each processor
	int rows_per_procs, row_remain;
	int current, last, index;
	int i, j, k, m;
	MPI_Init(&argc, &argv);
	MPI_Status stat;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	/* Clock information */
	double startTime, endTime, startComp, startComm, startWrite;

	
	rows_per_procs = iYmax / size;
	row_remain = iYmax % size;
	/* compute and write image data bytes to the file */
	if (rank == root) {
		// Get current clock time as the start time
		startTime = MPI_Wtime();
		/*create new file,give it a name and open it in binary mode  */
		fp = fopen(filename, "wb"); /* b -  binary mode */
		/*write ASCII header to the file (PPM file format)*/
		fprintf(fp,"P6\n %s\n %d\n %d\n %d\n", comment, iXmax, iYmax, MaxColorComponentValue);
		printf("File: %s successfully opened for writing.\n", filename);
		printf("Computing Mandelbrot Set. Please wait...\n");
		// create a dynamic array to hold the results
		colors = (unsigned char*) malloc((iYmax*iXmax*3)*sizeof(unsigned char));
		// Get current time as the start of the computation
		startComp = MPI_Wtime();
	} else if (rank < row_remain) {
		colors = (unsigned char*) malloc(((rows_per_procs+1)*iXmax*3)*sizeof(unsigned char));
	} else {
		colors = (unsigned char*) malloc((rows_per_procs*iXmax*3)*sizeof(unsigned char));
	}

	for(iY = rank; iY < iYmax; iY+=size)
	{
		Cy = CyMin + (iY * PixelHeight);
		if (fabs(Cy) < (PixelHeight / 2))
		{
			Cy = 0.0; /* Main antenna */
		}
		current = iY/size*iXmax*3;
		for(iX = 0; iX < iXmax; iX++)
		{
			Cx = CxMin + (iX * PixelWidth);
			/* initial value of orbit = critical point Z= 0 */
			Zx = 0.0;
			Zy = 0.0;
			Zx2 = Zx * Zx;
			Zy2 = Zy * Zy;

			/* */
			for(Iteration = 0; Iteration < IterationMax && ((Zx2 + Zy2) < ER2); Iteration++)
			{
				Zy = (2 * Zx * Zy) + Cy;
				Zx = Zx2 - Zy2 + Cx;
				Zx2 = Zx * Zx;
				Zy2 = Zy * Zy;
			};
			index = current+iX*3;
			/* compute  pixel color (24 bit = 3 bytes) */
			if (Iteration == IterationMax)
			{
				// Point within the set. Mark it as black
				colors[index] = 0;
				colors[index+1] = 0;
				colors[index+2] = 0;
			}
			else 
			{
				// Point outside the set. Mark it as white
				double c = 3*log((double)Iteration)/log((double)(IterationMax) - 1.0);
				if (c < 1)
				{
					colors[index] = 0;
					colors[index+1] = 0;
					colors[index+2] = 255*c;
				}
				else if (c < 2)
				{
					colors[index] = 0;
					colors[index+1] = 255*(c-1);
					colors[index+2] = 255;
				}
				else
				{
					colors[index] = 255*(c-2);
					colors[index+1] = 255;
					colors[index+2] = 255;
				}
			}
		}
	}
	if (rank == root) {
		// Get the clock current time as the start of the communication
		startComm = MPI_Wtime();
		current += iXmax*3;
		// master processor gathers data from all other processors
		for (i=1; i<size; i++) {
			if (i<row_remain) {
				MPI_Recv(colors+current, (rows_per_procs+1)*iXmax*3, MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD, &stat);
				current += (rows_per_procs+1)*iXmax*3;
			} else {
				MPI_Recv(colors+current, rows_per_procs*iXmax*3, MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD, &stat);
				current += rows_per_procs*iXmax*3;
			}
		}
		// Get the clock current time as the start of the writing
		startWrite = MPI_Wtime();
		current = 0;
		// retrieve data from the array and write into the file
		for (j=root; j<rows_per_procs+1; j++) {
			if (j == rows_per_procs) {
				last = row_remain;
			} else {
				last = size;
			}
			for (k=root; k<last; k++) {
				for (m=0;m<iXmax;m++) {
					fwrite(colors+(current+m*3), 1, 3, fp);
				}
				if (k < row_remain) {
					current += ((rows_per_procs+1)*iXmax*3);
				} else if (k == size-1) {
					current = (j+1)*3*iXmax;
				} else {
					current += (rows_per_procs*iXmax*3);
				}
			}
		}
		// Get the current time as the end of the program
		endTime = MPI_Wtime();
		printf("Completed Computing Mandelbrot Set.\n");
		printf("File: %s successfully closed.\n", filename);
		// overall time
		printf("Mandelbrot computational process time: %lf\n", endTime-startTime);
		// initiation time
		printf("Initiation process time: %lf\n", startComp-startTime);
		// parallel computing time
		printf("Parallel computing process time: %lf\n", startComm-startComp);
		// Communication time
		printf("Communication time: %lf\n", startWrite-startComm);
		// writing time
		printf("Writing process time: %lf\n", endTime-startWrite);
		fclose(fp);
		// other processors send the data back to the master processor
	} else if (rank < row_remain) {
		MPI_Send(colors, (rows_per_procs+1)*iXmax*3, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
	} else {
		MPI_Send(colors, rows_per_procs*iXmax*3, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
	}
	free(colors);
	MPI_Finalize();
	return 0;
 }