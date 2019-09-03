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
	const int IterationMax = 1000; // default

	/* bail-out value , radius of circle ;  */
	const double EscapeRadius = 400;
	double ER2 = EscapeRadius * EscapeRadius;

	int rank, size;
	int root = 0;
	int rows_per_procs, row_remain;
	int current, last, index;
	int i, j, k, m;
	MPI_Init(&argc, &argv);
	MPI_Status stat;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	/* Clock information */
	double startTime, endTime, parallelEnd;

	// Get current clock time.
	/* compute and write image data bytes to the file */
	rows_per_procs = iYmax / size;
	row_remain = iYmax % size;
	if (rank == root) {
		/*create new file,give it a name and open it in binary mode  */
		fp = fopen(filename, "wb"); /* b -  binary mode */
		/*write ASCII header to the file (PPM file format)*/
		fprintf(fp,"P6\n %s\n %d\n %d\n %d\n", comment, iXmax, iYmax, MaxColorComponentValue);
		printf("File: %s successfully opened for writing.\n", filename);
		printf("Computing Mandelbrot Set. Please wait...\n");
		startTime = MPI_Wtime();
		colors = (unsigned char*) malloc((iYmax*iXmax*3)*sizeof(unsigned char));
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
			// printf("%d,%d,%d\n", current, iX,rank);
			// printf("%d, %d, %d\n", colors[current+iX],colors[current+iX+1], colors[current+iX+2]);
		}
	}
	// for(i=0;i<rows_per_procs*iXmax*3;i++) {
	// 	printf("%d\n", colors[i]);
	// }
	// Get the clock current time again
	// Subtract end from start to get the CPU time used.
	if (rank == root) {
		parallelEnd = MPI_Wtime();
		current += iXmax*3;
		for (i=1; i<size; i++) {
			if (i<row_remain) {
				MPI_Recv(colors+current, (rows_per_procs+1)*iXmax*3, MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD, &stat);
				current += (rows_per_procs+1)*iXmax*3;
			} else {
				MPI_Recv(colors+current, rows_per_procs*iXmax*3, MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD, &stat);
				current += rows_per_procs*iXmax*3;
			}
		}
		endTime = MPI_Wtime();
		current = 0;
		for (j=root; j<rows_per_procs+1; j++) {
			if (j == rows_per_procs) {
				last = row_remain;
			} else {
				last = size;
			}
			for (k=root; k<last; k++) {
				for (m=0;m<iXmax;m++) {
					// printf("%d\n", colors[current+m]);
					fwrite(colors+(current+m*3), 1, 3, fp);
					// printf("%d,%d,%d\n", colors[current+m],colors[current+m+1],colors[current+m+2]);
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
		
		printf("Completed Computing Mandelbrot Set.\n");
		printf("File: %s successfully closed.\n", filename);
		printf("Mandelbrot computational process time: %lf\n", parallelEnd-startTime);
		printf("my time: %lf\n", endTime-parallelEnd);
		fclose(fp);
	} else if (rank < row_remain) {
		MPI_Send(colors, (rows_per_procs+1)*iXmax*3, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
	} else {
		MPI_Send(colors, rows_per_procs*iXmax*3, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
	}
	free(colors);
	MPI_Finalize();
	return 0;
 }