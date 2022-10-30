

#include "mpi.h"
#include <iostream>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <cmath>

float f(double x, double y, double z) {

	return x * y * y * z * z * z;
}



float g(double x, double y, double z) {
	if (!(z >= 0 && z <= x * y))
		return 0;
	if (!(x < 1 && x > 0))
		return 0;
	if (y > x)
		return 0;
	return f(x, y, z);
}
int main(int argc, char* argv[])
{

	int size, rank;
	int x, y;
	double epsilon;
	MPI_Status status;

	epsilon = atof(argv[1]);
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Bcast(&epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	double start_time, finish_time;

	const float ANALYTIC_INTEGRAL_VALUE = 1. / (13 * 28);
	const int BATCH_SIZE = 16;
	const int DIMENSIONALITY = 3;
	const int PATIENCE = 10;
	double PARALLELEPIPED_LOWER_BOUNDS[] = { 0, 0, 0 };
	double PARALLELEPIPED_UPPER_BOUNDS[] = { 2, 2, 2 };
	double local_point[DIMENSIONALITY];

	double monte_carlo_sum = 0;
	double monte_carlo_integral_value = 0;
	double monte_carlo_batch;
	double step_monte_carlo_sum = 0;

	int i, j;
	long overall_num_points = 0;
	int num_batch_processed_points = 0;
	int step_num_points = 0;
	int patience = 0;
	double random_double;


	double parallelepiped_volume = 1;
	for (i = 0; i < DIMENSIONALITY; i++)
		parallelepiped_volume *= PARALLELEPIPED_UPPER_BOUNDS[i] - PARALLELEPIPED_LOWER_BOUNDS[i];

	srand(rank);

	int k = 0;
	start_time = MPI_Wtime();
	while (fabs(ANALYTIC_INTEGRAL_VALUE - monte_carlo_integral_value) >= epsilon && patience < PATIENCE) {
		if (fabs(ANALYTIC_INTEGRAL_VALUE - monte_carlo_integral_value) < epsilon) {
			patience += 1;
		}
		else {
			patience = 0;
		}

		monte_carlo_batch = 0.;
		num_batch_processed_points = 0;
		for (i = 0; i < BATCH_SIZE; i++) {
			for (j = 0; j < DIMENSIONALITY; j++) {
				random_double = (double)rand() / RAND_MAX;
				random_double = PARALLELEPIPED_LOWER_BOUNDS[j] + random_double * (PARALLELEPIPED_UPPER_BOUNDS[j] - PARALLELEPIPED_LOWER_BOUNDS[j]);
				local_point[j] = random_double;

			}
			monte_carlo_batch += g(local_point[0], local_point[1], local_point[2]);
			num_batch_processed_points += 1;
		}
		k += 1;
		MPI_Reduce(&monte_carlo_batch, &step_monte_carlo_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&num_batch_processed_points, &step_num_points, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			overall_num_points += step_num_points;
			monte_carlo_sum += step_monte_carlo_sum;

			monte_carlo_integral_value = monte_carlo_sum * parallelepiped_volume / overall_num_points;

		}
		MPI_Bcast(&monte_carlo_integral_value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	}

	finish_time = MPI_Wtime();
	if (rank == 0) {
		double execution_time = finish_time - start_time;
		double error = fabs(ANALYTIC_INTEGRAL_VALUE - monte_carlo_integral_value);
		printf("%.12f %.16f %ld %.12f\n", monte_carlo_integral_value, error, overall_num_points, execution_time);
	}


	MPI_Finalize();
	return 0;
}
