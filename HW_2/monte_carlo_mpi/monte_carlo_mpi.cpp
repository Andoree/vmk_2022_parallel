// monte_carlo_mpi.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include "mpi.h"
#include <iostream>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
// TODO: add mpi

float f(double x, double y, double z) {
	// TODO
	// return 1;
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
int main(int argc, char** argv)
{	
	int size, rank;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double start_time, finish_time;
	// printf("Current process: %d, Total amount of processes: %d\n", rank, size);
	const float ANALYTIC_INTEGRAL_VALUE = 1. / (13 * 28); // TODO
	const int BATCH_SIZE = 1024;
	// const MPI_INT BATCH_SIZE = 1024;
	const int DIMENSIONALITY = 3;
	double PARALLELEPIPED_LOWER_BOUNDS[] = { 0, 0, 0 };
	double PARALLELEPIPED_UPPER_BOUNDS[] = { 1, 1, 1 };
	double local_point[DIMENSIONALITY];
	// TODO: Это не значение интеграла, это сумма
	double monte_carlo_sum = 0;
	double monte_carlo_integral_value = 0; // TODO: Это везде или только на нулевом процессе?
	double monte_carlo_batch;
	double step_monte_carlo_sum = 0;
	// TODO: Это считывать с консоли, пока захардкодю
	// double epsilon = 1.0e-9;// TODO atof(argv[0]);
	double epsilon = 0.8e-5;
	
	// double local_points[BATCH_SIZE * DIMENSIONALITY];// TODO: Не надо? Пусть точки обрабатываются потоково
	
	int i, j;
	int overall_num_points = 0;
	int num_batch_processed_points = 0;
	int step_num_points = 0;
	double random_double;


	double parallelepiped_volume = 1;
	for (i = 0; i < DIMENSIONALITY; i++)
		parallelepiped_volume *= PARALLELEPIPED_UPPER_BOUNDS[i] - PARALLELEPIPED_LOWER_BOUNDS[i];

	srand(rank);
	printf("Current process: %d, Total amount of processes: %d\n", rank, size);
	// TODO: В каждом процессе srand() с номером процессора

	int k = 0;
	start_time = MPI_Wtime();
	while (abs(ANALYTIC_INTEGRAL_VALUE - monte_carlo_integral_value) >= epsilon) {
	// for (int m = 0; m < 10;m++) {
		monte_carlo_batch = 0.;
		num_batch_processed_points = 0;
		// Создание и обработка случайных точек
		for (i = 0; i < BATCH_SIZE; i++) {
			for (j = 0; j < DIMENSIONALITY; j++) {
				random_double = (double)rand() / RAND_MAX;
				random_double = PARALLELEPIPED_LOWER_BOUNDS[j] + random_double * (PARALLELEPIPED_UPPER_BOUNDS[j] - PARALLELEPIPED_LOWER_BOUNDS[j]);
				// printf("random_double %f\n", random_double);
				local_point[j] = random_double;

			}
			monte_carlo_batch += g(local_point[0], local_point[1], local_point[2]);
			num_batch_processed_points += 1;
		}
		k += 1;
		MPI_Reduce(&monte_carlo_batch, &step_monte_carlo_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&num_batch_processed_points, &step_num_points, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rank == 0) {
			// printf("ANALYTIC_INTEGRAL_VALUE %f\n", ANALYTIC_INTEGRAL_VALUE);
			// printf("monte_carlo_integral_value %f\n", monte_carlo_integral_value);
			// printf("epsilon %f\n", epsilon);
			// printf("abs(ANALYTIC_INTEGRAL_VALUE - monte_carlo_integral_value) %f\n", abs(ANALYTIC_INTEGRAL_VALUE - monte_carlo_integral_value));
			// printf("abs(ANALYTIC_INTEGRAL_VALUE - monte_carlo_integral_value) >= epsilon %i\n---\n", abs(ANALYTIC_INTEGRAL_VALUE - monte_carlo_integral_value) >= epsilon);
			// printf("monte_carlo_step_sum %f\n", step_monte_carlo_sum);
			// printf("step_num_points %f\n", step_num_points);
			// printf("monte_carlo_integral_value %f\n", monte_carlo_integral_value);
			//printf("overall_num_points %i\n", overall_num_points);
			// printf("num_batch_processed_point %i\n", num_batch_processed_points);
			// TODO: Ошибка в том, что я затираю старые результаты. Надо переделать REDUCE.
			// printf("Step %i. Monte Carlo integral value: %f %f %f\n", k, monte_carlo_integral_value, ANALYTIC_INTEGRAL_VALUE, epsilon);
			overall_num_points += step_num_points;
			monte_carlo_sum += step_monte_carlo_sum;

			monte_carlo_integral_value = monte_carlo_sum * parallelepiped_volume / overall_num_points;
			

			// printf("Current process: %d, Total amount of processes: %d\n", nPr, size);
			
		}
		MPI_Bcast(&monte_carlo_integral_value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		printf("Step %i, rank %i. Monte Carlo integral value: %.8f %.8f %.8f\n", k, rank, monte_carlo_integral_value, ANALYTIC_INTEGRAL_VALUE, epsilon);
		
	}
	printf("--\n\nANALYTIC_INTEGRAL_VALUE %f\n", ANALYTIC_INTEGRAL_VALUE);
	printf("monte_carlo_integral_value %f\n", monte_carlo_integral_value);
	printf("epsilon %f\n", epsilon);
	printf("abs(ANALYTIC_INTEGRAL_VALUE - monte_carlo_integral_value) %f\n", abs(ANALYTIC_INTEGRAL_VALUE - monte_carlo_integral_value));
	printf("abs(ANALYTIC_INTEGRAL_VALUE - monte_carlo_integral_value) >= epsilon %i\n---\n", abs(ANALYTIC_INTEGRAL_VALUE - monte_carlo_integral_value) >= epsilon);


	finish_time = MPI_Wtime();
	if (rank == 0)
		printf("Execution time: %f\n", (finish_time - start_time) );
	
	
	
	MPI_Finalize();
	return 0;
}
