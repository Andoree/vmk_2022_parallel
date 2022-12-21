
#define _USE_MATH_DEFINES

#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <time.h>

// Вариант 4
// u_4 = 1 + cos(pi * x * y)
// k_3 = 4 + x + y
// q_0 = 0
// gamma_R - тип 3
// gamma_L - тип 1  
// gamma_T - тип 3  
// gamma_B - тип 1  
// П = [0, 2] x [0, 1]

// Функция u(x, y), вариант 4, функция u_4
double u_4(double x, double y) {
    return 1.0 + cos(M_PI * x * y);
}

// Функция k(x, y), вариант 4, функция k_3
double k_3(double x, double y) {
    return 4.0 + x + y;
}

// Функция q(x, y), вариант 4, функция q_0
double q_0(double x, double y) {
    return 0.;
}

// Выведенная аналитически функция F(x, y)
double F(double x, double y) {
    double dy = M_PI * x * (sin(M_PI * x * y) + k_3(x, y) * M_PI * x * cos(M_PI * x * y));
    double dx = M_PI * y * (sin(M_PI * x * y) + k_3(x, y) * M_PI * y * cos(M_PI * x * y));
    return dx + dy;
}

double phi(double x, double y) {
    return u_4(x, y);
}

double psi_right_border(double x, double y) {
    return u_4(x, y) - (k_3(x, y) * M_PI * y * sin(M_PI * x * y));
}

double psi_top_border(double x, double y) {
    return u_4(x, y) - (k_3(x, y) * M_PI * x * sin(M_PI * x * y));
}

double right_w_ij_difference(double** w, int i, int j, int M, double h_1) {
    return (w[i + 1][j] - w[i][j]) / h_1;
}

double left_w_ij_difference(double** w, int i, int j, int M, double h_1) {
    return -(w[i][j] - w[i - 1][j]) / h_1;
}

double top_w_ij_difference(double** w, int i, int j, int M, double h_2) {
    return (w[i][j + 1] - w[i][j]) / h_2;
}

double bottom_w_ij_difference(double** w, int i, int j, int M, double h_2) {
    return -(w[i][j] - w[i][j - 1]) / h_2;
}

double delta_h_w_ij_x_right_part(double** w, int i, int j, int M, double h_1, double h_2) {
    return (k_3(i * h_1 + 0.5 * h_1, j * h_2) * right_w_ij_difference(w, i, j, M, h_1));
}

double delta_h_w_ij_x_left_part(double** w, int i, int j, int M, double h_1, double h_2) {
    return  (k_3(i * h_1 - 0.5 * h_1, j * h_2) * left_w_ij_difference(w, i, j, M, h_1));
}

double delta_h_w_ij_y_top_part(double** w, int i, int j, int M, double h_1, double h_2) {
    return (k_3(i * h_1, j * h_2 + 0.5 * h_2) * top_w_ij_difference(w, i, j, M, h_2));
}

double delta_h_w_ij_y_bottom_part(double** w, int i, int j, int M, double h_1, double h_2) {
    return (k_3(i * h_1, j * h_2 - 0.5 * h_2) * bottom_w_ij_difference(w, i, j, M, h_2));
}

double delta_h_w_ij_x_part(double** w, int i, int j, int M, double h_1, double h_2) {
    double right_part = (k_3(i * h_1 + 0.5 * h_1, j * h_2) * right_w_ij_difference(w, i, j, M, h_1));
    double left_part = (k_3(i * h_1 - 0.5 * h_1, j * h_2) * left_w_ij_difference(w, i, j, M, h_1));

    return (1. / h_1) * (right_part + left_part);
}

double delta_h_w_ij_y_part(double** w, int i, int j, int M, double h_1, double h_2) {
    double top_part = (k_3(i * h_1, j * h_2 + 0.5 * h_2) * top_w_ij_difference(w, i, j, M, h_2));
    double bottom_part = (k_3(i * h_1, j * h_2 - 0.5 * h_2) * bottom_w_ij_difference(w, i, j, M, h_2));

    return (1. / h_2) * (top_part + bottom_part);
}

double delta_h_w_ij(double** w, int i, int j, int M, double h_1, double h_2) {
    double x_part = delta_h_w_ij_x_part(w, i, j, M, h_1, h_2);
    double y_part = delta_h_w_ij_y_part(w, i, j, M, h_1, h_2);

    return  x_part + y_part;
}

void equation_right_part(double** B, int M, int N, double h_1, double h_2) {
    // Обработка внутренних точек сетки
    for (int i = 1; i < M; i++) {
        for (int j = 1; j < N; j++) {
            B[i][j] = F(i * h_1, j * h_2);
        }
    }

    // Нижние и верхние границы
    for (int i = 1; i < M; i++) {
        // Нижняя граница - Тип 1 (phi = u)
        B[i][0] = phi(i * h_1, 0.);
        // Верхняя граница - Тип 3 
        B[i][N] = F(i * h_1, N * h_2) + (2. / h_2) * psi_top_border(i * h_1, N * h_2);
    }
    // Левая и правая границы без угловых точек.
    for (int j = 1; j < N; j++) {
        // Левая  граница - Тип 1 (phi = u)
        B[0][j] = u_4(0., h_2 * j);
        // Правая граница - Тип 3
        B[M][j] = F(M * h_1, j * h_2) + (2. / h_1) * psi_right_border(M * h_1, j * h_2);

    }

    // Левый нижний угол 
    B[0][0] = u_4(0., 0.);
    // Правый нижний угол
    B[M][0] = u_4(M * h_1, 0.);
    // Левый верхний угол
    B[0][N] = u_4(0., N * h_2);
    // Тип 3 + Тип 3, верхний правый угол
    B[M][N] = F(M * h_1, N * h_2) + (2. / h_1) * psi_right_border(M * h_1, N * h_2) 
    + (2. / h_2) * psi_top_border(M * h_1, N * h_2);
    
}


void equation_left_part(double** w, double** r, double h_1, double h_2, int M, int N) {
    // Обработка внутренних точек сетки
    for (int i = 1; i < M; i++) {
        for (int j = 1; j < N; j++) {
            r[i][j] = -delta_h_w_ij(w, i, j, M, h_1, h_2) + q_0(i * h_1, j * h_2) * w[i][j];
        }
    }

    // Нижние и верхние границы
    for (int i = 1; i < M; i++) {
        // Нижняя граница - Тип 1 (phi = u), j = 0 
        r[i][0] = u_4(h_1 * i, 0.);
        // Верхняя граница - Тип 3
        r[i][N] = -delta_h_w_ij_x_part(w, i, N, M, h_1, h_2) -
            +(2. / h_2) * delta_h_w_ij_y_bottom_part(w, i, N, M, h_1, h_2)
            + w[i][N] * ((2. / h_2));
    }
    // Левая и правая границы без угловых точек.
    for (int j = 1; j < N; j++) {
        // Левая  граница - Тип 1 (phi = u), i = 0
        r[0][j] = u_4(0., h_2 * j);
        // Правая граница - Тип 3
        r[M][j] = -delta_h_w_ij_y_part(w, M, j, M, h_1, h_2)
           - (2. / h_1) * delta_h_w_ij_x_left_part(w, M, j, M, h_1, h_2)
            + w[M][j] * ((2. / h_1));
    }

    // Левый нижний угол 
    r[0][0] = u_4(0., 0.);
    // Правый нижний угол
    r[M][0] = u_4(M * h_1, 0.);
    // Левый верхний угол
    r[0][N] = u_4(0., N * h_2);
    // Правый верхний угол - тип 3 + тип 3
    r[M][N] = -(2. / h_1) * delta_h_w_ij_x_left_part(w, M, N, M, h_1, h_2)
     -(2. / h_2) * delta_h_w_ij_y_bottom_part(w, M, N, M, h_1, h_2)
       + ((2. / h_1) + (2. / h_2)) * w[M][N];
    
    
}

double rho_x_i(int i, int M) {
    if (i >= 1 && i <= M - 1)
        return 1.;
    return 0.5;

}

double rho_y_j(int i, int N) {
    if (i >= 1 && i <= N - 1)
        return 1.;
    return 0.5;

}

// Коэффициент rho, нужный для вычисления скалярного произведения
double rho(int i, int j, int M, int N) {
    return rho_x_i(i, M) * rho_y_j(j, N);
}

// Скалярное произведение
double dot_product(double** u, double** v, double h_1, double h_2, int M, int N) {
    double dot_p = 0.;
    double inner_sum;
    for (int i = 0; i <= M; i++) {
        inner_sum = 0.0;
        for (int j = 0; j <= N; j++) {
            inner_sum += h_2 * rho(i, j, M, N) * u[i][j] * v[i][j];
        }
        dot_p += h_1 * inner_sum;
    }
    return dot_p;
}

// Евклидова норма
double euclidean_norm(double** u, double h_1, double h_2, int M, int N) {
    return sqrt(dot_product(u, u, h_1, h_2, M, N));
}

// Функция инициализации сеточной функции 
void init_grid(double** w, double h_1, double h_2, int M, int N) {
    for (int i = 0; i <= M; i++)
        for (int j = 0; j <= N; j++) {
            if (i == 0 || j == 0) {
                w[i][j] = u_4(h_1 * i, h_2 * j);
            }
            else {
                w[i][j] = 1.;
            }
        }
}


int main(int argc, char** argv) {

    int M = atoi(argv[1]);
    int N = atoi(argv[2]);
    // int M = 120;
    // int N = 120;
    

    const double EPS = 1e-6;
    const double h_1 = 2. / (double)M;
    const double h_2 = 1. / (double)N;
    double tau = 0.;
    double last_solution_update_norm;
    time_t begin, end;

    double** B = (double**)malloc((M + 1) * sizeof(double*));
    for (int i = 0; i <= M; i++)
        B[i] = (double*)malloc((N + 1) * sizeof(double));
    double** w = (double**)malloc((M + 1) * sizeof(double*));
    for (int i = 0; i <= M; i++)
        w[i] = (double*)malloc((N + 1) * sizeof(double));
    double** r = (double**)malloc((M + 1) * sizeof(double*));
    for (int i = 0; i <= M; i++)
        r[i] = (double*)malloc((N + 1) * sizeof(double));
    double** current_w = (double**)malloc((M + 1) * sizeof(double*));
    for (int i = 0; i <= M; i++)
        current_w[i] = (double*)malloc((N + 1) * sizeof(double));
    double** real_u = (double**)malloc((M + 1) * sizeof(double*));
    for (int i = 0; i <= M; i++)
        real_u[i] = (double*)malloc((N + 1) * sizeof(double));
    double** error_matrix = (double**)malloc((M + 1) * sizeof(double*));
    for (int i = 0; i <= M; i++)
        error_matrix[i] = (double*)malloc((N + 1) * sizeof(double));
    double** Ar = (double**)malloc((M + 1) * sizeof(double*));
    for (int i = 0; i <= M; i++)
        Ar[i] = (double*)malloc((N + 1) * sizeof(double));
    
    
    
    // Инициализация сетки
    init_grid(w, h_1, h_2, M, N);

    time(&begin);

    // Инициализация правой части уравнения
    equation_right_part(B, M, N, h_1, h_2);
    
    bool keep_calculation = true;
    while (keep_calculation) {
        // Вычисление невязки r
        equation_left_part(w, r, h_1, h_2, M, N);
        for (int i = 0; i <= M; i++)
            for (int j = 0; j <= N; j++) {
                r[i][j] -= B[i][j];
                current_w[i][j] = w[i][j];
            }
        // Вычисление Ar
        equation_left_part(r, Ar, h_1, h_2, M, N);
        // Вычисление tau
        tau = dot_product(Ar, r, h_1, h_2, M, N) / pow(euclidean_norm(Ar, h_1, h_2, M, N), 2);
        // Вычисление следующего значения w
        for (int i = 0; i <= M; i++)
            for (int j = 0; j <= N; j++)
                w[i][j] = w[i][j] - tau * r[i][j];
        for (int i = 0; i <= M; i++)
            for (int j = 0; j <= N; j++) {
                // printf("%lf,%lf\n", w[i][j], current_w[i][j]);
                current_w[i][j] = w[i][j] - current_w[i][j];
            }
        last_solution_update_norm = euclidean_norm(current_w, h_1, h_2, M, N);

        // printf("last_solution_update_norm %.10f\n", last_solution_update_norm);
        
        if (last_solution_update_norm < EPS) {
            keep_calculation = false;
        }

    }
    time(&end);
    time_t elapsed = end - begin;

    for (int i = 0; i <= M; i++) {
        for (int j = 0; j <= N; j++) {
            real_u[i][j] = u_4(i * h_1, j * h_2);
        }
    }
    /*
    printf("Real U\n");
    for (int i = 0; i <= M; i++) {
        for (int j = 0; j <= N; j++) {
            printf("%.7f ", real_u[i][j]);
        }
        printf("\n");
    }
    printf("Approximate solution W\n");
    for (int i = 0; i <= M; i++) {
        for (int j = 0; j <= N; j++) {
            printf("%.7f ", w[i][j]);
        }
        printf("\n");
    }
    */

    for (int i = 0; i <= M; i++)
        for (int j = 0; j <= N; j++)
            error_matrix[i][j] = real_u[i][j] - w[i][j];
    double error_scalar = euclidean_norm(error_matrix, h_1, h_2, M, N);

    for (int i = 0; i <= M; i++)
         for (int j = 0; j <= N; j++)
            printf("%.10f,%.10f\n", w[i][j], real_u[i][j]);
    printf("Execution time: %ld seconds.\n", elapsed);
    printf("Error value: %.10f\n", error_scalar);
    return 0;
}
