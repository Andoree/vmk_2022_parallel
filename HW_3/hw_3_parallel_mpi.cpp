
#define _USE_MATH_DEFINES

#include "mpi.h"
#include <iostream>
#include <stdio.h>      
#include <stdlib.h>    
#include <cmath>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

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

double right_w_ij_difference(double** w, int i, int j, double h_1) {
    return (w[i + 1][j] - w[i][j]) / h_1;
}

double left_w_ij_difference(double** w, int i, int j,  double h_1) {
    return -(w[i][j] - w[i - 1][j]) / h_1;
}

double top_w_ij_difference(double** w, int i, int j, double h_2) {
    return (w[i][j + 1] - w[i][j]) / h_2;
}

double bottom_w_ij_difference(double** w, int i, int j, double h_2) {
    return -(w[i][j] - w[i][j - 1]) / h_2;
}

double delta_h_w_ij_x_right_part(double** w, int i, int j, double global_x, double global_y, double h_1, double h_2) {
    return (1. / h_1) * (k_3(global_x + 0.5 * h_1, global_y) * right_w_ij_difference(w, i, j, h_1));
}

double delta_h_w_ij_x_left_part(double** w, int i, int j, double global_x, double global_y, double h_1, double h_2) {
    return  (1. / h_1) * (k_3(global_x - 0.5 * h_1, global_y) * left_w_ij_difference(w, i, j, h_1));
}

double delta_h_w_ij_y_top_part(double** w, int i, int j, double global_x, double global_y, double h_1, double h_2) {
    return (1. / h_2) * (k_3(global_x, global_y + 0.5 * h_2) * top_w_ij_difference(w, i, j, h_2));
}

double delta_h_w_ij_y_bottom_part(double** w, int i, int j, double global_x, double global_y, double h_1, double h_2) {
    return (1. / h_2) * (k_3(global_x, global_y - 0.5 * h_2) * bottom_w_ij_difference(w, i, j, h_2));
}


int log2_(int num) {
    if (num <= 0)
        return -1;
    int power = 0;
    while ((num & 1) == 0) {
        ++power;
        num = num >> 1;
    }
    if ((num >> 1) != 0)
        return -1;
    return power;
}

int split(int M, int N, int power) {
    double m = (double)M;
    double n = (double)N;
    int px = 0;
    for (int i = 0; i < power; i++) {
        if (m > n) {
            m /= 2.0;
            ++px;
        }
        else {
            n /= 2.0;
        }
    }
    return px;
}

// Функция инициализации сеточной функции 
void init_grid(double** w, double h_1, double h_2, int local_M, int local_N, int local_x_start, int local_y_start) {
    int global_x, global_y;
    for (int i = 0; i <= local_M; i++) {
        global_x = (local_x_start + i);
        for (int j = 0; j <= local_N; j++) {
            global_y = (local_y_start + j);
            if (global_x == 0 || global_y == 0) {
                w[i][j] = u_4(h_1 * global_x, h_2 * global_y);
            }
            else {
                w[i][j] = 1.;
            }
        }
    }
}

void equation_left_part(double** w, double* top_row_of_bottom_neighbor, double* bottom_row_of_top_neighbor,
    double* left_column_of_right_neighbor, double* right_column_of_left_neighbor, 
    int local_M, int local_N, int local_x_start, int local_y_start,
    int left_neighbor_proc_id, int right_neighbor_proc_id, int top_neighbor_proc_id, int bottom_neighbor_proc_id,
    double** r, double h_1, double h_2) {

    double global_x;
    double global_y;
    
    for (int i = 1; i < local_M; i++) {
        for (int j = 1; j < local_N; j++) {
            global_x = (local_x_start + i) * h_1;
            global_y = (local_y_start + j) * h_2;
            r[i][j] = -(delta_h_w_ij_x_right_part(w, i, j, global_x, global_y, h_1, h_2)
                + delta_h_w_ij_x_left_part(w, i, j, global_x, global_y, h_1, h_2)
                + delta_h_w_ij_y_top_part(w, i, j, global_x, global_y, h_1, h_2)
                + delta_h_w_ij_y_bottom_part(w, i, j, global_x, global_y, h_1, h_2));
        }
    }
    
    for (int i = 1; i < local_M; i++) {
        global_x = (local_x_start + i) * h_1;
        global_y = local_y_start * h_2;

        // Нижняя граница - Тип 1 (phi = u), j = 0 
        if (bottom_neighbor_proc_id == -1) {
            r[i][0] = phi(global_x, global_y);
        }
        else {
      
            r[i][0] = -(delta_h_w_ij_x_right_part(w, i, 0, global_x, global_y, h_1, h_2)
                + delta_h_w_ij_x_left_part(w, i, 0, global_x, global_y, h_1, h_2)
                + delta_h_w_ij_y_top_part(w, i, 0, global_x, global_y, h_1, h_2)
                - (1. / (h_2 * h_2) * k_3(global_x, global_y - h_2 * 0.5) * 
                    (w[i][0] - top_row_of_bottom_neighbor[i])));
            
        }
        global_x = (local_x_start + i) * h_1;
        global_y = (local_y_start + local_N) * h_2;

        // Глобальная верхняя граница, тип 3
        if (top_neighbor_proc_id == -1) {
            r[i][local_N] = -(delta_h_w_ij_x_right_part(w, i, local_N, global_x, global_y, h_1, h_2)
                + delta_h_w_ij_x_left_part(w, i, local_N, global_x, global_y, h_1, h_2)
                + 2. * delta_h_w_ij_y_bottom_part(w, i, local_N, global_x, global_y, h_1, h_2))
                 + w[i][local_N] * ((2. / h_2));
        }
        else {      
            r[i][local_N] = -(delta_h_w_ij_x_right_part(w, i, local_N, global_x, global_y, h_1, h_2)
                + delta_h_w_ij_x_left_part(w, i, local_N, global_x, global_y, h_1, h_2)
                + delta_h_w_ij_y_bottom_part(w, i, local_N, global_x, global_y, h_1, h_2)
                + (1. / (h_2 * h_2) * k_3(global_x, global_y + h_2 * 0.5) *
                    (bottom_row_of_top_neighbor[i] - w[i][local_N])));
        }
    }
    
    for (int j = 1; j < local_N; j++) {
        global_x = local_x_start * h_1;
        global_y = (local_y_start + j) * h_2;
        // Левая  граница - Тип 1 (phi = u), i = 0
        if (left_neighbor_proc_id == -1) {
            r[0][j] = phi(global_x, global_y);
        }
        
        else {

            
            r[0][j] = -(delta_h_w_ij_x_right_part(w, 0, j, global_x, global_y, h_1, h_2)
                + delta_h_w_ij_y_top_part(w, 0, j, global_x, global_y, h_1, h_2)
                + delta_h_w_ij_y_bottom_part(w, 0, j, global_x, global_y, h_1, h_2)
                - (1. / (h_1 * h_1) * k_3(global_x - h_1 * 0.5, global_y) *
                    (w[0][j] - right_column_of_left_neighbor[j])));
        }
        global_x = (local_x_start + local_M) * h_1;
        global_y = (local_y_start + j) * h_2;
        // Правая граница - тип 3
        if (right_neighbor_proc_id == -1) {
            r[local_M][j] = -(2. * delta_h_w_ij_x_left_part(w, local_M, j, global_x, global_y, h_1, h_2)
                + delta_h_w_ij_y_top_part(w, local_M, j, global_x, global_y, h_1, h_2)
                + delta_h_w_ij_y_bottom_part(w, local_M, j, global_x, global_y, h_1, h_2))
                + w[local_M][j] * ((2. / h_1));
            
        }   
        else {
            r[local_M][j] = -(delta_h_w_ij_x_left_part(w, local_M, j, global_x, global_y, h_1, h_2)
                + delta_h_w_ij_y_top_part(w, local_M, j, global_x, global_y, h_1, h_2)
                + delta_h_w_ij_y_bottom_part(w, local_M, j, global_x, global_y, h_1, h_2)
                + (1. / (h_1 * h_1) * k_3(global_x + h_1 * 0.5, global_y) *
                    (left_column_of_right_neighbor[j] - w[local_M][j])));
        }
        
    }
    // ОБРАБОТКА ЛЕВОГО ВЕРХНЕГО УГЛА ЛОКАЛЬНОГО ПРЯМОУГОЛЬНИКА
    global_x = local_x_start * h_1;
    global_y = local_y_start * h_2;
    // Если левый нижний угол находится на левой/нижней границе -> условие 1 (phi)
    if (bottom_neighbor_proc_id == -1 || left_neighbor_proc_id == -1) {

        r[0][0] = phi(local_x_start * h_1, local_y_start * h_2);
    }

    // Левый нижний угол локального не касается границ глобального прямоугольника
    if (bottom_neighbor_proc_id != -1 && left_neighbor_proc_id != -1) {
        // Левый нижний угол

        r[0][0] = -(delta_h_w_ij_x_right_part(w, 0, 0, global_x, global_y, h_1, h_2)
            + delta_h_w_ij_y_top_part(w, 0, 0, global_x, global_y, h_1, h_2)
            - (1. / (h_1 * h_1) * k_3(global_x - h_1 * 0.5, global_y)) * (w[0][0] - right_column_of_left_neighbor[0])
            - (1. / (h_2 * h_2) * k_3(global_x, global_y - h_2 * 0.5)) * (w[0][0] - top_row_of_bottom_neighbor[0]));
    
    }
    
    // ОБРАБОТКА ЛЕВОГО ВЕРХНЕГО УГЛА ЛОКАЛЬНОГО ПРЯМОУГОЛЬНИКА
    global_x = local_x_start * h_1;
    global_y = (local_y_start + local_N) * h_2;
    // Если процесс содержит левый верхний угол глобального прямоугольника
    if (top_neighbor_proc_id == -1 && left_neighbor_proc_id == -1) {

        r[0][local_N] = phi(local_x_start * h_1, (local_y_start + local_N) * h_2);
    }
    
    // Глобальная верхняя граница (тип 3), но сосед слева есть
    if (top_neighbor_proc_id == -1 && left_neighbor_proc_id != -1) {
        r[0][local_N] = -(delta_h_w_ij_x_right_part(w, 0, local_N, global_x, global_y, h_1, h_2)
            + delta_h_w_ij_y_bottom_part(w, 0, local_N, global_x, global_y, h_1, h_2)
            - (1. / (h_1 * h_1) * k_3(global_x - h_1 * 0.5, global_y)) * (w[0][local_N] - right_column_of_left_neighbor[local_N]))
            + w[0][local_N] * (2. / h_2);

    }
    // Глобальная левая граница (тип 1), но сосед сверху есть
    if (top_neighbor_proc_id != -1 && left_neighbor_proc_id == -1) {
        r[0][local_N] = phi(local_x_start * h_1, (local_y_start + local_N) * h_2);
    }
    // Левый верхний угол локального не касается границ глобального прямоугольника
    if (top_neighbor_proc_id != -1 && left_neighbor_proc_id != -1) {
        
        r[0][local_N] = -(delta_h_w_ij_x_right_part(w, 0, local_N, global_x, global_y, h_1, h_2)
            + delta_h_w_ij_y_bottom_part(w, 0, local_N, global_x, global_y, h_1, h_2)
            - (1. / (h_1 * h_1) * k_3(global_x - h_1 * 0.5, global_y)) * (w[0][local_N] - right_column_of_left_neighbor[local_N])
            + (1. / (h_2 * h_2) * k_3(global_x, global_y + h_2 * 0.5)) * (bottom_row_of_top_neighbor[0] - w[0][local_N]));
 
    }
    
    // ОБРАБОТКА ПРАВОГО НИЖНЕГО УГЛА ЛОКАЛЬНОГО ПРЯМОУГОЛЬНИКА
    global_x = (local_x_start + local_M) * h_1;
    global_y = local_y_start * h_2;
    // Если процесс содержит правый нижний угол глобального прямоугольника
    if (bottom_neighbor_proc_id == -1 && right_neighbor_proc_id == -1) {

        r[local_M][0] = phi((local_x_start + local_M) * h_1, local_y_start * h_2);
    }
    
    // Глобальная нижняя граница (тип 1), но сосед справа есть
    if (bottom_neighbor_proc_id == -1 && right_neighbor_proc_id != -1) {
        r[local_M][0] = phi((local_x_start + local_M) * h_1, local_y_start * h_2);

    }

    // Глобальная правая граница (тип 3), но сосед снизу есть
    if (bottom_neighbor_proc_id != -1 && right_neighbor_proc_id == -1) {
        r[local_M][0] = -(delta_h_w_ij_x_left_part(w, local_M, 0, global_x, global_y, h_1, h_2)
            + delta_h_w_ij_y_top_part(w, local_M, 0, global_x, global_y, h_1, h_2)
            - (1. / (h_2 * h_2) * k_3(global_x, global_y + h_2 * 0.5)) * (w[local_M][0] - top_row_of_bottom_neighbor[local_M]))
            + w[local_M][0] * (2. / h_1);

    }
    // Правый нижний угол локального не касается границ глобального прямоугольника
    if (bottom_neighbor_proc_id != -1 && right_neighbor_proc_id != -1) {
        
        r[local_M][0] = -(delta_h_w_ij_x_left_part(w, local_M, 0, global_x, global_y, h_1, h_2)
            + delta_h_w_ij_y_top_part(w, local_M, 0, global_x, global_y, h_1, h_2)
            + (1. / (h_1 * h_1) * k_3(global_x + h_1 * 0.5, global_y)) * (left_column_of_right_neighbor[0] - w[local_M][0])
            - (1. / (h_2 * h_2) * k_3(global_x, global_y + h_2 * 0.5)) * (w[local_M][0] - top_row_of_bottom_neighbor[local_M] ));

    }

    // ОБРАБОТКА ПРАВОГО ВЕРХНЕГО УГЛА ЛОКАЛЬНОГО ПРЯМОУГОЛЬНИКА
    // Если процесс содержит правый верхний угол глобального прямоугольника
    global_x = (local_x_start + local_M) * h_1;
    global_y = (local_y_start + local_N) * h_2;
    if (top_neighbor_proc_id == -1 && right_neighbor_proc_id == -1) {

        r[local_M][local_N] = -(2. * delta_h_w_ij_x_left_part(w, local_M, local_N, global_x, global_y, h_1, h_2)
            + 2. * delta_h_w_ij_y_bottom_part(w, local_M, local_N, global_x, global_y, h_1, h_2))
            + ((2. / h_1) + (2. / h_2)) * w[local_M][local_N];

    }
    
    // Глобальная верхняя граница (тип 3), но сосед справа есть
    if (top_neighbor_proc_id == -1 && right_neighbor_proc_id != -1) {
        r[local_M][local_N] = -(delta_h_w_ij_x_left_part(w, local_M, local_N, global_x, global_y, h_1, h_2)
            + delta_h_w_ij_y_bottom_part(w, local_M, local_N, global_x, global_y, h_1, h_2)
            + (1. / (h_1 * h_1) * k_3(global_x + h_1 * 0.5, global_y)) * (left_column_of_right_neighbor[local_N] - w[local_M][local_N]))
            + w[local_M][local_N] * (2. / h_2);
    }
    // Глобальная правая граница (тип 3), но сосед сверху есть
    if (top_neighbor_proc_id != -1 && right_neighbor_proc_id == -1) {
        r[local_M][local_N] = -(delta_h_w_ij_x_left_part(w, local_M, local_N, global_x, global_y, h_1, h_2)
            + delta_h_w_ij_y_bottom_part(w, local_M, local_N, global_x, global_y, h_1, h_2)
            + (1. / (h_2 * h_2) * k_3(global_x, global_y + h_2 * 0.5)) * (bottom_row_of_top_neighbor[local_M] - w[local_M][local_N]))
            + w[local_M][local_N] * (2. / h_1);
    }

    // Правый верхний угол локального не касается границ глобального прямоугольника
    if (top_neighbor_proc_id != -1 && right_neighbor_proc_id != -1) {
        r[local_M][local_N] = -(delta_h_w_ij_x_left_part(w, local_M, local_N, global_x, global_y, h_1, h_2)
            + delta_h_w_ij_y_bottom_part(w, local_M, local_N, global_x, global_y, h_1, h_2)
            + (1. / (h_1 * h_1) * k_3(global_x + h_1 * 0.5, global_y)) * (left_column_of_right_neighbor[local_N] - w[local_M][local_N])
            + (1. / (h_2 * h_2) * k_3(global_x, global_y + h_2 * 0.5)) * (bottom_row_of_top_neighbor[local_M] - w[local_M][local_N]));

    }
}


void equation_right_part(double** B, int local_M, int local_N, int local_x_start, int local_y_start,
    int left_neighbor_proc_id, int right_neighbor_proc_id, int top_neighbor_proc_id, int bottom_neighbor_proc_id,
     double h_1, double h_2) {

    for (int i = 0; i <= local_M; i++) {
        for (int j = 0; j <= local_N; j++) {
            B[i][j] = F((local_x_start + i) * h_1, (local_y_start + j) * h_2);
        }
    }
    // Если процесс содержит нижнюю границу глобального прямоугольника
    if (bottom_neighbor_proc_id == -1) {
        for (int i = 0; i <= local_M; i++) {
            B[i][0] = phi((local_x_start + i) * h_1, local_y_start * h_2);
        }
    }
    // Если процесс содержит верхнюю границу глобального прямоугольника
    if (top_neighbor_proc_id == -1) {
        for (int i = 0; i <= local_M; i++) {
            B[i][local_N] = F((local_x_start + i) * h_1, (local_y_start + local_N) * h_2) 
                + (2. / h_2) * psi_top_border((local_x_start + i) * h_1, (local_y_start + local_N) * h_2);
        }

    }
    // Если процесс содержит левую границу глобального прямоугольника
    if (left_neighbor_proc_id == -1) {
        for (int j = 0; j <= local_N; j++) {
            B[0][j] = phi(local_x_start * h_1, (local_y_start + j) * h_2);
        }

    }
    // Если процесс содержит правую границу глобального прямоугольника
    if (right_neighbor_proc_id == -1) {
        for (int j = 0; j <= local_N; j++) {
            B[local_M][j] = F((local_x_start + local_M) * h_1, (local_y_start + j) * h_2) + (2. / h_1) * psi_right_border((local_x_start + local_M) * h_1, (local_y_start + j) * h_2);
        }

    }
    // Если процесс содержит левый нижний угол глобального прямоугольника
    if (bottom_neighbor_proc_id == -1 && left_neighbor_proc_id == -1) {
        B[0][0] = phi(local_x_start * h_1, local_y_start * h_2);
    }
    // Если процесс содержит левый верхний угол глобального прямоугольника
    if (top_neighbor_proc_id == -1 && left_neighbor_proc_id == -1) {
        B[0][local_N] = phi(local_x_start * h_1, (local_y_start + local_N) * h_2);
    }
    // Если процесс содержит правый нижний угол глобального прямоугольника
    if (bottom_neighbor_proc_id == -1 && right_neighbor_proc_id == -1) {
        B[local_M][0] = phi((local_x_start + local_M) * h_1, local_y_start * h_2);
    }
    // Если процесс содержит правый верхний угол глобального прямоугольника
    if (top_neighbor_proc_id == -1 && right_neighbor_proc_id == -1) {
        // F(M * h_1, N * h_2) + (2. / h_1) * psi_right_border(M * h_1, N * h_2) + (2. / h_2) * psi_top_border(M * h_1, N * h_2);
        B[local_M][local_N] = F((local_x_start + local_M) * h_1, (local_y_start + local_N) * h_2) 
            + (2. / h_1) * psi_right_border((local_x_start + local_M) * h_1, (local_y_start + local_N) * h_2) 
            + (2. / h_2)  * psi_top_border((local_x_start + local_M) * h_1, (local_y_start + local_N) * h_2);
            // phi((local_x_start + local_M) * h_1, (local_y_start + local_N) * h_2);
    }

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

// Вычисление локального куска скалярного произведения
double local_dot_product(double** u, double** v, double h_1, double h_2, 
    int local_M, int local_N, int local_x_start, int local_y_start,
    int M, int N) {

    double dot_p = 0.0;
    double inner_sum;
    int global_x, global_y;
    for (int i = 0; i <= local_M; i++) {
        global_x = local_x_start + i;
        inner_sum = 0.;
        for (int j = 0; j <= local_N; j++) {
            global_y = local_y_start + j;
            inner_sum += h_2 * rho(global_x, global_y, M, N) * u[i][j] * v[i][j];
        }
        dot_p += h_1 * inner_sum;
    }

    return dot_p;
}


void borders_exchange(double** local_area,
    double* my_top_row, double* my_bottom_row, double* my_right_column, double* my_left_column,
    double* top_row_of_bottom_neighbor,  double* bottom_row_of_top_neighbor, double* left_column_of_right_neighbor, double* right_column_of_left_neighbor,
    int left_neighbor_proc_id, int right_neighbor_proc_id, int top_neighbor_proc_id, int bottom_neighbor_proc_id,
    int local_M, int local_N, MPI_Comm commutator) {

    for (int i = 0; i <= local_M; i++) {
        my_bottom_row[i] = local_area[i][0];
        my_top_row[i] = local_area[i][local_N];
    }
    for (int j = 0; j <= local_N; j++) {
        my_left_column[j] = local_area[0][j];
        my_right_column[j] = local_area[local_M][j];
    }
    MPI_Status status;

    if (top_neighbor_proc_id != -1) {
        MPI_Send(my_top_row, local_M + 1, MPI_DOUBLE, top_neighbor_proc_id, 0, commutator);
    }
    if (bottom_neighbor_proc_id != -1) {
        MPI_Send(my_bottom_row, local_M + 1, MPI_DOUBLE, bottom_neighbor_proc_id, 0, commutator);
    }
    if (left_neighbor_proc_id != -1) {
        MPI_Send(my_left_column, local_N + 1, MPI_DOUBLE, left_neighbor_proc_id, 0, commutator);
    }
    if (right_neighbor_proc_id != -1) {
        MPI_Send(my_right_column, local_N + 1, MPI_DOUBLE, right_neighbor_proc_id, 0, commutator);
    }
    if (top_neighbor_proc_id != -1) {
        MPI_Recv(bottom_row_of_top_neighbor, local_M + 1, MPI_DOUBLE, top_neighbor_proc_id, MPI_ANY_TAG, commutator, &status);
    }
    if (bottom_neighbor_proc_id != -1) {
        MPI_Recv(top_row_of_bottom_neighbor, local_M + 1, MPI_DOUBLE, bottom_neighbor_proc_id, MPI_ANY_TAG, commutator, &status);
    }
    if (left_neighbor_proc_id != -1) {
        MPI_Recv(right_column_of_left_neighbor, local_N + 1, MPI_DOUBLE, left_neighbor_proc_id, MPI_ANY_TAG, commutator, &status);
    }
    if (right_neighbor_proc_id != -1) {
        MPI_Recv(left_column_of_right_neighbor, local_N + 1, MPI_DOUBLE, right_neighbor_proc_id, MPI_ANY_TAG, commutator, &status);
    }

}




int main(int argc, char** argv) {

    const int M = atoi(argv[1]);
    const int N = atoi(argv[2]);

    const int M = 120;
    const int N = 120;


    double EPS = 1e-6;
    const double h_1 = 2. / (double)M;
    const double h_2 = 1. / (double)N;

    int left_neighbor_proc_id;
    int right_neighbor_proc_id;
    int bottom_neighbor_proc_id;
    int top_neighbor_proc_id;

    int local_M;
    int local_N;
    int rx;
    int ry;
    int coords[2];

    int local_x_start;
    int local_y_start;
    int local_x_end;
    int local_y_end;

    MPI_Comm commutator;
    MPI_Init(&argc, &argv);


    int rank, size;
    int power, px, py;
    int dims[2];

    const int ndims = 2;
    int periods[2] = { 0, 0 };

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if ((power = log2_(size)) == -1) {
        if (rank == 0)
            printf("error");
        MPI_Finalize();
    }

    px = split(M, N, power);
    py = power - px;

    dims[0] = pow(2, px);
    dims[1] = pow(2, py);

    local_M = M / dims[0];
    local_N = N / dims[1];
    rx = M + 1 - dims[0] * local_M;
    ry = N + 1 - dims[1] * local_N;

    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, 1, &commutator);

    MPI_Cart_coords(commutator, rank, ndims, coords);
    // Находим глобальные координаты границ локального прямоугольника
    local_x_start = MIN(rx, coords[0]) * (local_M + 1) + MAX(0, (coords[0] - rx)) * local_M;
    local_x_end = local_x_start + local_M + (coords[0] < rx ? 1 : 0);
    local_y_start = MIN(ry, coords[1]) * (local_N + 1) + MAX(0, (coords[1] - ry)) * local_N;
    local_y_end = local_y_start + local_N + (coords[1] < ry ? 1 : 0);

    local_M = local_x_end - local_x_start;
    local_N = local_y_end - local_y_start;

    

    
    rank = rank;
    coords[0] = coords[0];
    coords[1] = coords[1];

    local_x_start = local_x_start;
    local_y_start = local_y_start;
    size = dims[0] * dims[1];

    // Находим соседей вершины
    MPI_Cart_shift(commutator, 0, 1, &left_neighbor_proc_id, &right_neighbor_proc_id);
    MPI_Cart_shift(commutator, 1, -1, &top_neighbor_proc_id, &bottom_neighbor_proc_id);

    double start_time = MPI_Wtime();

    double tau = 0.;

    double** local_B = (double**)malloc((local_M + 1) * sizeof(double*));
    for (int i = 0; i <= local_M; i++)
        local_B[i] = (double*)malloc((local_N + 1) * sizeof(double));
    double** local_w = (double**)malloc((local_M + 1) * sizeof(double*));
    for (int i = 0; i <= local_M; i++)
        local_w[i] = (double*)malloc((local_N + 1) * sizeof(double));
    double** local_r = (double**)malloc((local_M + 1) * sizeof(double*));
    for (int i = 0; i <= local_M; i++)
        local_r[i] = (double*)malloc((local_N + 1) * sizeof(double));
    double** local_current_w = (double**)malloc((local_M + 1) * sizeof(double*));
    for (int i = 0; i <= local_M; i++)
        local_current_w[i] = (double*)malloc((local_N + 1) * sizeof(double));
    double** local_real_u = (double**)malloc((local_M + 1) * sizeof(double*));
    for (int i = 0; i <= local_M; i++)
        local_real_u[i] = (double*)malloc((local_N + 1) * sizeof(double));
    double** local_Ar = (double**)malloc((local_M + 1) * sizeof(double*));
    for (int i = 0; i <= local_M; i++)
        local_Ar[i] = (double*)malloc((local_N + 1) * sizeof(double));
    double** real_u = (double**)malloc((local_M + 1) * sizeof(double*));
    for (int i = 0; i <= local_M; i++)
        real_u[i] = (double*)malloc((local_N + 1) * sizeof(double));
    double** local_error_matrix = (double**)malloc((local_M + 1) * sizeof(double*));
    for (int i = 0; i <= local_M; i++)
        local_error_matrix[i] = (double*)malloc((local_N + 1) * sizeof(double));

    double* my_top_row = (double*)malloc((local_M + 1) * sizeof(double));
    double* my_bottom_row = (double*)malloc((local_M + 1) * sizeof(double));
    double* my_right_column = (double*)malloc((local_N + 1) * sizeof(double));
    double* my_left_column = (double*)malloc((local_N + 1) * sizeof(double));

    double* top_row_of_bottom_neighbor = (double*)malloc((local_M + 1) * sizeof(double));
    double* bottom_row_of_top_neighbor = (double*)malloc((local_M + 1) * sizeof(double));
    double* left_column_of_right_neighbor = (double*)malloc((local_N + 1) * sizeof(double));
    double* right_column_of_left_neighbor = (double*)malloc((local_N + 1) * sizeof(double));

    double local_Ar_r_dot_product_part = 0.;
    double global_Ar_r_dot_product = 0.;

    double local_Ar_Ar_dot_product_part = 0.;
    double global_Ar_Ar_dot_product = 0.;

    double local_last_solution_update_norm = 0.;
    double global_last_solution_update_norm = 0.;
    double execution_time;

    double difference_local = EPS;
    double difference_global = EPS;

    init_grid(local_w, h_1, h_2, local_M, local_N, local_x_start, local_y_start);



    equation_right_part(local_B, local_M, local_N, local_x_start, local_y_start,
        left_neighbor_proc_id, right_neighbor_proc_id, top_neighbor_proc_id, bottom_neighbor_proc_id,
        h_1, h_2);

    bool keep_calculation = true;
    while (keep_calculation)
    {
        borders_exchange(local_w,
            my_top_row, my_bottom_row, my_right_column, my_left_column,
            top_row_of_bottom_neighbor, bottom_row_of_top_neighbor, left_column_of_right_neighbor, right_column_of_left_neighbor,
            left_neighbor_proc_id, right_neighbor_proc_id, top_neighbor_proc_id, bottom_neighbor_proc_id,
            local_M, local_N, commutator);

        equation_left_part(local_w, top_row_of_bottom_neighbor, bottom_row_of_top_neighbor,
            left_column_of_right_neighbor, right_column_of_left_neighbor,
            local_M, local_N, local_x_start, local_y_start,
            left_neighbor_proc_id, right_neighbor_proc_id, top_neighbor_proc_id, bottom_neighbor_proc_id,
            local_r, h_1, h_2);



        for (int i = 0; i <= local_M; i++)
            for (int j = 0; j <= local_N; j++) {
                local_r[i][j] -= local_B[i][j];
                local_current_w[i][j] = local_w[i][j];
            }

        borders_exchange(local_r,
            my_top_row, my_bottom_row, my_right_column, my_left_column,
            top_row_of_bottom_neighbor, bottom_row_of_top_neighbor, left_column_of_right_neighbor, right_column_of_left_neighbor,
            left_neighbor_proc_id, right_neighbor_proc_id, top_neighbor_proc_id, bottom_neighbor_proc_id,
            local_M, local_N, commutator);
        equation_left_part(local_r, top_row_of_bottom_neighbor, bottom_row_of_top_neighbor,
            left_column_of_right_neighbor, right_column_of_left_neighbor,
            local_M, local_N, local_x_start, local_y_start,
            left_neighbor_proc_id, right_neighbor_proc_id, top_neighbor_proc_id, bottom_neighbor_proc_id,
            local_Ar, h_1, h_2);

        local_Ar_r_dot_product_part = local_dot_product(local_Ar, local_r, h_1, h_2,
            local_M, local_N, local_x_start, local_y_start, M, N);
        MPI_Reduce(&local_Ar_r_dot_product_part, &global_Ar_r_dot_product, 1, MPI_DOUBLE, MPI_SUM, 0, commutator);
        MPI_Bcast(&global_Ar_r_dot_product, 1, MPI_DOUBLE, 0, commutator);

        local_Ar_Ar_dot_product_part = local_dot_product(local_Ar, local_Ar, h_1, h_2,
            local_M, local_N, local_x_start, local_y_start, M, N);
        MPI_Reduce(&local_Ar_Ar_dot_product_part, &global_Ar_Ar_dot_product, 1, MPI_DOUBLE, MPI_SUM, 0, commutator);
        MPI_Bcast(&global_Ar_Ar_dot_product, 1, MPI_DOUBLE, 0, commutator);

        tau = global_Ar_r_dot_product / global_Ar_Ar_dot_product;

        // Вычисление следующего значения w
        for (int i = 0; i <= local_M; i++)
            for (int j = 0; j <= local_N; j++)
                local_w[i][j] = local_w[i][j] - tau * local_r[i][j];
        for (int i = 0; i <= local_M; i++)
            for (int j = 0; j <= local_N; j++) {
                local_current_w[i][j] = local_w[i][j] - local_current_w[i][j];
            }

        local_last_solution_update_norm = sqrt(local_dot_product(local_current_w, local_current_w, h_1, h_2,
            local_M, local_N, local_x_start, local_y_start, M, N));
   
        MPI_Reduce(&local_last_solution_update_norm, &global_last_solution_update_norm, 1, MPI_DOUBLE, MPI_SUM, 0, commutator);
        MPI_Bcast(&global_last_solution_update_norm, 1, MPI_DOUBLE, 0, commutator);


        if (global_last_solution_update_norm < EPS) {
            keep_calculation = false;
        }
    }

    execution_time = MPI_Wtime() - start_time;


    for (int i = 0; i <= local_M; i++) {
        for (int j = 0; j <= local_N; j++) {
            real_u[i][j] = u_4((local_x_start + i) * h_1, (local_y_start + j) * h_2);
        }
    }
    /*
    if (rank == 0) {
        printf("REAL U \n");
        printf("RANK %d, local_x_start %d, local_y_start %d, local_M %d, local_N %d\n", rank, local_x_start, local_y_start, local_M, local_N);
        for (int i = 0; i <= local_M; i++) {
            for (int j = 0; j <= local_N; j++) {
                printf("%.8f ", real_u[i][j]);
            }
            printf("\n");
        }

        printf("LOCAL W \n");
        for (int i = 0; i <= local_M; i++) {
            for (int j = 0; j <= local_N; j++) {
                printf("%.8f ", local_w[i][j]);
            }
            printf("\n");
        }
    }
    */

    for (int i = 0; i <= local_M; i++)
        for (int j = 0; j <= local_N; j++)
            local_error_matrix[i][j] = real_u[i][j] - local_w[i][j];
    double local_error_scalar = local_dot_product(local_error_matrix, local_error_matrix, h_1, h_2,
        local_M, local_N, local_x_start, local_y_start, M, N);

    double global_error;
    MPI_Reduce(&local_error_scalar, &global_error, 1, MPI_DOUBLE, MPI_SUM, 0, commutator);
    global_error = sqrt(global_error);

    if (rank == 0) {
        printf("ERROR %.10f\n", global_error);
        printf("Execution time: %.10f\n", execution_time);
        // printf("size ,  M , N   , time        , boost      , max_diff\n");
        // printf("%d   &  %d \\times %d & %.10f  & %.10f\n", size, M, N, global_time_diff, global_error);
    }

    MPI_Finalize();
}


