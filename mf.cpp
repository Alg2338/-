#include <iostream> 
#include <omp.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cmath>

typedef std::vector<double> base_key;
typedef std::vector<std::vector<double> > key;

void fill_linspace(base_key &x, double a, double b, int n) {
    double h = (b - a) / n;
    for (int i = 0; i <= n; ++i) {
        x.push_back(a + i * h);
    }
}

void fill_W(key &W, int M, int N) {
    
    for (int i = 0; i < M; ++i) {
        base_key tmp;
        W.push_back(tmp);
        for (int j = 0; j < N; ++j) {
            W[i].push_back(0);
        }
    }
    return;
}

void fill_A(key &A, const base_key &x, const base_key &y, int n, int m, double h1, double h2, double eps) {
    
    for (int i = 1; i <= m; ++i) {
        base_key tmp;
        A.push_back(tmp);
        for (int j = 1; j < n; ++j) {
            double cur_x = x[i] - h1 / 2;
            double min_y = y[j] - h2 / 2, max_y = y[j] + h2 / 2;

            if (max_y < -3 * (cur_x - 3)) {
                A[i-1].push_back(1);
            }
            else if (min_y > -3 * (cur_x - 3)) {
                A[i-1].push_back(1/eps);
            }
            else {
                double l = -3 * (cur_x - 3) - min_y;
                A[i-1].push_back(l / h2 + (1 - l / h2)/eps);
            }
            
        }
    }
}

void fill_B(key &B, const base_key &x, const base_key &y, int n, int m, double h1, double h2, double eps) {
    for (int i = 1; i < m; ++i) {
        base_key tmp;
        B.push_back(tmp);
        for (int j = 1; j <= n; ++j) {
            double cur_y = y[j] - h2 / 2;
            double min_x = x[i] - h1 / 2, max_x = x[i] + h1 / 2;

            if (min_x > 3-cur_y/3) {
                B[i-1].push_back(1/eps);
            }
            else if (max_x < 3 - cur_y / 3) {
                B[i-1].push_back(1);
            }
            else {
                double l = 3-cur_y/3 - min_x;
                B[i-1].push_back(l/h1 + (1-l/h1)/eps);
            }
        }
    }
}

void fill_F(key &F, const base_key &x, const base_key &y, int n, int m, double h1, double h2, double eps) {
    for (int i = 1; i < m; ++i) {
        base_key tmp;
        F.push_back(tmp);
        for (int j = 1; j < n; ++j) {
            double x1 = x[i] - h1/2, x2 = x[i] + h1/2;
            double y1 = y[j] - h2/2, y2 = y[j] + h2/2;

            if (y2 <= -3 * (x2-3)) 
                F[i-1].push_back(1);
            else if (y1 >= -3 * (x1 - 3))
                F[i-1].push_back(0);
            else if ((y2 > -3 * (x2 - 3)) && (y1 <= -3 * (x2 - 3)) && (y2 <= -3 * (x1 - 3))) {
                double al = y2 + 3 * (x2 - 3);
                double bet = x2 - (3 - y2 / 3);
                F[i-1].push_back(1 - al * bet / (2 * h1 * h2));
            }
            else if ((y2 > -3 * (x2 - 3)) && (y1 > -3 * (x2 - 3)) && (y2 <= -3 * (x1 - 3))) {
                double al = 3 - y2 / 3 - x1;
                double bet = 3 - y1 / 3 - x1;
                F[i-1].push_back((al + bet) / (2 * h1));
            }
            else if ((y2 > -3 * (x2 - 3)) && (y1 <= -3 * (x2 - 3)) && (y2 > -3 * (x1 - 3))) {
                double al = -3 * (x2 - 3) - y1;
                double bet = -3 * (x1 - 3) - y1;
                F[i-1].push_back((al + bet) / (2 * h2));
            }
            else {
                double al = 3 - y1/3 - x1;
                double bet = -3 * (x1 - 3) - y1;
                F[i-1].push_back(al * bet / (2 * h1 * h2));
            }
        }
    }
}

double sc_product(const key &a, const key &b, int m, int n, double h1, double h2) {
    double res = 0;
    #pragma omp parallel for reduction (+:res)
    for (int i = 0; i < m; ++i) 
        for (int j = 0; j < n; ++j) 
            res += a[i][j] * b[i][j] * h1 * h2;
    
    return res;
}

double norm(const key &a, int m, int n, double h1, double h2) {
    return std::sqrt(sc_product(a, a, m, n, h1, h2));
}

void LinOp(key &res, const key &A, const key &B, const key &w, int m, int n, double h1, double h2) {
    #pragma omp parallel for
    for (int i = 0; i < m - 1; ++i) {
        for (int j = 0; j < n - 1; ++j) {
            double dx1, dx2, dy1, dy2;
            if(i == 0)
                dx1 = w[i][j] / h1;
            else
                dx1 = (w[i][j] - w[i-1][j]) / h1;
            
            if(i == m - 1)
                dx2 = -w[i][j] / h1;
            else
                dx2 = (w[i+1][j] - w[i][j]) / h1;
            
            if(j == 0)
                dy1 = w[i][j] / h2;
            else
                dy1 = (w[i][j] - w[i][j-1]) / h2;
            
            if(j == n - 1)
                dy2 = -w[i][j] / h2;
            else
                dy2 = (w[i][j+1] - w[i][j]) / h2;

            res[i][j] = -(A[i+1][j]*dx2-A[i][j]*dx1) / h1 - (B[i][j+1]*dy2-B[i][j]*dy1) / h2;
        }
    }
}

void W_min_TauR(key &res, const key &w, const key &r, double tau, int m, int n) {
    #pragma omp parallel for
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            res[i][j] = w[i][j] - tau * r[i][j];
}



std::pair<int, double> IterAndTime(int M, int N, int num_threads, int max_iter, double delta, std::ofstream &myfile) {
    double a = 0, b = 3, c = 0, d = 3, tau; 

    double h1 = (b - a) / M, h2 = (d - c) / N;
    double eps = std::max(h1, h2) * std::max(h1, h2);

    base_key x, y;

    fill_linspace(x, a, b, M);

    fill_linspace(y, c, d, N);
    key A, B, F, W, r, Ar;

    fill_W(W, M, N);
    fill_W(r, M, N);
    fill_W(Ar, M, N);
    fill_A(A, x, y, N, M, h1, h2, eps);
    fill_B(B, x, y, N, M, h1, h2, eps);
    fill_F(F, x, y, N, M, h1, h2, eps);

    int iter = 0;
    omp_set_num_threads(num_threads);
    double start_time = omp_get_wtime();
    while ((iter == 0 || (tau*norm(r, M-1, N-1, h1, h2) > delta)) && iter <= max_iter) {
        iter++;
        LinOp(r, A, B, W, M, N, h1, h2);
        W_min_TauR(r, r, F, 1, M-1, N-1);
        LinOp(Ar, A, B, r, M, N, h1, h2);
        tau = sc_product(r, r, M-1, N-1, h1, h2) / sc_product(Ar, r, M-1, N-1, h1, h2);
        W_min_TauR(W, W, r, tau, M, N);
        // std::cout << iter << ' ' << tau*norm(r, M-1, N-1, h1, h2) << "\n";
    }

    // myfile << "\nW\n";

    // for (int i = 0; i < W.size(); i++) {
    //     for (int j = 0; j < W.size(); j++) {
    //         myfile << W[i][j] << " ";
    //     }
    //     myfile << "\n";
    // }

    double end_time = omp_get_wtime();
    return std::make_pair(iter, end_time - start_time);
}


int main(int argc, char *argv[]) {

    int max_iter = 10000000;
    double delta = 0.0000011;
    
    int m1 = 10, n1 = 10, m2 = 20, n2 = 20;
    int n_threads_arr1[] = {1,2,4,8,16};
    int n_threads_arr2[] = {1,4,8,16,32};

    int num_threads = 1;

    std::ofstream myfile;
    myfile.open ("table_of_iter_and_time_1_10_20.txt");
    myfile << "  t | M |  N | iter | time\n";

    int M = m1, N = n1;
    std::pair<int, double> a = IterAndTime(M, N, num_threads, max_iter, delta, myfile);
    myfile << num_threads / 10 << num_threads % 10 << ' ' << M << ' ' << N << ' ' << a.first << ' ' << a.second << "\n";

    // for (int i = 0; i < 5; ++i) {
    //     int M = m1, N = n1;
    //     int num_threads = n_threads_arr1[i];
    //     std::pair<int, double> a = IterAndTime(M, N, num_threads, max_iter, delta);
    //     myfile << num_threads / 10 << num_threads % 10 << ' ' << M << ' ' << N << ' ' << a.first << ' ' << a.second << "\n";
    // }

    M = m2;
    N = n2;
    delta = 0.0000001;
    std::pair<int, double> a2 = IterAndTime(M, N, num_threads, max_iter, delta, myfile);
    myfile << num_threads / 10 << num_threads % 10 << ' ' << M << ' ' << N << ' ' << a2.first << ' ' << a2.second << "\n";

    // for (int i = 0; i < 5; ++i) {
    //     int M = m2, N = n2;
    //     int num_threads = n_threads_arr2[i];
    //     std::pair<int, double> a = IterAndTime(M, N, num_threads, max_iter, delta);
    //     myfile << num_threads / 10 << num_threads % 10 << ' ' << M << ' ' << N << ' ' << a.first << ' ' << a.second << "\n";
    // }

    myfile.close();

    return 0;
}