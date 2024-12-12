#include <iostream> 
#include <omp.h>
#include <vector>
#include <mpi.h>
#include <algorithm>
#include <fstream>
#include <cmath>

typedef std::vector<double> base_key;
typedef std::vector<std::vector<double> > key;


void fill_linspace(base_key &x, double a, double b, int n, int block_size, int idx_start, int32_t rank) {
    double h = (b - a) / n;
    for (int i = 0; i <= n; ++i) {
        if ((i >= block_size * idx_start) && ((i < block_size * (idx_start + 1)) || i == n)) {
            x.push_back(a + i * h);
        }
    }
    // std::cout << "Rank " <<  rank << " block_size " << block_size << " idx_start " << idx_start << " size " << x.size() << std::endl;
    // std::cout << "rank " << rank << " start_idx " << block_size * idx_start << " start_val " << x[0] << " end_idx " << block_size * (idx_start + 1) << " end_val " << x[x.size() - 1] << std::endl;
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
    
    MPI_Allreduce(MPI_IN_PLACE, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return res;
}

double norm(const key &a, int m, int n, double h1, double h2) {
    return std::sqrt(sc_product(a, a, m, n, h1, h2));
}

void LinOp(key &res, const key &A, const key &B, const key &w, int m, int n, double h1, double h2, int32_t x_loc, int32_t y_loc, int32_t rank, int32_t dims[2],
            base_key &up_get, base_key &up_put, base_key &down_get, base_key &down_put, base_key &left_get, base_key &left_put, base_key &right_get, base_key &right_put) {
    MPI_Status status;
    // std::cout << "huy" << std::endl;
    // swap data
    if (x_loc % 2) {
        // put
        if (x_loc != 0) { // not left
            left_put = w[0];
            
            MPI_Send(left_put.data(), left_put.size(), MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            // std::cout << " Send size " << left_put.size() << " rank " << rank << "hhhhhHHHHHH" << std::endl;
            MPI_Recv(left_get.data(), left_get.size(), MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
            // std::cout << "KKKKKKKKKKK" << std::endl;
            // std::cout << "Murrr" << std::endl;
        }
        
        // std::cout << "hhhhh" << std::endl;


        if (x_loc != dims[0] - 1) { // not right
            MPI_Recv(right_get.data(), right_get.size(), MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);

            right_put = w[w.size()-1];
            MPI_Send(right_put.data(), right_put.size(), MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }

    } 
    else {
        if (x_loc != dims[0] - 1) { // not right
            // std::cout << " Get size " << right_get.size() << " rank " << rank << " Murrr " << std::endl;
            MPI_Recv(right_get.data(), right_get.size(), MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
            // std::cout << "Murrr" << std::endl;
            right_put = w[w.size()-1];
            MPI_Send(right_put.data(), right_put.size(), MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }

        if (x_loc != 0) { // not left
            left_put = w[0];
            
            MPI_Send(left_put.data(), left_put.size(), MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(left_get.data(), left_get.size(), MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
        }
        

    }

    // std::cout << "huuuuuuy" << std::endl;

    if (y_loc % 2) {
        // put
        if (y_loc != 0) { // not up
            // up_put
            for (int k = 0; k < w.size(); ++k) {
                up_put[k] = w[k][0];
            }
            
            MPI_Send(up_put.data(), up_put.size(), MPI_DOUBLE, rank - dims[0], 0, MPI_COMM_WORLD);
            MPI_Recv(up_get.data(), up_get.size(), MPI_DOUBLE, rank - dims[0], 0, MPI_COMM_WORLD, &status);
        }
        
        if (y_loc != dims[1] - 1) { // not right
            MPI_Recv(down_get.data(), down_get.size(), MPI_DOUBLE, rank + dims[0], 0, MPI_COMM_WORLD, &status);

            for (int k = 0; k < w.size(); ++k) {
                down_put[k] = w[k][w[k].size() - 1];
            }
            MPI_Send(down_put.data(), down_put.size(), MPI_DOUBLE, rank + dims[0], 0, MPI_COMM_WORLD);
        }

    } 
    else {
        if (y_loc != dims[1] - 1) { // not right
            MPI_Recv(down_get.data(), down_get.size(), MPI_DOUBLE, rank + dims[0], 0, MPI_COMM_WORLD, &status);

            for (int k = 0; k < w.size(); ++k) {
                down_put[k] = w[k][w[k].size() - 1];
            }
            MPI_Send(down_put.data(), down_put.size(), MPI_DOUBLE, rank + dims[0], 0, MPI_COMM_WORLD);
        }

        if (y_loc != 0) { // not up
            // up_put
            for (int k = 0; k < w.size(); ++k) {
                up_put[k] = w[k][0];
            }
            
            MPI_Send(up_put.data(), up_put.size(), MPI_DOUBLE, rank - dims[0], 0, MPI_COMM_WORLD);
            MPI_Recv(up_get.data(), up_get.size(), MPI_DOUBLE, rank - dims[0], 0, MPI_COMM_WORLD, &status);
        }
        
        

    }

    // std::cout << "huy" << std::endl;
    #pragma omp parallel for
    for (int i = 0; i < m - 1; ++i) {
        for (int j = 0; j < n - 1; ++j) {
            double dx1, dx2, dy1, dy2;
            if(i == 0) {
                dx1 = (w[i][j] - left_get[j])/ h1;
            }
            else
                dx1 = (w[i][j] - w[i-1][j]) / h1;
            
            if(i == m - 1)
                dx2 = (right_get[j]-w[i][j]) / h1;
            else
                dx2 = (w[i+1][j] - w[i][j]) / h1;
            
            if(j == 0)
                dy1 = (w[i][j] - up_get[i]) / h2;
            else
                dy1 = (w[i][j] - w[i][j-1]) / h2;
            
            if(j == n - 1)
                dy2 = (down_get[i]-w[i][j]) / h2;
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



std::pair<int, double> IterAndTime(int M, int N, int num_threads, int max_iter, double delta, std::ofstream &myfile, int32_t size, int32_t rank, int32_t dims[2]) {
    // std::cout << "csnlnsdv" << std::endl;
    double a = 0, b = 3, c = 0, d = 3, tau; 
    double time_solve;

    double h1 = (b - a) / M, h2 = (d - c) / N;
    double eps = std::max(h1, h2) * std::max(h1, h2);

    base_key x, y, up_get, up_put, down_get, down_put, left_get, left_put, right_get, right_put;

    int32_t x_loc = rank % dims[0], y_loc = rank / dims[0];


    int DX = (M + dims[0]) / dims[0];
    fill_linspace(x, a, b, M, DX, x_loc, rank);
    M = x.size() - 1;
    up_get.resize(M, 0);
    up_put.resize(M, 0);
    down_get.resize(M, 0);
    down_put.resize(M, 0);
    // std::cout << dims[0] << " " << dims[1] << " rank " << rank << " x_loc " << x_loc << " y_loc " << y_loc << " M " << M << std::endl;

    int DY = (N + dims[1]) / dims[1];
    fill_linspace(y, c, d, N, DY, y_loc, rank);
    N = y.size() - 1;
    left_get.resize(N, 0);
    left_put.resize(N, 0);
    right_get.resize(N, 0);
    right_put.resize(N, 0);

    key A, B, F, W, r, Ar;

    fill_W(W, M, N);
    fill_W(r, M, N);
    fill_W(Ar, M, N);
    fill_A(A, x, y, N, M, h1, h2, eps);
    fill_B(B, x, y, N, M, h1, h2, eps);
    fill_F(F, x, y, N, M, h1, h2, eps);


    MPI_Barrier(MPI_COMM_WORLD);
    time_solve = MPI_Wtime();

    int iter = 0;
    omp_set_num_threads(num_threads);
    // double start_time = omp_get_wtime();
    // std::cout << iter << ' ' << tau*norm(r, M-1, N-1, h1, h2) << std::endl;
    // std::cout << dims[0] << ' ' << dims[1] << ' ' << W.size() << ' ' << W[0].size() << std::endl;
    while ((iter == 0 || (tau*norm(r, M-1, N-1, h1, h2) > delta)) && iter <= max_iter) {
        iter++;
        // std::cout << "cdns" << std::endl;
        LinOp(r, A, B, W, M, N, h1, h2, x_loc, y_loc, rank, dims, up_get, up_put, down_get, down_put, left_get, left_put, right_get, right_put);
        W_min_TauR(r, r, F, 1, M-1, N-1);
        LinOp(Ar, A, B, r, M, N, h1, h2, x_loc, y_loc, rank, dims, up_get, up_put, down_get, down_put, left_get, left_put, right_get, right_put);
        tau = sc_product(r, r, M-1, N-1, h1, h2) / sc_product(Ar, r, M-1, N-1, h1, h2);
        W_min_TauR(W, W, r, tau, M, N);
        // std::cout << iter << ' ' << tau*norm(r, M-1, N-1, h1, h2) << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    time_solve = MPI_Wtime() - time_solve;
    return std::make_pair(iter, time_solve);
}


int main(int argc, char *argv[]) {

    int max_iter = 10000000;
    double delta = 1e-7;
    
    int m1 = 10, n1 = 10, m2 = 20, n2 = 20;
    int n_threads_arr1[] = {1,2,4,8,16};
    int n_threads_arr2[] = {1,4,8,16,32};

    int num_threads = 1;

    std::ofstream myfile;
    // myfile.open ("table_of_iter_and_time_1_10_20.txt");
    // myfile << "  t | M |  N | iter | time\n";

    // int M = m1, N = n1;

    if (argc != 3){
        std::cerr << "Grid sizes are not set!" << std::endl;
        return 1;
    }

    int32_t M = atoi(argv[1]), N = atoi(argv[2]);


    MPI_Init(&argc, &argv);

    int32_t size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int32_t dims[2] = {0, 0};

    if (MPI_Dims_create(size, 2, dims) != MPI_SUCCESS){
        if (rank == 0){
            std::cerr << "Process grid cannot be created!" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    std::pair<int, double> a = IterAndTime(M, N, num_threads, max_iter, delta, myfile, size, rank, dims);

    double sum_time;
    MPI_Allreduce(&a.second,  &sum_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << " iters: " << a.first << " time: " << a.second << " sum time: " << sum_time / size << std::endl;
    }

    MPI_Finalize();
    
    // myfile << num_threads / 10 << num_threads % 10 << ' ' << M << ' ' << N << ' ' << a.first << ' ' << a.second << "\n";

    // myfile.close();

    return 0;
}