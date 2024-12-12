Для запуска кода OMP:

Задаю в файле mf.cpp нужные значения delta, количества нитей, название файла, в который будет выводиться результат (в случае необходимости вывода итоговой матрицы W, раскомменчиваю соответствующий кусок кода)
g++ -fopenmp mf.cpp -o task
bsub < OpenMP_job.lsf

Для запуска кода MPI:
Компилировать файл task2/mf_mpi.cpp, задать в *.lsf названий выходного файла
mpicxx -o task2 task2/mf_mpi.cpp 
bsub < task2/MPI_JOB.lsf

Для запуска кода MPI+OpenMP:
Аналогично предыдущему пункту для файла task3/mf_mpi.cpp
mpicxx -o task2 task3/mf_mpi.cpp -fopenmp
bsub < OMP_MPI.lsf
