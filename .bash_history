#1730061273
 bjobs -l 3045 
#1730061288
 bjobs -l 3045
#1730061293
 bjobs -l
#1730061295
 bjobs 
#1730061496
 bjobs -all
#1730058769
ls
#1730058772
pwd
#1730058824
cp -r /gpfs/quickstart ~
#1730058832
ls /gpfs/automountdir/
#1730059672
scp /Users/aleksey/Desktop/Суперкомпьютеры/tasks/task1/mf.cpp edu-cmc-skmodel24-615-02@polus.hpc.cs.msu.ru/mf1.cpp
#1730059680
scp /Users/aleksey/Desktop/Суперкомпьютеры/tasks/task1/mf.cpp edu-cmc-skmodel24-615-02@polus.hpc.cs.msu.ru
#1730059941
ls
#1730059942
pwd
#1730059966
scp /Users/aleksey/Desktop/Суперкомпьютеры/tasks/task1/mf.cpp edu-cmc-skmodel24-615-02@polus.hpc.cs.msu.ru:/home_edu/edu-cmc-skmodel24-615/edu-cmc-skmodel24-615-02
#1730059972
ls
#1730059986
pwd
#1730059987
ls
#1730059999
scp /Users/aleksey/Desktop/Суперкомпьютеры/tasks/task1/mf.cpp edu-cmc-skmodel24-615-02@polus.hpc.cs.msu.ru:/home_edu/edu-cmc-skmodel24-615/edu-cmc-skmodel24-615-02/
#1730060002
ls
#1730060013
scp /Users/aleksey/Desktop/Суперкомпьютеры/tasks/task1/mf.cpp edu-cmc-skmodel24-615-02@polus.hpc.cs.msu.ru:/home_edu/edu-cmc-skmodel24-615/edu-cmc-skmodel24-615-02/task1
#1730060068
scp /Users/aleksey/Desktop/Суперкомпьютеры/tasks/task1/mf.cpp 
#1730060074
ls
#1730060096
scp /Users/aleksey/Desktop/Суперкомпьютеры/tasks/task1/mf.cpp edu-cmc-skmodel24-615-02@polus.hpc.cs.msu.ru:/home_edu/edu-cmc-skmodel24-615/edu-cmc-skmodel24-615-02/
#1730060111
vim mf.cpp
#1730060127
ls
#1730060176
g++ mf.cpp -fopenmp -o task1
#1730060187
g++ mf.cpp -fopenmp 1
#1730060190
ls
#1730060208
g++ mf.cpp
#1730060384
vim mf.cpp
#1730060398
rm mf.cpp 
#1730060404
vim mf1.cpp
#1730060420
g++ mf.cpp -fopenmp
#1730060469
g++ -fopenmp mf1.cpp 
#1730060489
vim mf1.cpp
#1730060693
rm mf1.cpp 
#1730060696
vim mf1.cpp
#1730060705
g++ -fopenmp mf1.cpp 
#1730060720
g++ -fopenmp mf1.cpp -o task1
#1730060723
ls
#1730060972
vim OpenMP_job.lsf
#1730061052
ls
#1730061060
chmod +x OpenMP_job.lsf 
#1730061063
./OpenMP_job.lsf 
#1730061538
bsub < OpenMP_job.lsf
#1730061589
mv task1 my_job
#1730061598
vim OpenMP_job.lsf 
#1730061616
bsub < OpenMP_job.lsf
#1730061627
vim OpenMP_job.lsf 
#1730061693
bsub < OpenMP_job.lsf
#1730062122
vim OpenMP_job.lsf 
#1730062139
bsub < OpenMP_job.lsf
#1730062200
bjobs -l 1177951
#1730062398
ls
#1730062404
bjobs -l 1177951
#1730062419
ls
#1730062426
cat table_of_iter_and_time.txt 
#1730063037
bjobs -l 1177951
#1730063046
cat table_of_iter_and_time.txt 
#1730063050
ls
#1730063057
blobs
#1730063067
bjobs -all
#1730063229
ls
#1730063235
cat table_of_iter_and_time.txt 
#1730063243

#1730063259
cat OpenMP_job.lsf 
#1730063274
bjobs -all
#1730063328
vim OpenMP_job.lsf 
#1730144765
scp polus:~/mf1.cpp /Users/aleksey/Desktop/Суперкомпьютеры/tasks/task1
#1730144776
ls
#1730144786
rm mf1.cpp 
#1730144790
scp polus:~/mf1.cpp /Users/aleksey/Desktop/Суперкомпьютеры/tasks/task1
#1730144793
ls
#1730144838
scp /Users/aleksey/Desktop/Суперкомпьютеры/tasks/task1 polus:~/mf1.cpp
#1730144841
ls
#1730144852
scp /Users/aleksey/Desktop/Суперкомпьютеры/tasks/task1 polus:~/
#1730144867
scp polus:~/ /Users/aleksey/Desktop/Суперкомпьютеры/tasks/task1
#1730144870
ls
#1730144887
pwd
#1730144905
scp polus:~/home_edu/edu-cmc-skmodel24-615/edu-cmc-skmodel24-615-02 /Users/aleksey/Desktop/Суперкомпьютеры/tasks/task1
#1730144998
scp /Users/aleksey/Desktop/Суперкомпьютеры/tasks/task1 polus:~/home_edu/edu-cmc-skmodel24-615/edu-cmc-skmodel24-615-02 
#1730145262
scp /Users/aleksey/Desktop/Суперкомпьютеры/tasks/task1 edu-cmc-skmodel24-615-02@polus.hpc.cs.msu.ru:~/home_edu/edu-cmc-skmodel24-615/edu-cmc-skmodel24-615-02 
#1730145278
scp /Users/aleksey/Desktop/Суперкомпьютеры/tasks/task1 polus.hpc.cs.msu.ru:~/home_edu/edu-cmc-skmodel24-615/edu-cmc-skmodel24-615-02 
#1730145512
scp /Users/aleksey/Desktop/Суперкомпьютеры/tasks/task1/mf.cpp polus.hpc.cs.msu.ru:~/home_edu/edu-cmc-skmodel24-615/edu-cmc-skmodel24-615-02
#1730145834
q
#1730145838
quit
#1730145878
s
#1730145881
ls
#1730146786
vim OpenMP_job.lsf 
#1730147109
g++ mf.cpp -o my_job -fopenmp
#1730147160
bsub < OpenMP_job.lsf
#1730147170
bjobs
#1730147185
bjobs -l 1178685
#1730238751
ls
#1730238758
cat table_of_iter_and_time.txt 
#1730238769
ls
#1730238781
vim mf.cpp 
#1730238946
bjobs -l 1178685
#1730238956
bjobs
#1730239058
rm table_of_iter_and_time.txt 
#1730239059
ls
#1730239241
cat my_job.1177951.err 
#1730239277
rm my_job.1177951.err 
#1730239281
rm my_job.1177951.out 
#1730239287
ls
#1730239302
rm my_job*
#1730239305
ls
#1730239310
rm a.out 
#1730239317
cat OpenMP_job.lsf 
#1730319393
ды
#1730319395
ls
#1730321537
rm mf.cpp 
#1730321568
ls
#1730321585
g++ -fopenmp mf.cpp
#1730321890

#1730321901
bjobs
#1730321905
ls
#1730321928

#1730321974
g++ -fopenmp mf.cpp -o my
#1730321980
g++ -fopenmp mf.cpp -o my_job
#1730322004
bsub < OpenMP_job.lsf
#1730322010
bjobs
#1730322710
ls
#1730322736
g++ -fopenmp mf.cpp -o my_job
#1730322738
bjobs
#1730322741
bsub < OpenMP_job.lsf
#1730322744
bjobs
#1730323156
g++ -fopenmp mf.cpp -o my_job
#1730323159
bsub < OpenMP_job.lsf
#1730323164
bjobs
#1730323778
g++ -fopenmp mf.cpp -o my_job
#1730323780
bsub < OpenMP_job.lsf
#1730323785
bjobs
#1730324400
g++ -fopenmp mf.cpp -o my_job
#1730324416
bsub < OpenMP_job.lsf
#1730324422
ls
#1730324436
cat table_of_iter_and_time_4.txt 
#1730325373
ls
#1730325381
rm *.txt
#1730325382
ls
#1730325401
rm my*
#1730325402
ls
#1730325418
g++ -fopenmp mf.cpp -o my_job
#1730325423
bsub < OpenMP_job.lsf
#1730325431
bjobs
#1730325678

#1730325684
vim OpenMP_job.lsf 
#1730325696
bsub < OpenMP_job.lsf
#1730325700
bjobs
#1730326397
ды
#1730326400
ls
#1730326407
cat table_of_iter_and_time_1.txt 
#1730326412
ls
#1730326416
ls -l
#1730326431
ls
#1730326473
bjobs
#1730326476
cat table_of_iter_and_time_1.txt 
#1730326547
bjobs
#1730326553
cat table_of_iter_and_time_2.txt 
#1730326567
ls -la
#1730326624
cat table_of_iter_and_time_1.txt 
#1730326704
g++ -fopenmp mf.cpp -o my_job
#1730326727
g++ -fopenmp mf.cpp -o my_job2
#1730326735
bsub < OpenMP_job.lsf
#1730326739
bjobs
#1730326747
ls
#1730326760
rm table_of_iter_and_time_1.txt 
#1730326762
ls
#1730326765
bjobs
#1730326768
ls
#1730326777
bjobs
#1730326789
bjobs -l 1185026
#1730326798
ls
#1730326808
bsub < OpenMP_job.lsf
#1730326813
ls
#1730326889
bjobs
#1730362803
ls
#1730362813
cat table_of_iter_and_time_1.txt 
#1730362824
ls
#1730362912
cat OpenMP_job.lsf
#1730388331
q
#1730388333
quit
#1730388342
exit
#1730374243
ls
#1730374285
g++ -fopenmp mf.cpp -o my_job2
#1730374295
bsub < OpenMP_job.lsf
#1730374299
cat OpenMP_job.lsf 
#1730375494
bjobs
#1730375496
ls
#1730375503
cat table_of_iter_and_time_4.txt 
#1730375555
ls
#1730375575
g++ -fopenmp mf.cpp -o my_job2
#1730375579
bsub < OpenMP_job.lsf
#1730375582
bjobs
#1730376379
ls
#1730376383
bjobs
#1730376391
cat table_of_iter_and_time_8.txt 
#1730376394
ls
#1730376413
g++ -fopenmp mf.cpp -o my_job2
#1730376416
bsub < OpenMP_job.lsf
#1730376419
bjobs
#1730376670
ls
#1730376674
bjobs
#1730376693
g++ -fopenmp mf.cpp -o my_job2
#1730376696
bsub < OpenMP_job.lsf
#1730377565
bjobs
#1730377569
ls
#1730377577
cat table_of_iter_and_time_32.txt 
#1730377645
g++  mf.cpp -o my_job2
#1730377655
g++ -fopenmp mf.cpp -o my_job2
#1730377660
bsub < OpenMP_job.lsf
#1730377663
bjobs
#1730380654
ls
#1730380661
ls -l
#1730380680
cat table_of_iter_and_time_1
#1730380686
cat table_of_iter_and_time_1.txt 
#1730390609
ls
#1730390624
cat table_of_iter_and_time_1.txt 
#1730391661
ls
#1730391669
g++ -fopenmp mf.cpp -o my_job2
#1730391673
bsub < OpenMP_job.lsf
#1730391677
bjobs
#1730397088
cat table_of_iter_and_time_1
#1730397091
cat table_of_iter_and_time_1.txt 
#1730397094
ls
#1730397104
g++ -fopenmp mf.cpp -o my_job2
#1730397110
bsub < OpenMP_job.lsf
#1730397112
ls
#1730397128
bjobs
#1730398124
g++ -fopenmp mf.cpp -o my_job2
#1730398129
ls
#1730398137
cat table_of_iter_and_time_1.txt 
#1730398143
bjobs
#1730398972
ls
#1730398980
cat table_of_iter_and_time_1.txt 
#1730398997
g++ -fopenmp mf.cpp -o my_job2
#1730399002
bsub < OpenMP_job.lsf
#1730399005
bjobs
#1730399034
bjobs -all -a | grep RUN
#1730399052
bjobs
#1730399158
ls
#1730399160
bjobs
#1730399461
g++ -fopenmp mf1_10_20.cpp -o my_job1_10_20
#1730399473
cp OpenMP_job.lsf OpenMP_job_1_10_20.lsf 
#1730399481
vim OpenMP_job_1_10_20.lsf 
#1730399511
bsub < OpenMP_job_1_10_20.lsf 
#1730399514
bjobs
#1730399552
ls
#1730399617
bjobs
#1730399736
bjobs -l 1188314
#1730399758
bjobs
#1730399820
g++ -fopenmp mf_4_40.cpp -o my_job_4_40
#1730399841
cp OpenMP_job.lsf OpenMP_job_4_40.lsf
#1730399854
vim OpenMP_job_4_40.lsf 
#1730399873
bsub < OpenMP_job_4_40.lsf 
#1730399876
bjobs
#1730399883
ls
#1730399890
cat table_of_iter_and_time_1.txt
#1730399910
vim mf.cpp 
#1730399942
g++ -fopenmp mf.cpp -o my_job2
#1730399952
bsub < OpenMP_job.lsf
#1730399955
bjobs
#1730399960
ls
#1730399968
cat table_of_iter_and_time_1_10_20.txt 
#1730399980
bjobs
#1730400024
иощиы
#1730400027
bjobs
#1730400033
ls
#1730400041
cat table_of_iter_and_time_4_40.txt 
#1730400099
cp OpenMP_job.lsf OpenMP_job_16_40.lsf
#1730400108
vim OpenMP_job_16_40.lsf 
#1730400122
ls
#1730400155
g++ -fopenmp mf_16_40.cpp -o my_job_16_40
#1730400170
bsub < OpenMP_job_16_40.lsf 
#1730400173
bjobs
#1730400178
ls
#1730400186
bjobs
#1730400190
ls
#1730400228
bjobs
#1730400253
cat table_of_iter_and_time_16_40.txt 
#1730400335
bjobs
#1730400365
vim OpenMP_job
#1730400371
vim OpenMP_job.lsf 
#1730400405
bjobs
#1730400408
ls
#1730400424
cat table_of_iter_and_time_1.txt 
#1730400938
ls
#1730400941
cat table_of_iter_and_time_1.txt 
#1730401041
ls *.txt
#1730401083
bjobs
#1730401139
bjobs -u all -a | grep RUN
#1730401246
ls *.txt
#1732736147
ls
#1732736153
mkdir task2
#1732736157
quit
#1732736164
exit
#1732736291
ls
#1732736294
cd task2/
#1732736295
ls
#1732736300
cat mf_mpi.cpp 
#1732818395
cd task2/
#1732818396
ls
#1732818526
mpicxx -o task2 mf_mpi.cpp
#1732818770
mpixlC -o task2 mf_mpi.cpp
#1732818846
mpicxx --мукышщт
#1732818853
mpicxx --version
#1732818947
module load SpectrumMPI
#1732818951
mpicxx --version
#1732818968
module load OpenMPI
#1732819240
mpicxx -o task2 mf_mpi.cpp
#1732819493
vim MPI_JOB.lsf
#1732819676
bsub < MPI_JOB.lsf ./task2 40 40
#1732819680
иощиы
#1732819683
bjobs
#1732819690
ls
#1732819736
bjobs
#1732819740
vim MPI_JOB.lsf 
#1732820196
ls
#1732820198
vim MPI_JOB.lsf 
#1732820216
bsub < MPI_JOB.lsf 
#1732820219
bjobs
#1732820222
ls
#1732820231
cat res_40_40.1237002.out
#1732820257
ls
#1732821521
rm res_40_40.1237002.*
#1732821522
ls
#1732821529
cat mf_mpi.cpp 
#1732821564
ls
#1732821569
rm mf_mpi.cpp 
#1732821572
rm task2 
#1732821681
ls
#1732821833
mpicxx -o task2 mf_mpi.cpp
#1732821852
vim MPI_JOB.lsf 
#1732821915
bsub < MPI_JOB.lsf 
#1732821917
ls
#1732821922
bjobs
#1732822015
bjobs -u all -a | grep RUN
#1732822061
bjobs
#1732822072
bjobs -u all -a | grep RUN
#1732822098
bjobs
#1732822187
ls
#1732822192
cat res_40_40.1
#1732822198
cat res_40_40.1.out 
#1732822223
vim MPI_JOB.lsf 
#1732822276
bsub < MPI_JOB.lsf 
#1732822284
bjobs
#1732822290
bjobs -u all -a | grep RUN
#1732822332
bjobs
#1732822335
ls
#1732822346
Ccat res_40_40.2.out 
#1732822351
cat res_40_40.2.out 
#1732822367
vim MPI_JOB.lsf 
#1732822383
bsub < MPI_JOB.lsf 
#1732822587
bjobs
#1733513195
ls
#1733513203
cd task2/
#1733513206
ls
#1733513213
cd ..
#1733513219
mkdir task3
#1733513230
cd task3/
#1733513232
ls
#1733513266
module load SpectrumMPI
#1733513295
module load OpenMPI
#1733513328
mpicxx -o task2 mf_mpi.cpp -fopenmp
#1733513330
ls
#1733513341
mpicxx -o task2 mf_mpi.cpp 
#1733513349
mpicxx -o task2 mf_mpi.cpp -fopenmp
#1733513634
vim OMP_MPI.lsf 
#1733514056
bsub < OMP_MPI.lsf 
#1733514059
bjobs
#1733514062
ls
#1733514070
cat my_job.40_40_1_1.out 
#1733514161
mpicxx -o task2 mf_mpi.cpp -fopenmp
#1733514169
bsub < OMP_MPI.lsf 
#1733514192
ls
#1733514197
bjobs
#1733514205
cat my_job.40_40_1_1.out 
#1733514218
ls
#1733514226
cat my_job.40_40_1_1.err 
#1733514282
bjobs
#1733514289
bjobs -u all -a | grep RUN
#1733514293
bjobs -u all -a 
#1733514296
bjobs -u all -a | grep RUN
#1733514306
ls
#1733514314
rm my_job.40_40_1_1.*
#1733514321
bsub < OMP_MPI.lsf 
#1733514324
bjobs
#1733514334
cat my_job.40_40_1_1.out
#1733514349
cat my_job.40_40_1_1.err 
#1733514373
vim OMP_MPI.lsf 
#1733514379
ls
#1733514381
vim OMP_MPI.lsf 
#1733514408
bsub < OMP_MPI.lsf 
#1733514412
bjobs
#1733514428
ls
#1733514433
cat my_job.40_40_1_1.out 
#1733514444
cat my_job.40_40_1_1.err 
#1733514501
vim OMP_MPI.lsf 
#1733514514
bsub < OMP_MPI.lsf 
#1733514517
bjobs
#1733514845
ls
#1733514849
cat my_job.40_40_1_1.out 
#1733515630
ls
#1733515636
vim OMP_MPI.lsf 
#1733515673
bsub < OMP_MPI.lsf 
#1733515676
bjobs
#1733515746
cat my_job.40_40_2_1.out 
#1733515756
cat my_job.40_40_1_1.out 
#1733515810
vim OMP_MPI.lsf 
#1733515839
bsub < OMP_MPI.lsf 
#1733515847
bjobs
#1733516277
ls
#1733516282
cat my_job.40_40_seq.
#1733516285
cat my_job.40_40_seq.out 
#1733516819
vim OMP_MPI.lsf 
#1733516852
bsub < OMP_MPI.lsf 
#1733516856
bjobs
#1733518147
cat OMP_MPI.lsf 
#1733518260
bjobs
#1733521203
ls
#1733521228
vim OMP_MPI.lsf 
#1733521256
bsub < OMP_MPI.lsf 
#1733521258
bjobs
#1733523623
vim OMP_MPI.lsf 
#1733523645
bsub < OMP_MPI.lsf 
#1733523647
bjobs
#1733524398
bsub < OMP_MPI.lsf 
#1733524405
bjobs
#1733572355
cd task3/
#1733572356
ls
#1733572362
vim OMP_MPI.lsf 
#1733572382
bsub < OMP_MPI.lsf 
#1733572385
bjobs
#1733572402
cat my_job.80_90_2_4.out 
#1733572445
cat my_job.80_90_2_2.out 
#1733572454
cat my_job.80_90_2_1.out 
#1733572484
ls
#1733572567
module load SpectrumMPI
#1733572574
module load OpenMPI
#1733572585
mpicxx -o task2 mf_mpi.cpp -fopenmp
#1733572589
cat Op
#1733572591
cat OMP_MPI.lsf 
#1733572601
bsub < OMP_MPI.lsf 
#1733572604
bjobs
#1733572620
cat my_job.40_40_1_1.out 
#1733572626
cat my_job.40_40_2_1.out 
#1733572630
cat my_job.40_40_2_2.out 
#1733572651
cat my_job.40_40_seq.out 
#1733572667
cat my_job.80_90_seq.out 
#1733572678
cat OMP_MPI.lsf 
#1733572680
bjobs
#1733573025
cat OMP_MPI.lsf 
#1733573034
cat my_job.80_90_2_4.out
#1733573042
vim OMP_MPI.lsf 
#1733573076
bsub < OMP_MPI.lsf 
#1733573078
bjobs
#1733573859
mpicxx -o task2 mf_mpi.cpp -fopenmp
#1733573862
bjobs
#1733573868
vim OMP_MPI.lsf 
#1733573901
bsub < OMP_MPI.lsf 
#1733573904
bjobs
#1733573920
cat OMP_MPI.lsf 
#1733573930
cat my_job.80_90_1_1_real.out
#1733573955
mpicxx -o task2 mf_mpi.cpp -fopenmp
#1733573961
bsub < OMP_MPI.lsf 
#1733573964
bjobs
#1733573975
cat my_job.80_90_1_1_real.out
#1733573992
mpicxx -o task2 mf_mpi.cpp -fopenmp
#1733573994
bsub < OMP_MPI.lsf 
#1733573998
bjobs
#1733574012
cat my_job.80_90_1_1.out
#1733574019
cat my_job.80_90_1_1_real.out
#1733574023
bjobs
#1733574225
mpicxx -o task2 mf_mpi.cpp -fopenmp
#1733574247
vim OMP_MPI.lsf 
#1733574264
bsub < OMP_MPI.lsf 
#1733574267
bjobs
#1733574276
cat my_job.80_90_1_1_real.out
#1733574296
mpicxx -o task2 mf_mpi.cpp -fopenmp
#1733574299
vim OMP_MPI.lsf 
#1733574325
bsub < OMP_MPI.lsf 
#1733574328
bjobs
#1733574406
cat OMP_MPI.lsf 
#1733574416
cat my_job.80_90_2_1_real.out
#1733574423
vim OMP_MPI.lsf 
#1733574436
bsub < OMP_MPI.lsf 
#1733574438
bjobs
#1733574660
cat OMP_MPI.lsf 
#1733574667
cat my_job.80_90_2_2_real.out
#1733574673
vim OMP_MPI.lsf 
#1733574694
bsub < OMP_MPI.lsf 
#1733574697
bjobs
#1733574850
cat OMP_MPI.lsf 
#1733574856
cat my_job.80_90_2_4_real.out
#1733574861
vim OMP_MPI.lsf 
#1733574893
bsub < OMP_MPI.lsf 
#1733574895
bjobs
#1733575136
cat OMP_MPI.lsf 
#1733575143
cat my_job.80_90_2_8_real.out
#1733575194
vim OMP_MPI.lsf 
#1733575234
bsub < OMP_MPI.lsf 
#1733575236
bjobs
#1733575309
cat OMP_MPI.lsf 
#1733575315
cat my_job.160_180_1_1_real.out
#1733575330
cat mf_mpi.cpp 
#1733575346
bjobs
#1733575365
mpicxx -o task2 mf_mpi.cpp -fopenmp
#1733575369
bsub < OMP_MPI.lsf 
#1733575372
bjobs
#1733575486
cat my_job.160_180_1_1_real.out
#1733575503
mpicxx -o task2 mf_mpi.cpp -fopenmp
#1733575505
bsub < OMP_MPI.lsf 
#1733576536
bjobs
#1733576541
cat OMP_MPI.lsf 
#1733576547
cat my_job.160_180_1_1_real.out
#1733576564
mpicxx -o task2 mf_mpi.cpp -fopenmp
#1733576567
bsub < OMP_MPI.lsf 
#1733576570
bjobs
#1733576578
cat my_job.160_180_1_1_real.out
#1733576591
mpicxx -o task2 mf_mpi.cpp -fopenmp
#1733576593
bsub < OMP_MPI.lsf 
#1733576596
bjobs
#1733577113
cat OMP_MPI.lsf 
#1733577122
cat my_job.160_180_1_1_real.out
#1733577128
vim OMP_MPI.lsf 
#1733577153
bsub < OMP_MPI.lsf 
#1733577155
bjobs
#1733577216
cat my_job.160_180_1_1_real.out
#1733577220
cat OMP_MPI.lsf 
#1733577226
cat "my_job.160_180_4_1_real.out


#1733577237
cat my_job.160_180_4_1_real.out
#1733577242
vim OMP_MPI.lsf 
#1733577260
bsub < OMP_MPI.lsf 
#1733577264
bjobs
#1733577861
cat OMP_MPI.lsf 
#1733577867
cat my_job.160_180_4_2_real.out
#1733577870
vim OMP_MPI.lsf 
#1733577884
bsub < OMP_MPI.lsf 
#1733578384
bjobs
#1733578387
cat OMP_MPI.lsf 
#1733578393
cat my_job.160_180_4_4_real.out
#1733578398
vim OMP_MPI.lsf 
#1733578418
bsub < OMP_MPI.lsf 
#1733578420
bjobs
#1733674741
ls
#1733689655
cd task3/
#1733689657
ls
#1733689741
ls *.out
#1733689788
ls *40*.out
#1733689795
cat my_job.40_40_seq.out
#1733689891
ls *40*.out
#1733689898
cat my_job.40_40_1_1.out
#1733690757
ls
#1733690780
ls *40*
#1733690785
ls *40*.out
#1733690836
cat my_job.40_40_seq.out
#1733690929
cat my_job.40_40_2_1.out 
#1733691095
cat my_job.4
#1733691204
ls *80*.out
#1733691218
ls *80*real*.out
#1733691225
ls *90*real*.out
#1733691243
cat my_job.80_90_1_1_real.out
#1733691309
cat my_job.80_90_2_1_real.out
#1733691746
ls *160*real*.out
#1733691770
cat my_job.160_180_4_2_real.out
#1733692009
ls *90*real*.out
#1733692031
cat my_job.80_90_1_1_real.out
#1733692046
ls *90*real*.out
#1733692055
cat my_job.80_90_2_1_real.out
#1733692072
ls *90*real*.out
#1733692086
cat  my_job.80_90_2_2_real.out
#1733692115
ls *90*real*.out
#1733692126
cat  my_job.80_90_2_4_real.out
#1733692137
ls *90*real*.out
#1733692152
cat my_job.80_90_2_8_real.out
#1733692191
ls *160*real*.out
#1733692199
cat my_job.160_180_1_1_real.out 
#1733692210
ls *160*real*.out
#1733692217
cat  my_job.160_180_4_1_real.out
#1733692669
ls *160*real*.out
#1733692683
cat  my_job.160_180_4_1_real.out
#1733692697
ls *160*real*.out
#1733692704
cat my_job.160_180_4_2_real.out
#1733692725
cat my_job.160_180_4_4_real.out
#1733692737
cat my_job.160_180_4_8_real.out
#1733695605
cd ..
#1733695612
ls
#1733695619
cat OpenMP_job_16_40.lsf 
#1733695625
vim OpenMP_job
#1733695645
vim OpenMP_job_16_40.lsf 
