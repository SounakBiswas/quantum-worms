#! /bin/sh
#ulimit -s unlimited
####### energy using column 4 of file 1
./op_maxerr_l.sh 4 1 
cp op_maxerr.dat l_opcalc.dat

####### sigmax_q_int using column 5 of file 1

./op_maxerr_l.sh 5 1
awk '{print $4,$5}' op_maxerr.dat > temp.dat
paste -d' ' l_opcalc.dat temp.dat > temp2.dat 
mv temp2.dat l_opcalc.dat
rm temp.dat

####### sigmaz_q_int

./op_maxerr_l.sh 6 1
awk '{print $4,$5}' op_maxerr.dat > temp.dat
paste -d' ' l_opcalc.dat temp.dat > temp2.dat 
mv temp2.dat l_opcalc.dat
rm temp.dat

####### sigmaz_q_eq

./op_maxerr_l.sh 7 1
awk '{print $4,$5}' op_maxerr.dat > temp.dat
paste -d' ' l_opcalc.dat temp.dat > temp2.dat 
mv temp2.dat l_opcalc.dat
rm temp.dat

####### sigmaz_0_int

./op_maxerr_l.sh 8 1
awk '{print $4,$5}' op_maxerr.dat > temp.dat
paste -d' ' l_opcalc.dat temp.dat > temp2.dat 
mv temp2.dat l_opcalc.dat
rm temp.dat

####### sigmaz_0_eq

./op_maxerr_l.sh 9 1
awk '{print $4,$5}' op_maxerr.dat > temp.dat
paste -d' ' l_opcalc.dat temp.dat > temp2.dat 
mv temp2.dat l_opcalc.dat
rm temp.dat

####### 

./op_maxerr_l.sh 4 2
awk '{print $4,$5}' op_maxerr.dat > temp.dat
paste -d' ' l_opcalc.dat temp.dat > temp2.dat 
mv temp2.dat l_opcalc.dat
rm temp.dat

####### 

./op_maxerr_l.sh 6 2
awk '{print $4,$5}' op_maxerr.dat > temp.dat
paste -d' ' l_opcalc.dat temp.dat > temp2.dat 
mv temp2.dat l_opcalc.dat
rm temp.dat

####### 

./op_maxerr_l.sh 7 2
awk '{print $4,$5}' op_maxerr.dat > temp.dat
paste -d' ' l_opcalc.dat temp.dat > temp2.dat 
mv temp2.dat l_opcalc.dat
rm temp.dat

####### 

./op_maxerr_l.sh 9 2
awk '{print $4,$5}' op_maxerr.dat > temp.dat
paste -d' ' l_opcalc.dat temp.dat > temp2.dat 
mv temp2.dat l_opcalc.dat
rm temp.dat

####### 

./op_maxerr_l.sh 10 2
awk '{print $4,$5}' op_maxerr.dat > temp.dat
paste -d' ' l_opcalc.dat temp.dat > temp2.dat 
mv temp2.dat l_opcalc.dat
rm temp.dat

####### 

./op_maxerr_l.sh 12 2
awk '{print $4,$5}' op_maxerr.dat > temp.dat
paste -d' ' l_opcalc.dat temp.dat > temp2.dat 
mv temp2.dat l_opcalc.dat
rm temp.dat

####### 

./op_maxerr_l.sh 13 2
awk '{print $4,$5}' op_maxerr.dat > temp.dat
paste -d' ' l_opcalc.dat temp.dat > temp2.dat 
mv temp2.dat l_opcalc.dat
rm temp.dat

####### 

./op_maxerr_l.sh 15 2
awk '{print $4,$5}' op_maxerr.dat > temp.dat
paste -d' ' l_opcalc.dat temp.dat > temp2.dat 
mv temp2.dat l_opcalc.dat
rm temp.dat

########## the end
rm op_maxerr.dat
rm binder_jackerr.dat
