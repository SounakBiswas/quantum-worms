#! /bin/sh
>op_temp.dat
>op_temp1.dat
>op_maxerr.dat
>op_in.dat


op=$1 #op is the column of the op for which binder is desired.
f_no=$2 #f_no is the file which contains the desired colunm


f1=binop
f2=binlocal

###################################################


for h in 1.0000 0.8000
do

for t in 1.0000 0.8000
do

for lat in 3
do


awk -v var1=$op '{print $var1}' $(eval "echo \$f"$f_no)L${lat}j1.0000t${t}h${h}  > op_in.dat

length=$(eval "wc -l < op_in.dat")

nbins=$((length))

sed -e "s/#NBINS#/$nbins/g" op_maxerr_gen.c > op_maxerr.c
##sed -e "s/#BINSIZE#/$binsize/g" temp.c > opcalc.c

#rm temp.c

gcc op_maxerr.c -lm

./a.out

rm op_maxerr.c

cat op_out.dat >> op_temp.dat

awk -v var1=$lat -v var2=$t -v var3=$h '{print var1,var2,var3,$1,$2}' op_temp.dat >> op_temp1.dat

rm op_in.dat
rm op_out.dat
rm op_temp.dat
done

cat op_temp1.dat >> op_temp2.dat
echo "" >> op_temp2.dat
echo "" >> op_temp2.dat
rm op_temp1.dat
done

cat op_temp2.dat >> op_maxerr.dat
#echo "" >> op_maxerr.dat
#echo "" >> op_maxerr.dat
rm op_temp2.dat
done




#cat << END_OF_PROGRAM > pl.p
#set terminal postscript eps size 3.5,2.62 enhanced color font 'Helvetica,20' linewidth 2
#set key title "binder_cumulant"
#set xlabel "no of bins"
#set ylabel "binder_c2dotc2"
#set output "binder.eps"
#pl "binder.dat" u 0:1:2 w errorlines title "L32", "binder.dat" u 0:3:4 w errorlines title "L48","binder.dat" u 0:5:6 w errorlines title "L64"
#END_OF_PROGRAM
#gnuplot pl.p
#convert -flatten -density 500 binder.eps binder.png
#rm foo.dat pl.p a.out

rm a.out

#echo 'done'


