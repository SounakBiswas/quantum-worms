#! /bin/sh

for t in 1.0000
#for t in 1.000
do
for h in 1.0000
#for h in 1.0000
do
for size in 3
do
for j in 1.0000
do


sed -e "s/#LX#/${size}/g" gen_global.h > foo1.h 
sed -e "s/#J#/${j}/g" foo1.h > foo2.h
sed -e "s/#T#/${t}/g" foo2.h > foo3.h
sed -e "s/#H#/${h}/g" foo3.h > global.h
rm foo1.h foo2.h foo3.h

sed -e "s/#LX#/${size}/g" gen_compile.sh > foo1.h 
sed -e "s/#J#/${j}/g" foo1.h > foo2.h
sed -e "s/#T#/${t}/g" foo2.h > foo3.h
sed -e "s/#H#/${h}/g" foo3.h > compile.sh
rm foo1.h foo2.h foo3.h 

chmod 777 compile.sh
./compile.sh

nohup ./L${size}j${j}t${t}h${h}.out &
#qsub L${size}j${j}t${t}h${h}.pbs


done
done
done
done

