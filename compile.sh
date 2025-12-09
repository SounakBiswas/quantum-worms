    job=L3j1.0000t1.0000h1.0000
gcc -O -o $job.out -fcommon link_list.c set_constant_values.c make_lattice.c do_initialisation.c mt19937ar.c nrutil.c do_updates.c bin.c measure.c save.c main.c autocorr.c -lm -g
