    job=L3j1.0000t0.1000h0.8000
gcc -o $job.out autocorr.c link_list.c set_constant_values.c make_lattice.c do_initialisation.c mt19937ar.c nrutil.c do_updates.c bin_min.c measure_min.c main.c -lm -g
