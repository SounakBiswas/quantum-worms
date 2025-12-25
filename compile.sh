    job=L3j1.0000t1.0000h1.0000
#gcc -O -o $job.out -fcommon link_list.c set_constant_values.c make_lattice.c do_initialisation.c mt19937ar.c  do_updates.c bin.c measure.c main.c autocorr.c utils_dual_graph.c construct_dual.c -lm -g
gcc -O -o $job.out -fcommon link_list.c set_constant_values.c make_lattice.c do_initialisation.c mt19937ar.c  do_updates.c bin_min.c measure_min.c main.c autocorr.c fat_utils.c construct_fat.c -lm -g
