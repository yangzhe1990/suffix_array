suffix_array
============

This program runs faster than sais. But I didn't minimize the memory comsumptions. I may try to compete kvark's implementation in future.

See compare.txt for the benchmark detail measured on a really old, P8600, single-channel 1066MHz memory laptop. i686|x86.benchmark.* shows the benchmark measured on a i7-3517U single-channel 1600MHz memory laptop. And the result shows that ka_reduced performs the best!

C.f. http://yangzhe1990.wordpress.com/2012/11/04/ka-and-induced-sorting-which-is-faster/ for the detail of the algorithm.
