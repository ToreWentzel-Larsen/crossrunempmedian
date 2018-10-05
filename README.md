# crossrunempmedian
Joint distribution of number of crossings and longest run with empirical median.

In statistical quality control the midline is in many cases the median in the sample itself. Then the calculations in the package crossrun do not apply. This case is addressed here. The calculations so far used have high complexity and are only practically feasible up to n=28. Work is ongoing to find another procedure,  closer to the basic idea behind the iterative procedure in the main function crossrunbin in the package crossrun with only linearly higher complexity than the procedure in crossrunbin. It this succeeds the procedure will probably be imnplemented as a function in crossrun. 
