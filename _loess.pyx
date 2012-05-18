#cdef struct in_t:
#    long    n
#    long    p
#    double  *y
#    double  *x
#    double  *weights
#
#cdef struct model_t:
#    double  span
#    long    degree
#    long    normalize
#    long    parametric[8]
#    long    drop_square[8]
#    char    *family
#
#cdef struct control_t:
#    char    *surface
#    char    *statistics
#    double  cell
#    char    *trace_hat
#    long    iterations
#
#cdef struct kd_tree_t:
#    long	*parameter
#    long	*a
#    double	*xi
#    double	*vert
#    double	*vval
#
#cdef struct out_t:
#    double	*fitted_values
#    double  *fitted_residuals
#    double  enp
#    double	s
#    double  one_delta
#    double	two_delta
#    double	*pseudovalues
#    double	trace_hat
#    double	*diagonal
#    double	*robust
#    double  *divisor
#
#cdef struct loess_struct_t:
#    struct in_t 
#    struct model_t
#    struct control_t
#    struct kd_tree_t
#    struct out_t

cdef extern from "loess.h":
    struct loess_struct

cdef extern from "loess.c":
    void loess_setup(double *x, double *y, long n, long p, loess_struct *lo)

    void loess(loess_struct *lo)

def test_setup():
    n = 100
    p = 2
    one_two = array([])
    response = array([])
