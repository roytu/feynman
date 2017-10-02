
""" Tests for things """

import sympy as sy
from sympy import *
from feynman import Amplitude
from term import Momentum

def test(result, name):
    if result:
        print("{0} passed".format(name))
    else:
        print("{0} failed".format(name))

def run_flip_variants_test():
    p_up = Momentum("p", "SLASHmu", 1)
    p_down = Momentum("p", "SLASHmu", 0)
    m = sy.Symbol("m")

    expr = Amplitude.flip_variant(p_up + m)
    test(expr == p_down + m, "run_flip_variants_test_1")

    expr = Amplitude.flip_variant(p_down)
    test(expr == p_up, "run_flip_variants_test_2")

def run_subtract_momenta_test():
    p = Momentum("p", "SLASHmu", 1)
    k = Momentum("k", "SLASHmu", 1)

    expr = p - k

    test(expr != sy.S.Zero, "run_subtract_momenta_test")

def run_highest_log_test():
    Lamb = sy.Symbol("Lamb")

    expr = 2 * sy.log(Lamb) + 1
    #expr = expr.subs(Lamb, sy.exp(Lamb))
    #expr = sy.O(expr, (Lamb, sy.oo)).args[0]
    import pdb; pdb.set_trace()
    expr = sy.O(expr, (Lamb, sy.oo)).args[0]
    #expr = expr.subs(Lamb, sy.log(Lamb))
    print(expr)

    test(expr == 2 * sy.log(Lamb), "run_highest_log_test")

if __name__ == "__main__":
    run_flip_variants_test()
    run_subtract_momenta_test()
    run_highest_log_test()
