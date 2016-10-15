
""" Given a Feynman diagram, perform integrals. """

from math import factorial
from term import Term, Frac, Integral, I, C, Sum, Sub, Super, Gamma, Prod, Metric
from latex import Latex

def fermion_propagator():
    """ Calculate expression for first-order corrected fermion propagator. """

    latex = Latex()

    # 1. Write out amplitude  TODO include regulator terms
    m = []

    phase = 0   # Power of i
    const = []  # List of constants (as a product)
    metrics = []  # List of metrics (as a product)
    integrals = []  # List of Integral objects
    numers = []  # List of Terms in numerator (as a sum of products)
    denoms = []  # List of Terms in denominator

    def latex_add():
        s = ""
        s += ["", "i", "-", "-i"][phase % 4]  # Phase
        s += " ".join([t.latex() for t in const])  # Const
        s += " ".join([t.latex() for t in metrics])  # Metrics
        s += " ".join([t.latex() for t in integrals])  # Integrals
        n = " + ".join([" ".join([t.latex() for t in p]) for p in numers])
        d = " ".join(["(" + t.latex() + ")" for t in denoms])  # Denoms
        s += "\\frac{{ 1 }}{{ {0} }}".format(d)
        s += "\\left( {0} \\right)".format(n)
        latex.add(s)

    def sop_mult(sop1, sop2):
        """ Multiply sum of products """
        if sop1 == []:
            return multiplicand
        else:
            lst = []
            for prod1 in sop1:
                for prod2 in sop2:
                    lst.append(prod1 + prod2)
            return lst

    e = C("e")

    # 1. Write out amplitude

    # S_F(p)
    phase += 1
    ind = "\\sigma_1"

    multiplicand = [[Gamma(ind), Sub("p", ind)], [C("m")]]
    numers = sop_mult(numers, multiplicand)
    denoms.append(Term("p^2 - m^2 + i\\epsilon"))

    # int dk
    integrals.append(Integral("k", norm=True))

    # ie gamma mu
    ind = "\\mu"
    phase += 1
    const.append(e)
    multiplicand = [[Gamma(ind)]]
    numers = sop_mult(numers, multiplicand)

    # S_F(p - k)
    phase += 1
    ind = "\\sigma_2"
    multiplicand = [[Gamma(ind), Sub("p", ind)], [Gamma(ind), Sub("-k", ind)], [C("m")]]
    numers = sop_mult(numers, multiplicand)
    denoms.append(Term("(p - k)^2 - m^2 + i\\epsilon"))

    # ie gamma nu
    ind = "\\nu"
    phase += 1
    const.append(e)
    multiplicand = [[Gamma(ind)]]
    numers = sop_mult(numers, multiplicand)

    # S_F(p)
    phase += 1
    ind = "\\sigma_1"
    multiplicand = [[Gamma(ind), Sub("p", ind)], [C("m")]]
    numers = sop_mult(numers, multiplicand)
    denoms.append(Term("p^2 - m^2 + i\\epsilon"))

    # D_{\mu, \nu}(k)
    phase += 3
    metrics.append(Metric("\\mu", "\\nu"))
    denoms.append(Term("k^2 + i\\epsilon"))

    latex_add()

    # Parse denominators
    n = len(denoms)
    const.append(Term(str(factorial(n - 1))))
    for i in range(n):
        a = Term("0")
        b = Term("1" + "".join([" - z_{0}".format(j+1) for j in range(i)]))
        integrals.append(Integral("z_{0}".format(i+1), norm=True, limits=(a, b)))
    denoms = [Term("\\left[" + " + ".join(["({0})z_{1}".format(d.latex(), i+1) for i, d in enumerate(denoms)]) + "\\right]^{0}".format(n))]

    latex_add()

    latex.render()

if __name__ == "__main__":
    fermion_propagator()
