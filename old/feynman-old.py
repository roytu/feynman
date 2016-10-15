
""" Given a Feynman diagram, perform integrals. """

from term import Term, Frac, Integral, I, C, Sum, Sub, Super, Gamma, Prod
from latex import Latex

def fermion_propagator():
    """ Calculate expression for first-order corrected fermion propagator. """

    latex = Latex()

    # 1. Write out amplitude  TODO include regulator terms
    m = []

    ind = "\sigma_1"
    sum_terms = [Sub("p", ind), C("m")]
    m.append(Frac(Prod([Sum(sum_terms), Gamma(ind)]), Term("p^2 - m^2 + i\\epsilon")))

    m.append(Integral("k", norm=True))
    m.append(I())
    m.append(C("e"))
    m.append(Term("\\gamma^\\mu"))

    ind = "\sigma_2"
    sum_terms = [Sub("p", ind), Sub("(-k)", ind), C("m")]
    m.append(Frac(Prod([Sum(sum_terms), Gamma(ind)]), Term("(p - k)^2 - m^2 + i\\epsilon")))
    m.append(I())
    m.append(C("e"))
    m.append(Term("\\gamma^\\nu"))
    m.append(Term("D_{\\mu \\nu}(k)"))

    ind = "\sigma_3"
    sum_terms = [Sub("p", ind), C("m")]
    m.append(Frac(Prod([Sum(sum_terms), Gamma(ind)]), Term("p^2 - m^2 + i\\epsilon")))

    latex.add(" ".join([t.latex() for t in m]))

    # 2. Use Feynman parameterization to rewrite propagators
    m_ = []
    ds = []
    for term in m:
        if term.is_fraction():
            # TODO
            m_.append(term.numer)
            ds.append(term.denom)
        else:
            m_.append(term)

    # Integrals
    for i, d in enumerate(ds):
        a = Term("0")
        b = Term("1" + "".join([" - z_{0}".format(j + 1) for j in range(i)]))
        m_.append(Integral("z_{0}".format(i+1), norm=True, limits=(a, b)))

    # Frac
    denom = ""
    denom += "\\left["
    ds_ = []
    for i, d in enumerate(ds):
        ds_.append("\\left(" + d.latex() + "\\right) " + "z_{0}".format(i + 1))
    denom += " + ".join(ds_)
    denom += "\\right]"
    denom += "^{0}".format(len(ds))
    frac = Frac(Term("1"), Term(denom))
    m_.append(frac)

    m = m_
    latex.add(" ".join([t.latex() for t in m]))

    # 3. Distribute numerators (p + m) ish terms
    # 4. Solve each internal momenta integral
    # 5. Move terms independent of z's to left
    # 6. Groan.  Solve.
    # 7. Solve remaining integrals, contract metrics, etc.
    # 8. (OPTIONAL) Evaluate traces for non-Abelian theories.

    latex.render()

if __name__ == "__main__":
    fermion_propagator()
