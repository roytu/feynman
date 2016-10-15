
""" Given a Feynman diagram, perform integrals. """

from math import factorial
import sympy as sy
from latex import Latex
from term import Gamma, Momentum, Metric, U, UBar, Gamma0, Gamma1, Gamma2, Gamma3, MatrixTerm

from itertools import product

RENDER_ALL = False

class Amplitude(object):
    """ Amplitude is a data structure for holding M expressions. """
    def __init__(self):
        self.expr = 1
        self.numer = 1
        self.denom = 1
        self.denom_z = 1
        self.metrics = []
        self.integrals_zs = []  # Format is (z, a, b)
        self.integrals_internal = []
        self.integrals_cutoffs = []  # Cutoff integral parameters
        self.spinors = []
        self.inners = []  # Format is (constant, expr)
        self.indices = set()

    def U(self, p, spin):
        u = U(p, spin)
        self.spinors.append(u)
        return u

    def UBar(self, p, spin):
        u = UBar(p, spin)
        self.spinors.append(u)
        return u

    def V(self, e, ind):
        """ e: string
            ind: string
        """
        self.expr *= sy.I
        self.expr *= sy.Symbol(e)
        self.numer *= Gamma(ind)
        self.indices.add(ind)

    def S_F(self, p, m, p2, m2, ind):
        """ p: expression
            m: expression
            p2: expression
            m2: expression
            ind: string
        """
        self.expr *= sy.I
        self.numer *= Gamma(ind) * p + m
        self.denom *= p2 - m2
        self.indices.add(ind)

    def D_F(self, k2, ind1, ind2, t, lamb, Lamb):
        """ k2: expression
            ind1: string
            ind2: string
            t: string
            Lamb: symbol
            lamb: symbol
        """
        self.expr *= -sy.I
        self.metrics.append((ind1, ind2))

        # Regulator
        self.expr *= -1
        self.integrals_cutoffs.append((sy.Symbol(t), lamb ** 2, Lamb ** 2))
        self.denom *= (k2 - sy.Symbol(t)) ** 2

        self.indices.add(ind1)
        self.indices.add(ind2)

    def latex_add(self, latex):
        """ Numer / denom format """

        expr_ = self.expr * (self.numer / self.denom) * (1 / self.denom_z)
        for (a, b) in self.metrics:
            expr_ *= Metric(a, b)

        s = ""
        for (z, a, b) in self.integrals_zs + self.integrals_cutoffs:
            if a is not None and b is not None:
                s += "\\int\\limits_{{ {1} }}^{{ {2} }} d{0}".format(z, a, b)
            else:
                s += "\\int d{0}".format(z)

        for (k, a, b) in self.integrals_internal:
            if a is not None and b is not None:
                s += "\\int\\limits_{{ {1} }}^{{ {2} }} \\frac{{d^4 {0} }}{{ (2\pi)^4 }}".format(k, a, b)
            else:
                s += "\\int\\frac{{d^4 {0} }}{{ (2\pi)^4 }}".format(k)
        s += latex.get(expr_)
        latex.add(s)

    def latex_add2(self, latex):
        """ Inner format with no integrals """
        expr_ = self.expr * sum([c_ * e_ for (c_, e_) in self.inners])
        for (a, b) in self.metrics:
            expr_ *= Metric(a, b)

        s = ""
        s += latex.get(expr_)
        latex.add(s)

def fermion_propagator():
    """ Calculate expression for first-order corrected fermion propagator. """

    latex = Latex()

    ################################################
    ########     CONSTRUCT AMPLITUDE      ##########
    ################################################
    e = "e"

    p = "p"
    p2 = "{{ p^2 }}"
    k = "k"
    k2 = "{{ k^2 }}"
    pk = "{{ ( pk ) }}"
    m = "m"
    m2 = "{{ m^2 }}"
    tensors = ["\\sigma_2"]
    Lamb = "\\Lambda"
    lamb = "\\lambda"
    t = "t"
    mu = "\\mu"
    nu = "\\nu"

    spin1 = 0  # Up
    spin2 = 0  # Up

    k_to_k2 = {}
    k_to_k2[k] = k2

    amp = Amplitude()

    amp.integrals_internal.append((k, None, None))
    amp.expr /= (2 * sy.pi) ** 4
    amp.numer *= amp.UBar(p, spin1)  # ubar
    amp.V(e, mu)  # ie gamma mu
    ind = tensors[0]  # \\sigma_2
    amp.S_F(Momentum(p, ind) - Momentum(k, ind),
            sy.Symbol(m),
            sy.Symbol(p2) + sy.Symbol(k2) - 2 * sy.Symbol(pk),
            sy.Symbol(m2),
            ind)
    amp.V(e, nu)  # ie gamma nu
    amp.numer *= amp.U(p, spin2)  # u
    amp.D_F(sy.Symbol(k), mu, nu, t, sy.Symbol(lamb), sy.Symbol(Lamb))

    # Render
    latex.add_text("Raw amplitude")
    amp.latex_add(latex)
    if RENDER_ALL:
        latex.render()

    ################################################
    ########      FEYNMAN'S TRICK         ##########
    ################################################

    # Denominator kill
    latex.add_text("Feynman parameterization")

    denom_ = []
    for arg in amp.denom.args:
        if type(arg) == sy.Pow:
            base, power = arg.args
            for _ in range(power):
                denom_.append(base)
        else:
            denom_.append(arg)

    n = len(denom_)
    amp.expr *= factorial(n - 1)

    zs = [sy.Symbol("{{ z_{{ {0} }} }}".format(i+1)) for i in range(n)]
    amp.denom = sum([d * z for (d, z) in zip(denom_, zs)]).expand() ** n
    for i, z in enumerate(zs):
        a = 0
        b = 1 - sum(zs[:i])
        amp.integrals_zs.append((z, a, b))

    amp.latex_add(latex)
    if RENDER_ALL:
        latex.render()

    # Numerator expansion
    amp.numer = amp.numer.expand()

    if RENDER_ALL:
        latex.add_text("Numerator expansion")
        amp.latex_add(latex)
        latex.render()

    # TODO evaluate internal momenta integrals
    # At this point we stop with numer and denom and combine them into
    # one expression, `inner`, which is a sum of fractions.

    ################################################
    ########     EVAL. INTERNAL MOMENTA   ##########
    ################################################

    for (k, _, _) in amp.integrals_internal:  # k = "k"
        k2 = k_to_k2[k]
        denom_nopow, b = amp.denom.args[0], amp.denom.args[1]

        [C, D_] = sy.Poly(denom_nopow, sy.Symbol(k2)).coeffs()

        D = D_ / C
        amp.denom_z *= C ** b

        amp.denom = (sy.Symbol(k2) + D) ** b

        for term in amp.numer.args:  # term = k_{\sigma_2} m^2 gamma
            prod = term.args
            
            # Get k-vectors
            ks = [p for p in prod if isinstance(p, Momentum) and p.args[0] == k]

            if len(ks) % 2 == 1:  # Ward identity for odd tensors
                pass
            else:
                if len(ks) > 0:
                    # TODO convert higher-order even tensor integral to
                    #      scalar integral
                    pass
                else:
                    c, a = term.as_coeff_exponent(sy.Symbol(k2))
                    # Golden integral
                    c_ = sy.I * factorial(b - a - 3) * factorial(a + 1)
                    c_ /= factorial(b - 1)
                    c_ *= sy.pi ** 2

                    # TODO solve matrices in c
                    #matrices = []
                    #for arg in c.args:
                    #    if isinstance(arg, MatrixTerm):
                    #        matrices.append(arg)

                    # Resolve all matrices

                    c_ *= 1  # CHANGE THIS

                    c_new = 0
                    n = len(amp.indices)
                    for indices_vals in product([0, 1, 2, 3], repeat=n):
                        # Construct substitution dictionary
                        ind2val = {}
                        for (k, v) in zip(amp.indices, indices_vals):
                            ind2val[k] = v
                    
                        # Limit metrics
                        do_break = False
                        for (a_, b_) in amp.metrics:
                            if ind2val[a_] != ind2val[b_]:
                                do_break = True
                                break
                        if do_break:
                            continue

                        c__ = c
                        # Construct a dictionary of matrices / spinors
                        # and their replacements
                        replace_dict = {}
                        for (ind, val) in ind2val.items():
                            # When in doubt, add more underscores!
                            def replace_helper(expr):
                                temp = expr.resolve(val)
                                # If has explicit representation, use that instead
                                print("temp: {0}".format(temp))
                                if callable(getattr(temp, "explicit", None)):
                                    temp = temp.explicit()

                                print("returning: {0}".format(temp))
                                return temp

                            # oh my god fuck this library
                            # stupid hack to fix the stupid matrix bugs

                            # replace the index aw fuck this
                            print("ind: {0}".format(ind))
                            print("c__ before: {0}".format(c__))
                            c__, m = c__.replace(lambda expr: ind in expr.args,
                                              lambda expr: expr.resolve(val),
                                              simultaneous=True, exact=True, map=True)
                            print("map: {0}".format(m))
                            print("c__ after: {0}".format(c__))

                            # please let this nightmare end

                            # Gets a dictionary of things that would normally be
                            # replaced
                            f = lambda expr: ind in expr.args
                            #id_ = lambda expr: expr
                            id_ = lambda expr: sy.Dummy()
                            _, repmap = c__.replace(f, id_, map=True, exact=True)

                            for rep in repmap.keys():
                                replace_dict[rep] = replace_helper(rep)

                        # Substitution for spinors, etc
                        f = lambda e: callable(getattr(e, "explicit", None))
                        id_ = lambda expr: sy.Dummy()
                        _, repmap = c__.replace(f, id_, map=True, exact=True)

                        for rep in repmap.keys():
                            replace_dict[rep] = rep.explicit()

                        # Finally do the replace
                        import pdb; pdb.set_trace()
                        c__ = c__.xreplace(replace_dict)

                        c_new += c__

                    c_ *= c_new

                    expr_ = 1
                    expr_ /= D ** (b - a - 2)

                    amp.inners.append((c_, expr_))

    #latex_add2()
    # Compress inners by term
    inners_dict = {}
    for (c_, expr_) in amp.inners:
        if expr_ not in inners_dict:
            inners_dict[expr_] = 0
        inners_dict[expr_] += c_
    amp.inners = [(v, k) for (k, v) in inners_dict.items()]
    amp.latex_add2(latex)

    ################################################
    ########  CUTOFF AND Z INTEGRATIONS   ##########
    ################################################

    # Integrate cutoffs, then z integrals (backwards!)
    for (z, a, b) in amp.integrals_cutoffs + amp.integrals_zs[::-1]:
        subs_dict = {}

        inners_ = []
        for (c_, expr_) in amp.inners:
            def add_symbol(expr___):
                """ Adds a new symbol to `subs_dict` with the value
                    `expr` and returns it. """
                n = len(subs_dict)
                symb = sy.Symbol("C_{{ {0} }}".format(n))
                subs_dict[symb] = expr___
                return symb

            def simplify(expr__):
                """ Take a list of args and replace things with subs_dict. """
                if isinstance(expr__, sy.Add):
                    const = 0
                    args = []
                    # Check if guy has z
                    for arg in expr__.collect(z).args:
                        if z not in arg.free_symbols:
                            const += arg
                        else:
                            args.append(arg)
                    C = add_symbol(const)
                    res = 0
                    for arg in args:
                        res += simplify(arg)
                    res += C
                    return res
                elif isinstance(expr__, sy.Mul):
                    const = 1
                    args = []
                    # Check if guy has z
                    for arg in expr__.collect(z).args:
                        if z not in arg.free_symbols:
                            const *= arg
                        else:
                            args.append(arg)
                    res = 1
                    for arg in args:
                        res *= simplify(arg)

                    if const != 1:
                        C = add_symbol(const)
                        res *= C
                    return res
                else:
                    if expr__.args:
                        return (expr__.__class__)(*[simplify(arg) for arg in expr__.args])
                    return expr__

            latex.add_text("Starting...")

            # UV collect
            latex.add(latex.get(expr_))
            latex.add_text("UV collecting...")
            if RENDER_ALL:
                latex.render()

            uv = sy.Symbol(Lamb)
            if uv in expr_.free_symbols:
                expr_ = sy.O(expr_, (uv, sy.oo)).args[0]

            # Simplifying
            latex.add(latex.get(expr_))
            latex.add_text("Simplifying...")
            if RENDER_ALL:
                latex.render()

            expr_ = simplify(expr_)

            # Integrating
            latex.add(latex.get(expr_))
            latex.add_text("Integrating wrt ${0}$...".format(z))

            if RENDER_ALL:
                latex.render()

            expr_ = sy.integrate(expr_, (z, a, b))

            latex.add(latex.get(expr_))
            inners_.append((c_, expr_))

        amp.inners = inners_

        # Perform substitutions
        inners_ = []
        for (c_, expr_) in amp.inners:
            expr_ = expr_.subs(subs_dict)
            inners_.append((c_, expr_))
        amp.inners = inners_

        latex.add_text("Substituting constants...")
        amp.latex_add2(latex)

        latex.add_text("Simplifying...")
        expr_ = simplify(expr_)
        amp.latex_add2(latex)

        if RENDER_ALL:
            latex.render()

    ################################################
    ########  EVALUATE SPINS AND GAMMAS   ##########
    ################################################

    # TODO explicit spins (currently sum over all)

    #n = len(amp.indices)

    #inners_ = []
    #for (c_, expr_) in amp.inners:
    #    c_ = c_ * expr_
    #    expr_ = 1

    #    c_new = 0
    #    for indices_vals in product([0, 1, 2, 3], repeat=n):
    #        # Construct substitution dictionary
    #        ind2val = {}
    #        for (k, v) in zip(amp.indices, indices_vals):
    #            ind2val[k] = v

    #        # Limit metrics
    #        do_break = False
    #        for (a, b) in amp.metrics:
    #            if ind2val[a] != ind2val[b]:
    #                do_break = True
    #                break
    #        if do_break:
    #            continue

    #        c__ = c_
    #        for (ind, val) in ind2val.items():
    #            c__ = c__.replace(lambda expr: ind in expr.args,
    #                              lambda expr: expr.resolve(val))

    #        # TODO explicity convert c__ into a matrix expression and evaluate it
    #        import pdb; pdb.set_trace()
    #        latex.add(latex.get(c__))
    #        c_new += c__

    #    inners_.append((c_new, expr_))
    #amp.inners = inners_

    ################################################
    ########            RENDER            ##########
    ################################################

    amp.latex_add2(latex)
    latex.render()

if __name__ == "__main__":
    fermion_propagator()
