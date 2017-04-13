
""" Given a Feynman diagram, perform integrals. """

from math import factorial
import sympy as sy
from latex import Latex
from term import GammaFactory, UFactory, UBarFactory, MomentumFactory, \
        MetricFactory, Momentum, Metric, MatrixTerm

from itertools import product
from Queue import Queue

from util import find_first

RENDER_ALL = True

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

    def U(self, p):
        u = UFactory(p)
        self.spinors.append(u)
        return u

    def UBar(self, p):
        u = UBarFactory(p)
        self.spinors.append(u)
        return u

    def V(self, e, ind):
        """ e: string
            ind: string
        """
        self.expr *= sy.I
        self.expr *= sy.Symbol(e)
        self.numer *= GammaFactory(ind)
        self.indices.add(ind)

    def S_F(self, p, m, p2, m2, ind):
        """ p: expression
            m: expression
            p2: expression
            m2: expression
            ind: string
        """
        self.expr *= sy.I
        self.numer *= p * GammaFactory(ind) + m
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
        return  # TODO fuck this

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
            expr_ *= MetricFactory(a, b)

        s = ""
        s += latex.get(expr_)
        latex.add(s)

    def feynmans_trick(self, latex):
        # Denominator kill
        latex.add_text("Feynman parameterization")

        denom_ = []
        for arg in self.denom.args:
            if type(arg) == sy.Pow:
                base, power = arg.args
                for _ in range(power):
                    denom_.append(base)
            else:
                denom_.append(arg)

        n = len(denom_)
        self.expr *= factorial(n - 1)

        zs = [sy.Symbol("{{ z_{{ {0} }} }}".format(i+1)) for i in range(n)]
        self.denom = sum([d * z for (d, z) in zip(denom_, zs)]).expand() ** n
        for i, z in enumerate(zs):
            a = 0
            b = 1 - sum(zs[:i])
            self.integrals_zs.append((z, a, b))

        self.latex_add(latex)
        if RENDER_ALL:
            latex.render()

        # Numerator expansion
        self.numer = self.numer.expand()

        if RENDER_ALL:
            latex.add_text("Numerator expansion")
            self.latex_add(latex)
            latex.render()

    def internal_momenta(self, k_to_k2, latex):
        for (k, _, _) in self.integrals_internal:  # k = "k"
            k2 = k_to_k2[k]
            denom_nopow, b = self.denom.args[0], self.denom.args[1]

            [C, D_] = sy.Poly(denom_nopow, sy.Symbol(k2)).coeffs()

            D = D_ / C
            self.denom_z *= C ** b

            self.denom = (sy.Symbol(k2) + D) ** b

            self.latex_add(latex)
            latex.render()
            for term in self.numer.args:  # term = k_{\sigma_2} m^2 gamma
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
                        c_ *= c

                        expr_ = 1
                        expr_ /= D ** (b - a - 2)

                        self.inners.append((c_, expr_))

        # Compress inners by term
        inners_dict = {}
        for (c_, expr_) in self.inners:
            if expr_ not in inners_dict:
                inners_dict[expr_] = 0
            inners_dict[expr_] += c_
        self.inners = [(v, k) for (k, v) in inners_dict.items()]
        self.latex_add2(latex)

    def cutoffs(self, latex, Lamb):
        # Integrate cutoffs, then z integrals (backwards!)
        for (z, a, b) in self.integrals_cutoffs + self.integrals_zs[::-1]:
            subs_dict = {}

            inners_ = []
            for (c_, expr_) in self.inners:
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

            self.inners = inners_

            # Perform substitutions
            inners_ = []
            for (c_, expr_) in self.inners:
                expr_ = expr_.subs(subs_dict)
                inners_.append((c_, expr_))
            self.inners = inners_

            latex.add_text("Substituting constants...")
            self.latex_add2(latex)

            latex.add_text("Simplifying...")
            expr_ = simplify(expr_)
            self.latex_add2(latex)

            if RENDER_ALL:
                latex.render()

    @staticmethod
    def _split_into_matrices(term):
        """ Split a term such as:

            ie pi p_{sigma_2} bar{u}(p) gamma^mu gamma^{sigma_2} gamma^nu u(p)

            Returns:
                (const, list of metrics, list of momenta, list of matrix terms)
        """
        const = 1
        metrics = []
        momentas = []
        matrix_terms = []

        factors = []
        for prod in term.as_ordered_factors():
            # First write out all exponents
            if isinstance(prod, sy.Pow):
                base, exponent = prod.args
                for _ in range(exponent):
                    factors.append(base)
            else:
                factors.append(prod)

        # Distribute factors into correct lists
        for factor in factors:
            if isinstance(factor, Metric):
                metrics.append(factor)
            elif isinstance(factor, Momentum):
                momentas.append(factor)
            elif isinstance(factor, MatrixTerm):
                matrix_terms.append(factor)
            else:
                # Constant??
                const *= factor
        return (const, metrics, momentas, matrix_terms)

    @staticmethod
    def _merge_matrices(const, metrics, momentas, matrix_terms):
        """ Does the opposite of _split_into_matrices(term)

            Returns:
                term
        """
        result = const
        for factor in metrics + momentas + matrix_terms:
            result *= factor
        return result

    @staticmethod
    def _simplify_metrics(metrics):
        """ Takes a list of metrics and returns a list of metrics where
            terms like [(\mu, \nu) (\nu, \sigma)] -> [(\mu, \sigma)]
        """
        # TODO
        return metrics

    @staticmethod
    def _find_metric_contraction(A, B, metrics):
        """ Takes two gamma matrices and a list of metrics and checks if there's a hit.
            Returns the index of the metric, or None if no metric.

            Args:
                A: gamma matrix
                B: gamma matrix
                metrics: metrics list

            Returns:
                index (or None)
        """
        for i, metric in enumerate(metrics):
            metric1, metric2 = metric.args[:2]
            if (metric1, metric2) == (A_ind, B_ind) or (metric1, metric2) == (B_ind, A_ind):
                return i
        return None

    @staticmethod
    def _find_momenta(ind, momentas):
        """ Find the momentum associated with an index

            Args:
                ind: index
                momentas: list of momentum variables

            Returns:
                momentum
        """
        # TODO do
        pass

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

    k_to_k2 = {}
    k_to_k2[k] = k2

    p_to_m = {}
    p_to_m[p] = m

    amp = Amplitude()

    amp.integrals_internal.append((k, None, None))
    amp.expr /= (2 * sy.pi) ** 4
    amp.numer *= amp.UBar(p)  # ubar
    amp.V(e, mu)  # ie gamma mu
    ind = tensors[0]  # \\sigma_2
    amp.S_F(MomentumFactory(p, ind) - MomentumFactory(k, ind),
            sy.Symbol(m),
            sy.Symbol(p2) + sy.Symbol(k2) - 2 * sy.Symbol(pk),
            sy.Symbol(m2),
            ind)
    amp.V(e, nu)  # ie gamma nu
    amp.numer *= amp.U(p)  # u
    amp.D_F(sy.Symbol(k), mu, nu, t, sy.Symbol(lamb), sy.Symbol(Lamb))

    # Render
    latex.add_text("Raw amplitude")
    amp.latex_add(latex)
    if RENDER_ALL:
        latex.render()

    import pdb; pdb.set_trace()

    ################################################
    ########      FEYNMAN'S TRICK         ##########
    ################################################

    amp.feynmans_trick(latex)

    # TODO evaluate internal momenta integrals
    # At this point we stop with numer and denom and combine them into
    # one expression, `inner`, which is a sum of fractions.

    ################################################
    ########     EVAL. INTERNAL MOMENTA   ##########
    ################################################

    amp.internal_momenta(k_to_k2, latex)

    ################################################
    ########  CUTOFF AND Z INTEGRATIONS   ##########
    ################################################

    #amp.cutoffs(latex, Lamb)

    ################################################
    ########  EVALUATE SPINS AND GAMMAS   ##########
    ################################################

    # TODO explicit spins (currently sum over all)

    # Multiply metrics into terms
    #metrics = 1
    #for (ind1, ind2) in amp.metrics:
    #    metrics *= MetricFactory(ind1, ind2)
    #amp.metrics = []

    ## terms is a queue now
    #terms = Queue()
    #for arg in amp.numer.args:
    #    terms.put(arg * metrics)

    ## Iterate through terms queue till simplified
    #result = 0
    #while not terms.empty():
    #    # Pick first term
    #    term = terms.get(block=False)
    #    (const, metrics, momentas, matrix_terms) = amp._split_into_matrices(term)

    #    # Continue until term is simplified
    #    while len(matrix_terms) > 0:
    #        # Simplify metrics
    #        metrics = simplify_metrics(metrics)

    #        # Simplify UBar(p) U(p)
    #        if len(matrix_terms) >= 2:
    #            [A, B] = matrix_terms[:2]
    #            A_arg, B_arg = A.args[0], B.args[1]
    #            if isinstance(A, UBar) and isinstance(B, U) and A_arg == B_arg:
    #                const *= -p_to_m[A_arg]  # TODO check sign?
    #                matrix_terms = matrix_terms[2:]

    #                # Resolve term
    #                result += const * metrics * momentas * matrix_terms
    #                break

    #        # Find first gamma matrix
    #        A, i = find_first(matrix_terms, \
    #                          lambda expr: isinstance(expr, Gamma))
    #        A_ind = A.args[0]

    #        # Find the term after
    #        B = matrix_terms[i+1]

    #        if isinstance(B, Gamma):
    #            B_ind = B.args[0]

    #            j = Amplitude._find_metric_contraction(A, B, metrics)
    #            if j is not None:
    #                # Contract gammas with metric
    #                const *= 4
    #                metrics = metrics[:j] + metrics[j+1:]
    #                matrix_terms = matrix_terms[:i] + matrix_terms[i+2:]
    #            else:
    #                # New term
    #                new_matrix_terms = matrix_terms[:i] + matrix_terms[i+2:]
    #                new_const = const * (-2)
    #                new_metrics = metrics + [MetricFactory(A_ind, B_ind)]
    #                terms.put(Amplitude._merge_matrices(new_const, new_metrics, momentas, new_matrix_terms))

    #                # Swap gammas
    #                const *= -1
    #                matrix_terms[i:i+2] = [matrix_terms[i+1], matrix_terms[i]]
    #        elif isinstance(B, U):
    #            B_arg = B.args[0]

    #            # Get A's momenta

    #            if A_ind == B_ind:
    #                # TODO get correct mass


    ################################################
    ########            RENDER            ##########
    ################################################

    #amp.latex_add2(latex)
    latex.render()

if __name__ == "__main__":
    fermion_propagator()
