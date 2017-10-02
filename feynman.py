
""" Given a Feynman diagram, perform integrals. """

import re

from math import factorial
import sympy as sy
import numpy as np
from scipy.special import gamma
from latex import Latex
from term import GammaFactory, UFactory, UBarFactory, EFactory, EBarFactory, \
        Momentum, MetricFactory, Momentum, Metric, MatrixTerm

from itertools import product

from util import find_first

#RENDER_ALL = True
RENDER_ALL = False

EPS = 1e-10

def get_highest_log_term(expr, uv):
    """ Return the highest order term of uv of expr

        Args:
            expr: sympy expression
            uv: ultraviolet cutoff
    """
    # If UV is not here, return expr
    if not expr.has(uv):
        return expr

    if expr.func is sy.Add:
        # Only arguments with a value can live
        args_ = [get_highest_log_term(arg, uv) for arg in expr.args if arg.has(uv)]
        return sy.Add(*args_)
    if expr.func in (sy.Mul, sy.log, sy.Pow):
        # Call recursively
        return expr.func(*[get_highest_log_term(arg, uv) for arg in expr.args])
    else:
        # TODO finish
        #import pdb; pdb.set_trace()
        return expr

#def rewrite_out_sqrt(expr):
#    """ Return expr, where all forms x^(a/2) are rewritten as sqrt(x^a)
#    """
#    if not expr.has(sy.Pow)
#        return expr
#
#    if expr.func is sy.Pow:
#        expr.

class TranslatorException(Exception):
    pass

class Translator(object):
    """ A Translator translates symbols into each other.
    
        For instance, it provides the mapping between a four-momenta and its
        scalar symbol:

            k |-> k^2

        and a four-momenta with an index into an indexed four-momenta:

            k, \mu |-> k^\mu
    """
    def __init__(self):
        self.four_momenta = []

    def add_four_momenta(self, k, k2, k_ind_func):
        """ Adds a four-momenta to the translator.

            Args:
                k: four-momenta symbol
                k2: scalar symbol
                k_ind_func: function (k, ind -> k^ind)
        """
        self.four_momenta.append((k, k2, k_ind_func))

    def k_to_k2(self, k):
        for (k_, k2_, k_ind_func) in self.four_momenta:
            if k == k_:
                return k2_
        raise TranslatorException

    def k2_to_k(self, k2):
        for (k_, k2_, k_ind_func) in self.four_momenta:
            if k2 == k2_:
                return k_
        raise TranslatorException

    def k_to_k_ind(self, k, ind):
        for (k_, k2_, k_ind_func) in self.four_momenta:
            if k == k_:
                return k_ind_func(ind)
        raise TranslatorException

class Amplitude(object):
    """ Amplitude is a data structure for holding M expressions. """
    def __init__(self):
        self.const = 1
        self.numer = 1
        self.denom = 1
        self.metrics = []

        # Integrals
        self.integrals_zs = []  # Format is (z, a, b)
        self.integrals_internal = []  # Internal momenta integrals
        self.integrals_cutoffs = []  # Cutoff integral parameters (dt)
        self.spinors = []

        # TODO sloppy
        self.qs = []

        # ??
        self.indices = set()

        self.denom_z = 1
        self.inners = []  # Format is (constant, expr)

    def copy(self):
        """ Returns a deep copy of itself """
        amp = Amplitude()

        amp.const = self.const
        amp.numer = self.numer
        amp.denom = self.denom
        amp.metrics = self.metrics[:]
        amp.integrals_zs = self.integrals_zs[:]
        amp.integrals_internal = self.integrals_internal[:]
        amp.integrals_cutoffs = self.integrals_cutoffs[:]
        amp.spinors = self.spinors[:]
        amp.indices = self.indices.copy()

        # ??
        amp.denom_z = self.denom_z
        amp.inners = self.inners[:]

        return amp

    def U(self, p_name):
        u = UFactory(p_name)
        self.spinors.append(u)
        self.numer *= u
        return u

    def UBar(self, p_name):
        u = UBarFactory(p_name)
        self.spinors.append(u)
        self.numer *= u
        return u

    def E(self, p_name, ind):
        e_ = EFactory(p_name, ind)
        self.spinors.append(e_)
        self.numer *= e_
        return e_

    def EBar(self, p_name, ind):
        e_ = EBarFactory(p_name, ind)
        self.spinors.append(e_)
        self.numer *= e_
        return e_

    def V(self, e, ind):
        """ e: string
            ind: string
        """
        self.const *= sy.I
        self.const *= e
        self.numer *= GammaFactory(ind)
        self.indices.add(ind)

    def S_F(self, p_up, p_down, p_up_dummy, p_down_dummy, m, ind):
        """ Fermionic propagator: i (\slash{p} + m} / (p^2 - m^2)
        
            p_up: contravariant momentum
            p_down: contravariant momentum
            m: symbol
            ind: string
        """
        self.const *= sy.I
        self.numer *= p_down * GammaFactory(ind) + m
        self.denom *= p_down_dummy * p_up_dummy - m ** 2
        self.indices.add(ind)

    def D_F(self, k_up_dummy, k_down_dummy, ind1, ind2, t, lamb, Lamb):
        """ Regulated fermionic propagator
        
            k_up: contravariant dummy symbol
            k_down: covariant dummy symbol
            ind1: string
            ind2: string
            t: symbol
            Lamb: symbol
            lamb: symbol
        """
        self.const *= -sy.I
        self.metrics.append((ind1, ind2))

        # Regulator
        self.const *= -1
        self.integrals_cutoffs.append(
                (t, lamb ** 2, Lamb ** 2)
            )
        self.denom *= (k_down_dummy * k_up_dummy - t) ** 2

        self.indices.add(ind1)
        self.indices.add(ind2)

    def get_latex(self, latex):
        """ Numer / denom format """
        s = ""
        s += latex.get(self.const)

        metrics = 1
        for (a, b) in self.metrics:
            metrics *= MetricFactory(a, b)
        s += latex.get(metrics)

        for (z, a, b) in self.integrals_zs + self.integrals_cutoffs:
            if a is not None and b is not None:
                s += "\\int\\limits_{{ {1} }}^{{ {2} }} d{0}".format(z, latex.get(a), latex.get(b))
            else:
                s += "\\int d{0}".format(z)

        for (k, a, b) in self.integrals_internal:
            if a is not None and b is not None:
                s += "\\int\\limits_{{ {1} }}^{{ {2} }} \\frac{{d^4 {0} }}{{ (2\pi)^4 }}".format(k, latex.get(a), latex.get(b))
            else:
                s += "\\int\\frac{{d^d {0} }}{{ (2\pi)^4 }}".format(k)

        s += "\\left("
        s += latex.get(self.numer / self.denom)
        s += "\\right)"
        return s

    def latex_add(self, latex):
        latex.add(self.get_latex(latex))

    def latex_add2(self, latex):
        """ Inner format with no integrals """
        expr_ = self.expr * sum([c_ * e_ for (c_, e_) in self.inners])
        for (a, b) in self.metrics:
            expr_ *= MetricFactory(a, b)

        s = ""
        s += latex.get(expr_)
        latex.add(s)

    @staticmethod
    def dummy_momentum_up(k):
        """ Make a contravariant momentum with a fake index.
            
            Args:
                k: symbol
        """
        return Momentum(k, "SLASHeta", 1)

    @staticmethod
    def dummy_momentum_down(k):
        """ Make a covariant momentum with a fake index.
            
            Args:
                k: symbol
        """
        return Momentum(k, "SLASHeta", 1)

    @staticmethod
    def square(k):
        """ Make a covariant momentum with a fake index.
            
            Args:
                k: symbol
        """
        return Amplitude.dummy_momentum_down(k) * Amplitude.dummy_momentum_up(k)

    @staticmethod
    def flip_variant(expr):
        """ Change every momentum in expr from contravariant to
            covariant or vice versa and returns the new expression.
        """
        a = sy.Wild("a")
        b = sy.Wild("b")
        c = sy.Wild("c")
        return expr.replace(Momentum(a, b, c), Momentum(a, b, 1 - c))

    @staticmethod
    def replace_momentum(expr, k_up, new_expr_up, wild_ind):
        """ Change every momentum in expr from k1 to k2, preserving
            indices.
        """
        any_k_up = Momentum(k_up.args[0], wild_ind, 1)
        any_k_down = Momentum(k_up.args[0], wild_ind, 0)

        any_expr_up = new_expr_up
        any_expr_down = Amplitude.flip_variant(new_expr_up)

        return expr.replace(any_k_up, any_expr_up) \
                   .replace(any_k_down, any_expr_down)

def calculate(config_str, internal_momenta):
    """ Calculate expression for configuration. """

    latex = Latex()

    ################################################
    ########     CONSTRUCT AMPLITUDE      ##########
    ################################################

    config, amp = make_amplitude(config_str, internal_momenta)

    # Render
    latex.add_text("\\section*{Raw amplitude}")
    amp.latex_add(latex)
    if RENDER_ALL:
        latex.render()

    ################################################
    ########   SIMPLIFY NUMERATOR         ##########
    ################################################

    amp.numer = amp.numer.expand()

    # Render
    latex.add_text("\\section*{Simplified numerator}")
    amp.latex_add(latex)
    if RENDER_ALL:
        latex.render()

    ################################################
    ########      FEYNMAN'S TRICK         ##########
    ################################################

    denom_ = []
    for arg in amp.denom.args:
        if type(arg) == sy.Pow:
            base, power = arg.args
            for _ in range(power):
                denom_.append(base)
        else:
            denom_.append(arg)

    n = len(denom_)
    amp.const *= gamma(n)

    zs = [sy.Symbol("{{ z_{{ {0} }} }}".format(i+1)) for i in range(n)]
    amp.denom = sum([d * z for (d, z) in zip(denom_, zs)]).expand() ** n
    for i, z in enumerate(zs):
        a = 0
        b = 1 - sum(zs[:i])
        amp.integrals_zs.append((z, a, b))

    # Render
    latex.add_text("\\section*{Feynman parameterization}")
    latex.add_text("Here, we perform the following expansion:")
    latex.add_text("""$$
    \\frac{1}{A_1} \\cdots \\frac{1}{A_n} = (n-1)! \\int\\limits_0^1 dz_1
                                                   \\int\\limits_0^{1-z_1} dz_2
                                                   \\cdots
                                                   \\int\\limits_0^{1-z_1-\\cdots-z_{n-1}} dz_n
                                                   \\frac{1}{(z_1 A_1 + \\cdots + z_n A_n)^n}
    $$""")
    latex.add_text("We use this form because a single denominator raised to a power can be simplified with the Golden Integral.")
    amp.latex_add(latex)
    if RENDER_ALL:
        latex.render()

    ################################################
    ######## SPLIT NUMERATOR INTO TERMS   ##########
    ################################################
    amps = []
    for numer_ in sy.Add.make_args(amp.numer):
        amp_ = amp.copy()
        amp_.numer = numer_
        amps.append(amp_)

    # Render
    latex.add_text("\\section*{Expanded numerator}")
    latex.add_text("We split the numerator into additive terms, to process individually.  The following is a list of such terms:")
    for amp_ in amps:
        amp_.latex_add(latex)
    if RENDER_ALL:
        latex.render()

    ################################################
    ########     EVAL. INTERNAL MOMENTA   ##########
    ################################################

    # TODO evaluate internal momenta integrals
    # At this point we stop with numer and denom and combine them into
    # one expression, `inner`, which is a sum of fractions.
    # TODO update this comment


    # Render an explanation
    latex.add_text("\\section*{{Golden Integral}}")
    latex.add_text("We resolve internal momentas with this transformation:")
    latex.add_text("""
    $$\\int \\frac{d^d q}{(2 \pi)^d} \\frac{(q^2)^a}{(q^2 + D)^b} = i \\frac{\\Gamma (b-a-\\frac{1}{2}d) \\Gamma (a + \\frac{1}{2} d)}{(4 \\pi)^{d/2} \\Gamma(b) \\Gamma(\\frac{1}{2}d)} D^{-(b-a-d/2)}$$
    """)
    latex.add_text("After this section, all internal momenta should disappear.  We will now resolve each term in a queue.  Each term may produce additional terms, which are pushed to the back of the queue and resolved later.")

    integrated_amps = []

    #for i, amp_ in enumerate(amps):
    i = 0
    while len(amps) > 0:
        # Pop off one amplitude
        amp_ = amps[0]
        amps = amps[1:]
        i += 1

        latex.add_text("\\section*{{Evaluating internal momenta in this term ({0} terms left)}}".format(len(amps)))
        amp_.latex_add(latex)
        if RENDER_ALL:
            latex.render()

        # Find an internal momenta
        if len(amp_.integrals_internal) > 0:
            (k, _, _) = amp_.integrals_internal[0]
            latex.add_text("Integrating over ${0}$\\\\".format(k))

            # TODO cleanup weird namespacing
            k_down_dummy = Momentum(k, "DUMMY", 0)
            k_up_dummy = Momentum(k, "DUMMY", 1)
            k2_dummy = k_down_dummy * k_up_dummy

            # Decompose denominator
            # denom = denom_nopow ^ b
            denom_nopow, b = amp_.denom.args[0], amp_.denom.args[1]

            # Completing the square
            # Denominator is always quadratic in momenta
            G = denom_nopow  # aliasing for convenience

            latex.add_text("Completing the square\\")

            A = G.collect(k2_dummy).coeff(k2_dummy)
            G = sy.simplify(G - A * k2_dummy)
            B_up = G.collect(k_down_dummy).coeff(k_down_dummy)
            B_down = G.collect(k_up_dummy).coeff(k_up_dummy)
            B = B_up + Amplitude.flip_variant(B_down)
            G = sy.simplify(G - B_up * k_down_dummy - B_down * k_up_dummy)
            C = G
            D = - (B_up * Amplitude.flip_variant(B_up)) / (4 * A) + C

            latex.add("A = " + latex.get(A))
            latex.add("B = " + latex.get(B))
            latex.add("C = " + latex.get(C))

            """ 
                The denominator is in the form:

                    A k^2 + Bk + C

                We define a new variable, q, such that

                    q = A^(1/2) k + B / (2 A^(1/2))

                and replace k:

                    k = q / A^(1/2) - B / (2A)
                    d^d k = (1 / A^(1/2)) d^d q

                The substitution k -> q yields:

                    A k^2 + Bk + C |-> q^2 + D

                where we define D = C - B^2 / (4A)
            """

            # Prepare to replace numerator
            k_name = k
            q_name = "q_{0}".format(len(amp_.qs) + 1)  # TODO sloppy af
            amp_.qs.append(q_name)

            q_up = Momentum(q_name, "DUMMY", 1)
            q_down = Momentum(q_name, "DUMMY", 0)

            any_name = sy.Wild("a")
            any_ind = sy.Wild("b")
            any_variant = sy.Wild("c")

            any_B_up = B_up.replace(Momentum(any_name, "DUMMY", any_variant),
                                    Momentum(any_name, any_ind, any_variant))
            any_B_down = Amplitude.flip_variant(any_B_up)

            # Actually replace numerator
            amp_.numer = amp_.numer.replace(
                Momentum(k_name, any_ind, 1),
                Momentum(q_name, any_ind, 1) / (A ** 0.5) - any_B_up / (2 * A)
            )

            amp_.numer = amp_.numer.replace(
                Momentum(k_name, any_ind, 0),
                Momentum(q_name, any_ind, 0) / (A ** 0.5) - any_B_down / (2 * A)
            )

            amp_.numer = sy.simplify(amp_.numer)

            # Replace denominator
            amp_.denom = (q_down * q_up + D) ** b

            # Replace integral
            # TODO replace integral
            amp_.integrals_internal[0] = (q_name, _, _)
            amp_.numer /= A ** 0.5

            latex.add_text("After ${0} \\to {1}$ substitutions".format(k_name, q_name))
            amp_.latex_add(latex)
            if RENDER_ALL:
                latex.render()

            # Expand the numerator into different amplitudes and multiply them
            # back in the end
            amps__ = []
            for numer in sy.Add.make_args(amp_.numer.expand()):
                new_amp = amp_.copy()
                new_amp.numer = numer
                amps__.append(new_amp)

            # Render
            latex.add_text("\\section*{{Expanding numerator into {0} term(s)}}".format(len(amps__)))
            for amp__ in amps__:
                amp__.latex_add(latex)
            if RENDER_ALL:
                latex.render()

            # Finish q_name integration for each amplitude separately
            for amp__ in amps__:
                prod = sy.Mul.make_args(amp__.numer)

                # Get k-vectors
                # TODO Figure out a way to collect qs nicely
                qs = [q for q in prod if isinstance(q, Momentum) and q.args[0].name == q_name]

                latex.add_text("\\subsection*{Integrating this term:}")
                amp__.latex_add(latex)
                latex.add_text("Found {0} q-vector terms in the numerator.\\\\".format(len(qs)))

                # Simplify q vectors

                # Ward identity for odd tensors
                if len(qs) % 2 == 1:
                    # TODO integral evaluates to zero
                    latex.add_text("Term vanishes due to Ward identity\\\\")
                    amp__.const = 0
                    continue

                if len(qs) > 0:
                    # TODO convert higher-order even tensor integral to
                    #      scalar integral
                    pass
                else:
                    # TODO Assume a = 0 for now
                    # This is obviously wrong in general but will be easier
                    # to fix with a good test case
                    a = 0

                    # Golden integral
                    #c, a = term.as_coeff_exponent(sy.Symbol(k2))
                    c_ = sy.I * gamma(np.float64(b - a - (4 - EPS)/2)) * gamma(np.float64(a + (4 + EPS) / 2))
                    c_ /= gamma(np.float64(b))
                    c_ /= (4 * sy.pi) ** 2
                    # Part of the Golden integral d^q factor
                    c_ *= (2 * sy.pi) ** 4
                    amp__.const *= c_

                    amp__.denom = D ** (b - a - 2)  # TODO generalize to d-dimensions
                                                    # with 2 -> D / 2

                    # Add to amps if nonzero
                    amps.append(amp__)

                    # Remove internal integral
                    amp__.integrals_internal = amp__.integrals_internal[1:]

                    # Render
                    latex.add_text("Apply golden integral")
                    amp__.latex_add(latex)
                if RENDER_ALL:
                    latex.render()
        else:
            integrated_amps.append(amp_)

        ## Compress inners by term
        #inners_dict = {}
        #for (c_, expr_) in amp.inners:
        #    if expr_ not in inners_dict:
        #        inners_dict[expr_] = 0
        #    inners_dict[expr_] += c_
        #amp.inners = [(v, k) for (k, v) in inners_dict.items()]
        #amp.latex_add2(latex)

    amps = integrated_amps

    latex.add_text("\\section*{Final amplitudes after momenta integration}")
    for amp_ in amps:
        amp_.latex_add(latex)

    if RENDER_ALL:
        latex.render()

    ##########################################
    ########  CUTOFF INTEGRATIONS   ##########
    ##########################################
    latex.add_text("\\section*{Integrating cutoffs}")
    latex.add_text("Here we integrate all $t$-variables, which represent the upper and lower cutoffs.")
    uv = sy.Symbol(config["Lamb"])

    integrated_amps = []

    for amp_ in amps:
        latex.add_text("\\subsection*{Integrating this term}")
        amp_.latex_add(latex)

        expr_ = 1 / amp_.denom

        latex.add_text("Denominator only")
        latex.add(latex.get(expr_))
        if RENDER_ALL:
            latex.render()

        # Integrate cutoffs
        for (t, a, b) in amp_.integrals_cutoffs:
            # Integrate w.r.t. cutoff
            expr_ = sy.integrate(expr_, (t, a, b))

            latex.add_text("Integrating wrt ${0}$...".format(t))
            latex.add(latex.get(expr_))
            if RENDER_ALL:
                latex.render()

            # Collecting highest order term
            old_expr = expr_
            while True:
                old_expr = expr_

                expr_ = sy.expand_log(expr_, force=True)
                expr_ = get_highest_log_term(expr_, uv)
                expr_ = sy.simplify(expr_)

                if old_expr == expr_:
                    # No more changes
                    break

            latex.add_text("Keeping only highest order term...")
            latex.add(latex.get(expr_))
            if RENDER_ALL:
                latex.render()

            amp_.numer *= expr_
            amp_.denom = 1

        amp_.integrals_cutoffs = []
        integrated_amps.append(amp_)

    amps = integrated_amps

    ######################################
    ########   Z INTEGRATIONS   ##########
    ######################################
    latex.add_text("\\section*{Integrating $z$-variables}")
    latex.add_text("Here we integrate all $z$-variables, the Feynman parameters.")

    integrated_amps = []

    for amp_ in amps:
        latex.add_text("\\subsection*{Integrating this term}")
        amp_.latex_add(latex)

        # Integrate cutoffs
        expr_ = amp_.numer

        for (z, a, b) in amp_.integrals_zs[::-1]:
            # Integrate w.r.t. cutoff
            # Rationalize decimal powers first, or sympy breaks
            expr_ = sy.nsimplify(expr_, tolerance=0.001, rational=True)
            expr_ = sy.integrate(expr_, (z, a, b))

            latex.add_text("Integrating wrt ${0}$...".format(z))
            latex.add(latex.get(expr_))
            if RENDER_ALL:
                latex.render()

        amp_.numer = expr_
        amp_.integrals_zs = []

        integrated_amps.append(amp_)

    amps = integrated_amps

    latex.add_text("\\section*{Final amplitudes after $z$-variable integration}")
    for amp_ in amps:
        amp_.latex_add(latex)

    if RENDER_ALL:
        latex.render()

    ################################################
    ########  EVALUATE SPINS AND GAMMAS   ##########
    ################################################

    latex.add_text("\\section*{Evaluating spins and gamma matrices}")
    latex.add_text("TODO jk do it yourself you slags, here's the sum, have fun.  Don't forget to take traces/multiply by -1 for internal fermion loops.")

    latex.add("\\; + \\;".join([amp_.get_latex(latex) for amp_ in amps]))

    if RENDER_ALL:
        latex.render()

    ################################################
    ########            RENDER            ##########
    ################################################

    #amp.latex_add2(latex)
    latex.render()

def make_amplitude(config_str, internal_momenta):
    """ Make amplitude from config string and internal momenta list. """
    # Replace slashes
    config_str = config_str.replace("\\", "SLASH")

    # Parse terms
    raw_terms = config_str.split()
    terms = []
    for raw_term in raw_terms:
        i = raw_term.find("(")
        j = raw_term.find(")")
        func = raw_term[:i]
        args = raw_term[i+1:j].split(",")
        terms.append((func, args))

    # Start building config
    config = {}
    config["e"] = "e"
    config["m"] = "m"
    config["Lamb"] = "SLASHLambda"
    config["lamb"] = "SLASHlambda"
    config["t"] = "t"

    # Internal momenta
    config["internal_momenta"] = internal_momenta
    config["momenta"] = {}

    ##############################
    ###    Parse amplitude     ###
    ##############################
    amp = Amplitude()
    amp.const /= (2 * sy.pi) ** 4

    # Add internal integrals for all internal momenta
    for k in config["internal_momenta"]:
        amp.integrals_internal.append((k, None, None))

    # Go through terms
    for (func, args) in terms:
        if func == "UBar":
            p_name = args[0]
            amp.UBar(p_name)
        elif func == "U":
            p_name = args[0]
            amp.U(p_name)
        elif func == "EBar":
            p_name = args[0]
            ind = args[1]
            amp.EBar(p_name, ind)
        elif func == "E":
            p_name = args[0]
            ind = args[1]
            amp.E(p_name, ind)
        elif func == "V":
            mu_name = args[0]
            amp.V(sy.Symbol(config["e"]), mu_name)
        elif func == "S":
            [ps, ind] = args

            # Tokenize
            p_names = re.split("\+|-", ps)
            ops = [c for c in ps if c in "+-"]

            # Make new momentums if they don't exist
            for p_name in p_names:
                if p_name not in config["momenta"]:
                    p_ = config["momenta"][p_name] = {}
                    p_["up"] = Momentum(p_name, ind, 1)
                    p_["down"] = Momentum(p_name, ind, 0)
                    p_["up_dummy"] = Momentum(p_name, "DUMMY", 1)
                    p_["down_dummy"] = Momentum(p_name, "DUMMY", 0)

            p_ = config["momenta"][p_names[0]]
            up_sum = p_["up"]
            down_sum = p_["down"]
            up_dummy_sum = p_["up_dummy"]
            down_dummy_sum = p_["down_dummy"]

            for p_name, op in zip(p_names[1:], ops):
                p_ = config["momenta"][p_name]
                if op == "+":
                    up_sum += p_["up"]
                    down_sum += p_["down"]
                    up_dummy_sum += p_["up_dummy"]
                    down_dummy_sum += p_["down_dummy"]
                elif op == "-":
                    up_sum -= p_["up"]
                    down_sum -= p_["down"]
                    up_dummy_sum -= p_["up_dummy"]
                    down_dummy_sum -= p_["down_dummy"]

            amp.S_F(up_sum, down_sum, up_dummy_sum, down_dummy_sum,
                    sy.Symbol(config["m"]), ind)
        elif func == "D":
            [ps, mu, nu] = args

            # Tokenize
            p_names = re.split("\+|-", ps)
            ops = [c for c in ps if c in "+-"]

            # Make new momentums if they don't exist
            for p_name in p_names:
                if p_name not in config["momenta"]:
                    p_ = config["momenta"][p_name] = {}
                    p_["up"] = Momentum(p_name, ind, 1)
                    p_["down"] = Momentum(p_name, ind, 0)
                    p_["up_dummy"] = Momentum(p_name, "DUMMY", 1)
                    p_["down_dummy"] = Momentum(p_name, "DUMMY", 0)

            p_ = config["momenta"][p_names[0]]
            up_dummy_sum = p_["up_dummy"]
            down_dummy_sum = p_["down_dummy"]

            for p_name, op in zip(p_names[1:], ops):
                p_ = config["momenta"][p_name]
                if op == "+":
                    up_dummy_sum += p_["up_dummy"]
                    down_dummy_sum += p_["down_dummy"]
                elif op == "-":
                    up_dummy_sum -= p_["up_dummy"]
                    down_dummy_sum -= p_["down_dummy"]

            amp.D_F(up_dummy_sum, down_dummy_sum, mu, nu,
                    sy.Symbol(config["t"]),
                    sy.Symbol(config["lamb"]),
                    sy.Symbol(config["Lamb"]))

    return (config, amp)

if __name__ == "__main__":
    # Electron self-energy correction
    #config_str = """UBar(p) V(\\mu) S(p-k,\\sigma_2) V(\\nu) U(p) D(k,\\mu,\\nu)"""
    #internal_momenta = ["k"]
    #calculate(config_str, internal_momenta)

    # Photon self-energy correction
    config_str = """EBar(k,\\mu) V(\\mu) S(k+p,\\sigma_1) S(p,\\sigma_2) V(\\nu) E(k,\\nu)"""
    internal_momenta = ["p"]
    calculate(config_str, internal_momenta)
