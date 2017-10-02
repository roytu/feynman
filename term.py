
import sympy as sy

""" NAMING CONVENTIONS

    All symbols are named with the following style:
    `name`_`index1`_`index2`_...

    where `name` does not include any underscores.

    This is so functions like `replace` work correctly.  If special printing
    behavior is needed the __str__ function is overridden.
"""

class MatrixTerm(object):
    pass

def GammaFactory(ind):
    s = sy.Symbol.__new__(Gamma, "\\gamma_{0}".format(ind), commutative=False)
    s._args += (ind,)
    return s

#def MomentumFactory(name, ind, contravariant=True):
#    #if contravariant:
#    #    #s = sy.Symbol.__new__(Momentum, "{{ {{ {0} }}^{{ {1} }} }}".format(name, ind), commutative=False)
#    #    s = sy.Function(Momentum, "M")
#    #else:
#    #    #s = sy.Symbol.__new__(Momentum, "{{ {{ {0} }}_{{ {1} }} }}".format(name, ind), commutative=False)
#    #    s = sy.Function(Momentum, "M")
#    ##s = Momentum("{0}".format(name))
#    ##s._args += (name, ind, contravariant)
#    return Momentum(name, ind, contravariant)
#    #return s

def EFactory(name, ind):
    s = sy.Symbol.__new__(E, "\\epsilon_{1}({0})".format(name, ind), commutative=True)
    #s._args += (name, spin)
    s._args += (name, ind)
    return s

def EBarFactory(name, ind):
    s = sy.Symbol.__new__(EBar, "\\bar{{\\epsilon}}_{1}({0})".format(name, ind), commutative=True)
    #s._args += (name, spin)
    s._args += (name, ind)
    return s

def UFactory(name):
    s = sy.Symbol.__new__(U, "u({0})".format(name), commutative=False)
    #s._args += (name, spin)
    s._args += (name,)
    return s

def UBarFactory(name):
    s = sy.Symbol.__new__(UBar, "\\bar{{u}}({0})".format(name), commutative=False)
    #s._args += (name, spin)
    s._args += (name,)
    return s

def MetricFactory(ind1, ind2):
    s = sy.Symbol.__new__(Metric, "g_{0}_{1}".format(ind1, ind2))
    s._args += (ind1, ind2)
    return s

# These are actually factories I guess
class Gamma(sy.Symbol, MatrixTerm):
    @property
    def is_rational_function(self):
        return True  # ??

#class Momentum(sy.Symbol):
#    pass

#class Momentum(sy.Symbol, sy.Function):
class Momentum(sy.Function):
    # Inspired by
    # http://docs.sympy.org/0.7.2/_modules/sympy/core/power.html#Pow

    nargs = (3,)

    @property
    def is_commutative(self):
        return True

    @property
    def is_rational_function(self):
        return True  # ??

class U(sy.Symbol, MatrixTerm):
    @property
    def is_rational_function(self):
        return True  # ??

class UBar(sy.Symbol, MatrixTerm):
    @property
    def is_rational_function(self):
        return True  # ??

class E(sy.Symbol):
    @property
    def is_rational_function(self):
        return True  # ??

class EBar(sy.Symbol):
    @property
    def is_rational_function(self):
        return True  # ??

class Metric(sy.Symbol):
    @property
    def is_rational_function(self):
        return True  # ??
