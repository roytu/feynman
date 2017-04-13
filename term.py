
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

def MomentumFactory(name, ind):
    if ind:
        s = sy.Symbol.__new__(Momentum, "{0}_{1}".format(name, ind))
    else:
        s = sy.Symbol.__new__(Momentum, "{0}".format(name))
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
    pass

class Momentum(sy.Symbol):
    pass

class U(sy.Symbol, MatrixTerm):
    pass

class UBar(sy.Symbol, MatrixTerm):
    pass

class Metric(sy.Symbol):
    pass

# doit dummy classes
class Gamma0(sy.MatrixSymbol):
    def __new__(cls):
        s = sy.MatrixSymbol.__new__(cls, "{{ \\gamma^0 }}", 4, 4)
        #s._assumptions["commutative"] = False
        return s

    @classmethod
    def explicit(self):
        return sy.Matrix([[1, 0,  0,  0],
                          [0, 1,  0,  0],
                          [0, 0, -1,  0],
                          [0, 0,  0, -1]])

class Gamma1(sy.MatrixSymbol):
    def __new__(cls):
        s = sy.MatrixSymbol.__new__(cls, "{{ \\gamma^1 }}", 4, 4)
        #s._assumptions["commutative"] = False
        return s

    @classmethod
    def explicit(self):
        return sy.Matrix([[0,  0, 0, 1],
                          [0,  0, 1, 0],
                          [0, -1, 0, 0],
                          [-1, 0, 0, 0]])

class Gamma2(sy.MatrixSymbol):
    def __new__(cls):
        s = sy.MatrixSymbol.__new__(cls, "{{ \\gamma^2 }}", 4, 4)
        #s._assumptions["commutative"] = False
        return s

    @classmethod
    def explicit(self):
        i = sy.I
        return sy.Matrix([[0,  0,  0, -i],
                          [0,  0,  i,  0],
                          [0,  i,  0,  0],
                          [-i, 0,  0,  0]])

class Gamma3(sy.MatrixSymbol):
    def __new__(cls):
        s = sy.MatrixSymbol.__new__(cls, "{{ \\gamma^3 }}", 4, 4)
        #s._assumptions["commutative"] = False
        return s

    @classmethod
    def explicit(self):
        return sy.Matrix([[0,  0,  1,  0],
                          [0,  0,  0, -1],
                          [1,  0,  0,  0],
                          [0, -1,  0,  0]])
