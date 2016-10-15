
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

class Gamma(sy.MatrixSymbol, MatrixTerm):
    def __new__(cls, ind):
        s = sy.Symbol.__new__(cls, "\\gamma_{0}".format(ind))
        s._args += (ind,)
        s._assumptions["commutative"] = False
        return s

    def resolve(self, ind):
        if ind == 0:
            return Gamma0()
        elif ind == 1:
            return Gamma1()
        elif ind == 2:
            return Gamma2()
        elif ind == 3:
            return Gamma3()
        else:
            raise Exception("WEIRD")

    def doit(self):
        ret = None
        if self.args[0] == 0:
            ret = Gamma0.explicit()
        elif self.args[0] == 1:
            ret = Gamma1.explicit()
        elif self.args[0] == 2:
            ret = Gamma2.explicit()
        elif self.args[0] == 3:
            ret = Gamma3.explicit()
        else:
            raise Exception("Uncontracted index! {0}".format(self.args[0]))
        print("Returning {0}".format(ret))
        return ret

class Momentum(sy.Symbol):
    def __new__(cls, name_, ind=None):
        if ind:
            s = sy.Symbol.__new__(cls, "{0}_{1}".format(name_, ind))
        else:
            s = sy.Symbol.__new__(cls, "{0}".format(name_))
        s._args += (name_,)
        if ind:
            s._args += (ind,)
        return s

    def resolve(self, new):
        """ This overrides the default behavior of subs. """
        # TODO make better?
        return sy.Symbol("{{ {{ {0} }}_{{ {1} }} }}".format(self.args[0], new))

class U(sy.MatrixSymbol, MatrixTerm):
    def __new__(cls, name_, spin):
        s = sy.Symbol.__new__(cls, "u({0})".format(name_))
        s._args += (name_, spin)  # Spin is either 0: up or 1: down
        s._assumptions["commutative"] = False
        return s

    def explicit(self):
        p = Momentum(self.args[0])
        E = p.resolve(0)
        p1 = p.resolve(1)
        p2 = p.resolve(2)
        p3 = p.resolve(3)

        m = sy.sqrt(E ** 2 - p1 ** 2 - p2 ** 2 - p3 ** 2)

        N = sy.sqrt((E + m) / (2 * m))
        if self.args[1] == 0:
            return sy.Matrix([[1],
                              [0],
                              [p3 / ( E + m )],
                              [(p1 + sy.I * p2) / (E + m)]])

            return sy.Matrix([[1],
                              [0],
                              [p3 / ( E + m )],
                              [(p1 + sy.I * p2) / (E + m)]])
        elif self.args[1] == 1:
            return sy.Matrix([[0],
                              [1],
                              [(p1 - sy.I * p2) / (E + m)],
                              [-p3 / ( E + m )]])
        else:
            raise Exception("Error!")

class UBar(sy.MatrixSymbol, MatrixTerm):
    def __new__(cls, name_, spin):
        s = sy.Symbol.__new__(cls, "\\bar{{u}}({0})".format(name_))
        s._args += (name_, spin)  # Spin is either 0: up or 1: down
        s._assumptions["commutative"] = False
        return s

    def explicit(self):
        p = Momentum(self.args[0])
        E = p.resolve(0)
        p1 = p.resolve(1)
        p2 = p.resolve(2)
        p3 = p.resolve(3)

        m = sy.sqrt(E ** 2 - p1 ** 2 - p2 ** 2 - p3 ** 2)

        N = sy.sqrt((E + m) / (2 * m))
        if self.args[1] == 0:
            return sy.Matrix([[1,
                               0,
                               p3 / (E + m),
                               (p1 - sy.I * p2) / (E + m)]])
        elif self.args[1] == 1:
            return sy.Matrix([[0,
                               1,
                               (p1 + sy.I * p2) / (E + m),
                               -p3 / (E + m)]])
        else:
            raise Exception("Error!")

class Metric(sy.Symbol):
    def __new__(cls, ind1, ind2):
        s = sy.Symbol.__new__(cls, "g_{0}_{1}".format(ind1, ind2))
        s._args += (ind1, ind2)
        return s

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
