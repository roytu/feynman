class Term(object):
    """ Term serves as an interface and generic Term class with
        default values for all functions.
    """
    def __init__(self, latex):
        self.__latex = latex

    def is_fraction(self):
        return False

    def latex(self):
        return self.__latex

    def __str__(self):
        return self.latex()

    def copy(self):
        return Term(self.__latex)

    def multiply(self, multiplicand):
        return Prod([self.copy(), multiplicand.copy()])

    def distribute(self):
        return Sum([self.copy()])

    def flatten(self):
        return Sum([self.copy()])

class I(Term):
    def __init__(self):
        pass

    def latex(self):
        return "i"

    def copy(self):
        return I()

class C(Term):
    def __init__(self, latex):
        Term.__init__(self, latex)

    def copy(self):
        return C(self.latex)

class Frac(Term):
    def __init__(self, numer, denom):
        self.numer = numer
        self.denom = denom

    def is_fraction(self):
        return True

    def latex(self):
        numer_latex = self.numer.latex()
        denom_latex = self.denom.latex()
        return "\\frac{{ {0} }}{{ {1} }}".format(numer_latex, denom_latex)

    def copy(self):
        return Frac(self.numer, self.denom)

class Sub(Term):
    """ Sub is a term with a four-vector subscript.
    """
    def __init__(self, latex, index):
        self.__latex = latex
        self.__index = index

    def latex(self):
        return self.__latex + "_{{ {0} }}".format(self.__index)

    def copy(self):
        return Sub(self.__latex, self.__index)

class Super(Term):
    """ Super is a term with a four-vector superscript.
    """
    def __init__(self, latex, index):
        self.__latex = latex
        self.__index = index

    def latex(self):
        return self.__latex + "^{{ {0} }}".format(self.__index)

    def copy(self):
        return Super(self.__latex, self.__index)

class Metric(Term):
    """ Metric is a Lorentz tensor with two indices.
    """
    def __init__(self, index1, index2):
        self.__index1 = index1
        self.__index2 = index2

    def latex(self):
        return "g_{{ {0} {1} }}".format(self.__index1, self.__index2)

    def copy(self):
        return Metric(self.__index1, self.__index2)

class Gamma(Term):
    """ Gamma is a gamma matrix with a four-vector superscript.
    """
    def __init__(self, index):
        self.__index = index

    def latex(self):
        return "\\gamma^{{ {0} }}".format(self.__index)

    def copy(self):
        return Gamma(self.__index)

class Sum(Term):
    def __init__(self, terms):
        self.terms = terms

    def multiply(self, multiplicand):
        """ Multiply each term in the sum by the new term.
            (a + bc) * d => ad + bcd
        """
        return Sum([t.multiply(multiplicand) for t in self.terms])

    def latex(self):
        return "\\left(" + " + ".join([t.latex() for t in self.terms]) + "\\right)"

    def __str__(self):
        return "(" + " + ".join([str(t) for t in self.terms]) + ")"

    def distribute(self):
        return Sum([t.distribute() for t in self.terms])

    def copy(self):
        return Sum([t.copy() for t in self.terms])

    def flatten(self):
        lst = []
        for t in self.terms:
            lst += t.flatten().terms
        return Sum(lst)

class Prod(Term):
    def __init__(self, terms=None):
        if terms:
            self.terms = terms
        else:
            self.terms = []

    def multiply(self, multiplicand):
        return Prod(self.terms + [multiplicand])

    def distribute(self):
        """ Return a Sum object that acts as a distributed Prod. e.g.

            a(b + c) => ab + ac
            a(b + c(d + e)) => ab + acd + ace

                Prod([Term("a"), Sum([Term("b") + Term("c")])])
            => Sum([Prod([Term("a"), Term("b")]), Prod([Term("a"), Term("c")])])
        """
        return Prod([t.distribute() for t in self.terms]).flatten()

    def add(self, term):
        """ Add a term to the products list. """
        self.terms.append(term)

    def add_(self, term):
        """ Add a term to the products list, creating a new Prod object. """
        p = Prod(self.terms[:])
        p.add(term)
        return p

    def latex(self):
        return "\\left(" + " ".join([t.latex() for t in self.terms]) + "\\right)"

    def __str__(self):
        return "(" + " ".join([str(t) for t in self.terms]) + ")"

    def copy(self):
        return Prod([t.copy() for t in self.terms])

    def flatten(self):
        return Sum([Prod([t.flatten() for t in self.terms])])

class Integral(Term):
    def __init__(self, param, norm=False, limits=None):
        self.param = param
        self.norm = norm
        self.limits = limits

    def latex(self):
        if self.limits:
            a, b = self.limits
            limits = "\\limits_{{ {0} }}^{{ {1} }}".format(a.latex(), b.latex())
        else:
            limits = ""

        if self.norm:
            return "\\int {1} \\frac{{ d^4 {0} }}{{ (2 \pi)^4 }}".format(self.param, limits)
        else:
            return "\\int {1} d^4 {0}".format(self.param, limits)

if __name__ == "__main__":
    # multiply()
    print("multiply()")
    print("=" * 30)

    s = Sum([Term("a"), Term("b")])
    print("{0} c => {1}".format(s, s.multiply(Term("c"))))

    s = Sum([Term("a"), Prod([Term("b"), Term("c")])])
    print("{0} d => {1}".format(s, s.multiply(Term("d"))))

    # flatten()
    print("\n")
    print("flatten()")
    print("=" * 30)

    s = Sum([Term("a")])
    print("{0} => {1}".format(s, s.flatten()))

    s = Sum([Term("a"), Term("b")])
    print("{0} => {1}".format(s, s.flatten()))

    s = Sum([Term("a"), Sum([Term("b"), Term("c")])])
    print("{0} => {1}".format(s, s.flatten()))

    s = Prod([Sum([Term("a"), Term("b")]), Term("c")])
    print("{0} => {1}".format(s, s.flatten()))

    s = Sum([Term("a"), Prod([Sum([Term("b"), Term("c")]), Term("d")])])
    print("{0} => {1}".format(s, s.flatten()))

    # distribute()
    print("\n")
    print("distribute()")
    print("=" * 30)

    p = Sum([Term("a"), Term("b"), Term("c")])
    print("{0} => {1}".format(p, p.distribute()))

    p = Prod([Term("a"), Sum([Term("b"), Term("c")])])
    print("{0} => {1}".format(p, p.distribute()))

    p = Prod([Sum([Term("a"), Term("b")]), Term("c")])
    print("{0} => {1}".format(p, p.distribute()))

    ## TODO pass test
    #p = Prod([Term("a"), Sum([Term("b"), Prod([Term(["c", "d"])])])])
    #print("{0} => {1}".format(p, p.distribute()))
