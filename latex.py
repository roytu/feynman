
""" Constructs a report with LaTeX and compiles it. """

import sympy as sy
from sympy.printing.latex import LatexPrinter
from subprocess import call

import os

class Latex(object):
    def __init__(self):
        self.eqs = []
        self._dummies = 0
    
    def get(self, expr):
        return CustomLatexPrinter().doprint(expr)

    def add(self, s):
        latex = ""
        latex += "\\begin{dmath}"
        latex += s
        latex += "\\end{dmath}"
        latex += "\n"

        self.eqs.append(latex)

    def add_text(self, s):
        latex = ""
        latex += s
        latex += "\n"

        self.eqs.append(latex)

    def render(self, fname="report.tex"):
        # Construct equation string
        s = "".join(self.eqs)

        # Replace special words
        s = s.replace("SLASH", "\\")
        #s = s.replace("DUMMY", ":)")
        s = s.replace("DUMMY", "\\eta")

        with open("reports/include.tex", "w") as f:
            f.write(s)

        #os.system("cd reports; lualatex --halt-on-error > /dev/null {0}".format(fname))
        os.system("cd reports; lualatex {0}".format(fname))

class CustomLatexPrinter(LatexPrinter):
    def _print_Gamma(self, expr):
        return "{{ \\gamma^{{ {0} }} }}".format(expr.args[0])

    def _print_Momentum(self, *args):
        expr = args[0]
        if expr.args[2]:
            return "{{ {{ {0} }}^{{ {1} }} }}".format(expr.args[0], expr.args[1])
        else:
            return "{{ {{ {0} }}_{{ {1} }} }}".format(expr.args[0], expr.args[1])

    def _print_U(self, expr):
        return "u({{ {0} }})".format(expr.args[0])

    def _print_UBar(self, expr):
        return "{{ \\bar{{u}}({0}) }}".format(expr.args[0])

    def _print_E(self, expr):
        return "\\epsilon_{{ {1} }}({{ {0} }})".format(expr.args[0], expr.args[1])

    def _print_EBar(self, expr):
        return "\\bar{{\\epsilon}}_{{ {1} }}({{ {0} }})".format(expr.args[0], expr.args[1])

    def _print_Metric(self, expr):
        return "g_{{ {0} {1} }}".format(expr.args[0], expr.args[1])
