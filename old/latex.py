
""" Constructs a report with LaTeX and compiles it. """

import sympy as sy
from subprocess import call

import os

class Latex(object):
    def __init__(self):
        self.eqs = []
        self._dummies = 0
    
    def add(self, s):
        latex = ""
        latex += "\\begin{dmath}"
        latex += s
        latex += "\\end{dmath}"
        latex += "\n"

        self.eqs.append(latex)

    def render(self, fname="report.tex"):
        # Construct equation string
        s = "".join(self.eqs)
        with open("reports/include.tex", "w") as f:
            f.write(s)

        os.system("cd reports; lualatex --halt-on-error > /dev/null {0}".format(fname))
        #os.system("cd reports; lualatex {0}".format(fname))
