import re
import os
import math
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt


def delt(s):
    return r"\Delta" + str(s) if '\\' in str(s) else r"\Delta\!" + str(s)


def ax_name(tab, qty):
    s = f"${sp.latex(sp.sympify(preSymp(qty)))}$"
    if not str(tab.units[qty]) == '1':
        s += r" ( $\rm {:}$ )".format(tab.units[qty])
    return s


def preSymp(string):
    """
    returns a string ready to be passed to Sympy.

    Sympy sometimes assumes the meaning of certain variables, this fuction
    prevents that.

    ex:
        >>> import sympy as sp
        >>> sp.sympify('N**2')
        ... TypeError: unsupported operand type(s) for ** or pow(): 'function'
         and'Integer'

        Here sympy assumed that N was it's own function sympy.N. In the case
        where you would want to use N as a variable instead, you could use:

        >>> import sympy as sp
        >>> from tablpy.extra_funcs import preSymp
        >>> sp.sympify(preSymp('N**2'))
        ... N**2
    """
    string = string.replace("{", "__0__").replace("}", "__1__")
    pat = (r"(?!ln|log|sin|cos|tan|exp|atan|sqrt)"
           r"\\?[a-zA-Z]+[_0-9a-zA-Z({{)(}}),]*")
    a = [f"Symbol('{i}')" for i in re.findall(pat, string)]
    return re.sub(pat, r"{}", string).format(*a).replace("__0__", "{")\
                                                .replace("__1__", "}")


def roundUp(num):  # arrondis vers le haut à la bonne décimale
    pow = math.floor(sp.log(num, 10).evalf())  # Position de la décimla
    snum = str(num).replace(".", "")
    x = 0
    while snum[x] == "0":  # Trouve le premier chiffre non-nul
        x += 1
    # Revoie arrodis vers le haut
    if "0" * (len(snum) - x - 1) != snum[x + 1:] and x < len(snum):
        return ((int(snum[x]) + 1) * 10**sp.sympify(pow)).evalf(chop=True)
    else:
        return num  # pas de correction à faire!


# Détermine la formule d'incertitude pour une expression donnée
def formule_incertitude(eq):
    r"""
    Returns a sympy expression of the uncertainty associated with a formula

    ex:
        >>> from tablpy.extra_funcs import formule_incertitude
        >>> import sympy as sp

        >>> eq = sp.sympify("(a+b)*c")
        >>> formule_incertitude(eq)
        ... sqrt(\Delta\!a**2*c**2 + \Delta\!b**2*c**2
                 + \Delta\!c**2*(a + b)**2)
    """
    eq = sp.sympify(eq)
    variables = list(eq.free_symbols)  # liste de tout les variables
    # liste de touts les incertitudes asssociées au variables
    uncertain = [sp.Symbol(delt(x)) for x in variables]
    fIncert = sp.sqrt(sum(
        [(sp.diff(eq, variables[i]) * uncertain[i])**2
         for i in range(len(variables))]))
    return sp.simplify(fIncert)


def clearfile(file):
    for i in os.listdir(file):
        if os.path.isdir(file + "/" + i):
            clearfile(file + "/" + i)
            os.rmdir(file + "/" + i)
        else:
            os.remove(file + "/" + i)


def extrem(x):
    """return minimum, maximum of a list"""
    return float(min(x)), float(max(x))


def isUncertain(string):
    marqeurs = ['delta',
                'Delta',
                'Δ',
                'Î”']
    return any([marqeur in str(string) for marqeur in marqeurs])


def getTexIncert(eq):
    print(sp.latex(formule_incertitude(eq)))


def fastPlot(file):
    """plots the two first colums of a data file"""
    array = np.genfromtxt(file, skip_header=True, delimiter="\t")
    plt.plot(array[:, 0], array[:, 1])


def rejectOutliers(data, n=2):
    return data[abs(data - np.mean(data)) < n * np.std(data)]


def roll(values, min, max):
    return ((values-min) % (max-min))+min
