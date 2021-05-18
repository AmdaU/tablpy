# %%

from . import extra_funcs as ef
#import extra_funcs as ef
from .units import unit
#from units import unit
from copy import deepcopy
import inspect
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy.optimize import curve_fit
import sympy as sp

exval, exUnis = {}, {}

showUncertain = True

# %%


class table:
    """\n

    The table class is the core of this package. Each table object contains
    columns of data most often imported from a .xlsx file of .csv like file.

    To initialize as table object, simply give it the name of the file you wish
    to import data from.

    ex:
        **having a file named foo.csv id the same directory**
        >>> import tablypy as tp
        >>> dt = tp.table('foo')
        >>> dt
        ...     bar   baz
        ... 0   0.0  1.00
        ... 1   1.0  0.50
        ... 2   6.0  0.90
        ... 3   7.0  0.60
        ... 4   8.0  0.01
        ... 5   9.0  3.00
        ... 6  12.0  6.30
        ... 7  15.0  8.90

    Parameters
    ----------
    (datname) the name or path to the file you want to import from, the
        extention is not required (only .csv, .txt and .xlsx are currently
        supported)

    (AutoInsert) Tables always excpect there to be an uncertainty column with
        every data column. If this paremeter is on, missing uncertainty column
        will be automatically identified and added (full of 0). If you set this
        parmeter to False, it will be assumed that every second column in your
        data is an uncertainty one

    (units) You can associate units with each of your columns unsint the
        dt.giveUnits() mehtod later or you can give the staight away using
        the units parameter. This sould be a dictionnary with the keys being
        the column name and the value being the unit in the form of a string

    (sheet) This is used to specify wich sheet of a .xlsx is to be used

    (data) if you do not want to import from a data file, you can pass a
        np.array to this argument. The datname you specified will be ignored.
    """

    def __init__(self, datname, AutoInsert=True, units={}, sheet=None,
                 data=None, delimiter=',', skiprows=None, noNames=False):
        self.units = units
        if data is not None:
            if isinstance(data, pd.DataFrame):
                self.data = data
            else:
                self.data = pd.DataFrame(data)
        else:
            paths_parts = datname.split('/')
            if '.' in paths_parts[-1]:
                parts = datname.split('.')
                name = '.'.join(parts[:-1])
                ext = parts[-1]
            else:
                name = datname
                possible_ext = ['xlsx', 'csv', 'txt', 'dat']
                found = False
                i = 0
                while not found and i < len(possible_ext):
                    if os.path.isfile(datname + '.' + possible_ext[i]):
                        found = True
                        ext = possible_ext[i]
                    i += 1
                if not found:
                    print("File '{}' was not found :(".format(datname))
                    return

            if ext == 'xlsx':
                self.data = pd.read_excel(name + '.' + ext, sheet_name=sheet)
            else:
                self.data = pd.read_csv(name + '.' + ext, delimiter=delimiter,
                                        skiprows=skiprows)
                if noNames:
                    names = list(map(str, range(len(self.data.columns))))
                    self.data = pd.read_csv(name + '.' + ext,
                                            delimiter=delimiter,
                                            skiprows=skiprows, header=0,
                                            names=names)

        self.formulas = {}

        if AutoInsert:
            lol = True
            if not ef.isUncertain(*self.data.columns[0:2]):
                self.data.insert(1, ef.delt(
                    self.data.columns[0]), [0 for i in range(len(self.data))])

            while lol:
                lol = False
                for i in range(1, len(self.data.columns) - 1):
                    if not ef.isUncertain(*self.data.columns[i:i+2]) and\
                       not ef.isUncertain(*self.data.columns[i-1:i+1]):
                        self.data.insert(i + 1, ef.delt(self.data.columns[i]),
                                         [0 for i in range(len(self.data))])
                        lol = True

            if not ef.isUncertain(*self.data.columns[-2:]):
                self.data.insert(len(self.data.columns), ef.delt(
                    self.data.columns[-1]), [0 for i in range(len(self.data))])
        self.data.columns = list(map(str, self.data.columns))
        for i in self.data.columns:
            self.units[i] = unit('1')
        self.giveUnits(units)

    def __repr__(self):
        dat = [c for c in list(self.data) if 'Delta' not in c]
        return str(self.data if showUncertain else self.data[dat])

    def __getitem__(self, x):
        return np.array(self.data[x])

    def __len__(self):
        return len(self.data)

    def __setitem__(self, c, v):
        self.data[c] = v

    def newCol(self, name, expression, extra=None, pos=None, NoIncert=False,
               units=None):
        """\n

        This methode allows you to add a new columns to your table given a
        certain formula. The uncertainty on the new values and units will be
        automatically taken care of.

        ex:
        #Context

        - you have a file named foo.csv in the same directory
        - you have set showUncertain = False

        >>> dt = table("foo")
        >>> dt
        ...     bar   baz
        ... 0   0.0  1.00
        ... 1   1.0  0.50
        ... 2   6.0  0.90
        ... 3   7.0  0.60
        ... 4   8.0  0.01
        ... 5   9.0  3.00
        ... 6  12.0  6.30
        ... 7  15.0  8.90

        >>> dt.newCol('qux', 'bar + baz*2')
        >>> dt
        ...     bar   baz    qux
        ... 0   0.0  1.00   2.00
        ... 1   1.0  0.50   2.00
        ... 2   6.0  0.90   7.80
        ... 3   7.0  0.60   8.20
        ... 4   8.0  0.01   8.02
        ... 5   9.0  3.00  15.00
        ... 6  12.0  6.30  24.60
        ... 7  15.0  8.90  32.80
        """
        if pos is None:
            pos = len(self.data.columns)
        expression = sp.sympify(ef.preSymp(expression))
        self.formulas[name] = expression

        out = np.empty((len(self.data), 2))

        pre_dict = dict(zip([f"dummy{x}" for x in range(
            len(expression.free_symbols))], map(str, expression.free_symbols)))
        expression = expression.subs({v: k for k, v in pre_dict.items()})
        lamb = sp.lambdify(tuple(expression.free_symbols), expression)
        V_dict = self.data[:].to_dict("records")
        if not NoIncert:
            expr_incert = ef.formule_incertitude(self.formulas[name])
            pre_dict_i = dict(zip(
                [f"dummy{x}" for x in range(len(expr_incert.free_symbols))],
                map(str, expr_incert.free_symbols)
            ))
            expr_incert = expr_incert.subs(
                {v: k for k, v in pre_dict_i.items()})
            lamb_incert = sp.lambdify(tuple(expr_incert.free_symbols),
                                      expr_incert)
        # Le cacul est fait ligne par ligne, probablement très optimisable
        for i in range(len(self.data)):
            vals = V_dict[i]
            vals.update(exval)
            nvals = {name: vals[pre_dict[name]] for name in pre_dict}
            delt_vals = {name: vals[pre_dict_i[name]] for name in pre_dict_i}
            out[i] = [
                lamb(**nvals), 0 if NoIncert else lamb_incert(**delt_vals)]
        # nameme les colonnes et les ajoute au tableau
        if units:
            self.units[name] = unit(units)
        else:
            self.units[name] = lamb(*(
                self.units[str(i)]
                if str(i) in self.units
                else exUnis[str(i)]for i in self.formulas[name].free_symbols)
            )
        self.data.insert(pos, name, out[:, 0])
        self.data.insert(pos + 1, ef.delt(name), out[:, 1])
        self.fixUnits()

    def delCol(self, names):
        """Delets specified column of columns from the table"""
        if not isinstance(names, list):
            names = [names]
        for name in names:
            self.data = self.data.drop(columns=[name, ef.delt(name)])

    def giveUnits(self, units):
        """\n

        This allows you to associate units to your columns

        Ex:
        >>> dt = table('foo')
        >>> dt.giveUnits({bar: "GHz"})

        This useful because not only will it show up in graph and .tex talbe
        you output, but will also be used to evaluate the units of a column
        created with dt.newCol().

        You can also use in in combination with the dt.changeUnits method to
        make for easy unit conversion. (Only works with SI units at the moment)
        """
        for col in self.data.columns:
            try:
                if isinstance(units[col], str):
                    units[col] = unit(units[col])
                self.units[col] = units[col]
            except:
                try:
                    self.units[col] = self.units[col]
                except:
                    self.units[col] = unit("1")

    def changeUnits(self, units):
        """\n

        Change the units and changing the values of the column accordingly.
        If you want to change te units without changing the values, use the
        .giveUnits method instead.

        Ex.


        """
        for col in units:
            fact = self.units[col].to(units[col])
            self.data[col] *= float(fact)
            self.data[ef.delt(col)] *= float(fact)

    def renameCols(self, names):
        """\n

        renames the columns of your dataTable

        (names) can either be a string with all the name separeted by spaces,
            a dictionary mapping old name to new names or a list of new names.
        """
        # vars = list(sp.symbols(noms))
        if isinstance(names, str):
            names = names.split(" ")
            for i in range(len(names)):
                self.units[names[i]] = self.units[self.data.columns[::2][i]]
            self.data.columns = [n for var in names
                                 for n in (var, ef.delt(var))]

        if isinstance(names, dict):
            for key in list(names.keys()):
                self.units[names[key]] = self.units[key]
            keys = [n for var in names.keys()
                    for n in (var, ef.delt(var))]
            vals = [n for var in names.values()
                    for n in (var, ef.delt(var))]
            self.data = self.data.rename(columns=dict(zip(keys, vals)))

        if isinstance(names, list):
            for i in range(len(names)):
                self.units[names[i]] = self.units[self.data.columns[::2][i]]
            self.data.columns =\
                [n for var in names for n in (var, ef.delt(var))]

    def squish(self):
        """Combines data and uncertainty collumns. Used with makeGoodTable"""
        out = pd.DataFrame([])
        for col in [str(i) for i in self.data if "Delta" not in str(i)]:
            outCol = []
            # prendre toutes les valeurs + leurs delta individuellement
            for x in range(len(self.data[col])):
                # calculer le nombre de nombre de chiffres significatif requis
                val = abs(self.data[col][x])
                d = self.data[ef.delt(col)][x]
                # formater la colonne en fonction des résultats
                if d == 0:
                    s = str(val)
                    while s.endswith('0') or s.endswith('.'):
                        s = s[:-1]
                    outCol.append("$" + s + "$")
                else:
                    miam = "{{:.{}g}}"\
                        .format(- np.ceil(-sp.log(abs(val), 10))
                                + np.ceil(-sp.log(d, 10)) + 1).format(val)
                    outCol.append(
                        "$" + miam + " \\pm " + "{:.1g}$"
                        .format(ef.roundUp(self.data[ef.delt(col)][x]))
                    )
            out["$" + col + "$"] = outCol
        return(out)

    # Export un ficher .tex en faisant toutes les modifications nécessaire
    def makeGoodTable(self, name, style='default'):
        """\n

        Outputs a latex version of itself in directory named tables

        (name) name of the exported file.
        """
        if 'tables' not in os.listdir():
            os.mkdir("tables")

        self.fixUnits()
        exp = self.squish()
        names = []
        for col in self.data.columns[::2]:
            if str(self.units[col].SIval) in ["1", '1.00000000000000']:
                names.append("${}$".format(col))
            else:
                names.append("${}$ {{$\rm({})$}}".format(
                    sp.latex(sp.sympify(ef.preSymp(col))),
                    str(self.units[col])))
        exp.columns = names
        latex = exp.to_latex(index=False)\
            .replace("\\textbackslash ", "\\")\
            .replace("\\_", "_")\
            .replace("\\\\", "\\\\ \\hline")\
            .replace("\\$", "$").replace("e+", "e")\
            .replace("\\toprule", "\\hline")\
            .replace("\\midrule", "")\
            .replace("\\bottomrule", "")\
            .replace('\\textasciicircum', '^')\
            .replace(r'\{', '{')\
            .replace(r'\}', '}')\
            .replace('newton', "N")\
            .replace('$tau', r'$\tau')\
            .replace("pourcent", r"\%")\
            .replace("$omega", r"$\omega")
        if style == 'default':
            latex = latex\
                .replace("l" * len(exp.columns), "|" + "c|" * len(exp.columns))
        elif style == 'clean':
            latex = latex\
                .replace("l" * len(exp.columns), "c" * len(exp.columns))\
                .replace(r'\hline', '', 1)\
                .replace(r'\hline', r'\hline\hline', 1)
            latex = (latex[::-1].replace('enilh\\', '', 1))[::-1]
        with open(r"tables\{}.tex".format(name), "w+") as final:
            final.write(latex)

    def fixUnits(self):
        for col in self.data.columns[::2]:
            if self.units[col] == 0:
                self.units[col] = unit("1")
            const = self.units[col].exctractConstant()
            self.data[col] *= float(const)
            self.data[ef.delt(col)] *= float(const)
            self.units[col] = unit(str(self.units[col].symb / const))

    def importCol(self, name, dt, index=None):
        """\n

        imports a column from another table into itself

        (name) name of the column to import

        (dt) talbe to import from

        (index) index where you want to insert the column
        """
        if index is None:
            index = (len(self.data) - 1) / 2
        self.data.insert(index*2, name, dt[name], False)
        self.data.insert(index*2+1, ef.delt(name), dt[ef.delt(name)], False)

    def plot(self, xn, yn, fig_ax=None, polar=False, NoErrorBar=False,
             **kwargs):
        if fig_ax is None:
            fig_ax = plt.figure(), plt.gca()
        fig, ax = fig_ax

        dt = self.copy()
        if xn not in self.data.columns:
            dt.newCol(xn, xn)
        if yn not in self.data.columns:
            dt.newCol(yn, yn)
        if sum(dt[ef.delt(xn)]) == sum(dt[ef.delt(yn)]) == 0 or NoErrorBar:
            ax.plot(*dt[[xn, yn]].T, ".", **kwargs)
        else:
            ax.errorbar(*dt[[xn, yn, ef.delt(yn), ef.delt(xn)]].T, ".",
                        **kwargs)
        ax.set_xlabel(ef.ax_name(dt, xn))
        ax.set_ylabel(ef.ax_name(dt, yn))

    def fit(self, func, xn, yn, fig_ax=None, show=True, p0=None, maxfev=1000,
            fit_label='fit', fit_color=None, **kwargs):
        """\n

        return optimal parameters of a function fitted on data

        Applies a fit on givens colums, interpreted as a 2D graphics

        Parameters
        ----------
        (func) The function you want to fit your function to

        (xn) The data representing the x axis.

        (yn) The data representing the y axis.

        (show) If True, shows the data and the fit

        (p0) Initial parameters

        (maxfev) Maximum number of itteration to find optimals parameters
        """
        if fig_ax is None:
            fig_ax = plt.figure(), plt.gca()
        fig, ax = fig_ax
        dt = self.copy()
        if xn not in self.data.columns:
            dt.newCol(xn, xn)
        if yn not in self.data.columns:
            dt.newCol(yn, yn)
        popt, cov = curve_fit(func, *dt[[xn, yn]].T, p0=p0, maxfev=maxfev)
        res = np.array([np.array([popt, np.sqrt(cov.diagonal())]).T.flatten()])
        out = table('', data=res, AutoInsert=False)
        noms = list(inspect.signature(func).parameters.keys())
        noms = list(map(lambda x: sp.latex(sp.sympify(ef.preSymp(x))), noms))
        out.renameCols(noms[1:])
        lines = inspect.getsourcelines(func)

        if len(lines) == 2 and fit_label == 'fit':
            expr = sp.sympify(ef.preSymp(lines[0][-1].replace('return ', '')))
            fit_label = "$"+sp.latex(expr)+"$"
        if show:
            x = np.linspace(*ef.extrem(dt[xn]), 1000)
            ax.plot(x, func(x, *popt), label=fit_label, color=fit_color)
            dt.plot(xn, yn, fig_ax=fig_ax, **kwargs)
            ax.legend()
        return out

    def append(self, name, data, incert=None, pos=None, units="1"):
        if pos is None:
            pos = len(self.data.columns)
        if incert is None:
            incert = np.zeros(len(self))
        self.units[name] = unit(units)
        self.data.insert(pos, name, data)
        self.data.insert(pos + 1, ef.delt(name), incert)
        self.fixUnits()

    def copy(self):
        return deepcopy(self)

    def getStr(self, cols):
        """Combines data and uncertainty collumns. Used with makeGoodTable"""
        cols = list(cols)
        out = dict()
        for col in [str(i) for i in cols]:
            outCol = []
            # prendre toutes les valeurs + leurs delta individuellement
            for x in range(len(self.data[col])):
                # calculer le nombre de nombre de chiffres significatif requis
                val = self.data[col][x]
                d = self.data[ef.delt(col)][x]
                # formater la colonne en fonction des résultats
                if d == 0:
                    s = str(val)
                    while s.endswith('0') or s.endswith('.'):
                        s = s[:-1]
                    outCol.append("$" + s + "$")
                else:
                    outCol.append(
                        "$" + ef.precisionStr(val, d) + " \\pm " + "{:.1g}$"
                        .format(ef.roundUp(self.data[ef.delt(col)][x]))
                    )
            out[col] = outCol
        if len(list(out.values())[0]) == 1:
            out = {key: out[key][0] for key in list(out.keys())}
        if len(cols) == 1:
            out = out[col[0]]
        return(out)


def defVal(vals):
    for i in vals:
        if isinstance(i, str):
            a = vals[i].split(" ")
        else:
            a = vals[i]
        exval[i] = float(a[0])
        exval[ef.delt(i)] = float(a[1])
        exUnis[i] = unit(a[2])


def genVal(name, formula, units=None, NoUncert=False):
    expression = sp.sympify(ef.preSymp(formula))
    expression_orig = deepcopy(expression)
    pre_dict = dict(zip([f"dummy{x}" for x in range(
        len(expression.free_symbols))], map(str, expression.free_symbols)))
    expression = expression.subs({v: k for k, v in pre_dict.items()})
    lamb = sp.lambdify(tuple(expression.free_symbols), expression)
    if not NoUncert:
        expr_incert = ef.formule_incertitude(formula)
        pre_dict_i = dict(zip(
            [f"dummy{x}" for x in range(len(expr_incert.free_symbols))],
            map(str, expr_incert.free_symbols)
        ))
        expr_incert = expr_incert.subs(
            {v: k for k, v in pre_dict_i.items()})
        lamb_incert = sp.lambdify(tuple(expr_incert.free_symbols),
                                      expr_incert)
    # Le cacul est fait ligne par ligne, probablement très optimisable
    vals = exval
    nvals = {name: vals[pre_dict[name]] for name in pre_dict}
    delt_vals = {name: vals[pre_dict_i[name]] for name in pre_dict_i}
    exval[name], exval[ef.delt(name)] =\
        [lamb(**nvals), 0 if NoUncert else lamb_incert(**delt_vals)]

    if units:
            exUnis[name] = unit(units)
    else:
        exUnis[name] = lamb(*(
            exUnis[str(i)] for i in expression_orig.free_symbols))


def getVal(name):
    return exval[name], exval[ef.delt(name)], exUnis[name]


def getStrVal(name, latex=True):
    val = getVal(name)
    if latex:
        s =  f"${ef.valueStr(*val[0:2])}$"
        if not str(val[2]) == '1':
            s += r" ( $\rm {:}$ )".format(exUnis[name])
        return s
    else:
        s = f"{ef.valueStr(*val[0,1])}"
        if not str(val[2]) == '1':
            s += r" ( {:} )".format(tab.units[qty])
        return s
# %%
