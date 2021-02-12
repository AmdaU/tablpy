import pandas as pd
import sympy as sp
import numpy as np
import math
import re
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

delt = lambda s: r"\Delta" + str(s) if '\\' in str(s) else r"\Delta\!" + str(s)

_Data_ = {}

prefixs = {"y": 1e-24, "z": 1e-21, "a": 1e-18, "f": 1e-15, "p": 1e-12, "n": 1e-9,
            "\\mu": 1e-6, "μ": 1e-6, "u": 1e-6, "m": 1e-3, "c": 1e-2, "d": 1e-1,
           "": 1, "h": 1e2, "k": 1e3, "M": 1e6, "G": 1e9, "T": 1e12, "P": 1e15,
            "B": 1e15, "E": 1e18, "Z": 1e21, "Y": 1e24}
SI = sp.symbols("g s m K A mol cd")
SIB = sp.symbols("Hz newton Pa J W C V F ohms S H Wb T H °C lm lx Bq Gy Sv kat pourcent")
defs = ["1/s", "1000*g*(m/s**2)", "1000*g/m/s**2", "1000*g*m**2/s**2", "1000*g*m**2/s**3", "s*A", "1000*g*m**2/s**3/A",
        "A**2*s**4/(1000*g)/m**2", "1000*g*m**2/s**3/A**2", "A**2/(1000*g)/m**2*s**3", "m**2*(1000*g)/s**2/A**2"]

SIBT = dict(zip(SIB, sp.sympify(defs)))

exval,exUnis = {},{}

class unit:
    def __init__(self, s):
        self.str = s
        s = s.replace("N", "newton").replace(
            "Ω", "ohms").replace("\Omega", "ohms").replace("%","pourcent")
        self.symb = sp.sympify(s)
        self.SIval = self.symb
        self.prefix = None
        for var in self.SIval.free_symbols:
            if len(str(var)[1:]):
                if sp.symbols(str(var)[1:]) in SIB:
                    self.SIval *= sp.symbols(str(var)
                                             [1:]) / var * prefixs[str(var)[0]]
        for var in self.SIval.free_symbols:
            if var in SIBT:
                self.SIval = self.SIval.subs(var, SIBT[var])
        for var in self.SIval.free_symbols:
            if len(str(var)[1:]):
                if sp.symbols(str(var)[1:]) in SI:
                    self.SIval *= sp.symbols(str(var)
                                             [1:]) / var * prefixs[str(var)[0]]

    def __add__(self, other):
        if type(other) == unit:
            if self.SIval == other.SIval:
                return self
        print("Cant add these two, boi")

    def __sub__(self, other):
        return self.__add__(other)

    def __neg__(self):
        return self

    def __mul__(self, other):
        if type(other) == unit:
            return unit(str(self.symb * other.symb))
        elif type(other) in (int, float):
            return self

    def __rmul__(self, other):
        if type(other) == unit:
            return unit(str(self.symb * other.symb))
        elif type(other) in (int, float):
            return self

    def __truediv__(self, other):
        if type(other) == unit:
            return unit(str(self.SIval / other.SIval))
        elif type(other) in (int, float):
            return self

    def __rtruediv__(self, other):
        if type(other) == unit:
            return unit(str(self.symb / other.symb))
        elif type(other) in (int, float):
            return unit(str(1 / self.symb))

    def __str__(self):
        return sp.latex(self.str)

    def __repr__(self):
        return str(self.symb)

    def __pow__(self, other):
        return unit(str(self.symb**other))

    def to(self, nunit):
        """
        Convertis des unités en d'autres unités
            Agit sur l'objet lui-même et retourne le facteur de converstion

            ex:
            >>> a = unit("m/s**2")
            >>> print(a)
            ...m/s**2

            >>> a.to("N/g")
            ... 1000
            >>> print(a)
            >>> N/kg
        """
        if type(nunit) == type(" "):
            nunit = unit(nunit)
        factor = nunit.SIval / self.SIval
        if not len(factor.free_symbols):
            self.SIval = nunit.SIval
            self.str = nunit.str
            self.symb = nunit.symb
        else:
            print(f"La converstion d'unité à échoué car les {self.str} et les {nunit.str} sont incompatibles")
            factor = 1
        return 1 / factor

    def exctractConstant(self):
        const = self.symb
        for var in const.free_symbols:
            const = const.subs(var, 1)
        return const

# %%

class tableau:
    def __init__(self, datname, units={}, formules={}, sheet=None, AutoInsert = False):
        data = None
        if type(dataname) != str:
            data = datname
        self.units = units
        if np.array(data != None).any:
            self.data = pd.DataFrame(data)
        elif os.path.isfile(datname + ".xlsx"):
            self.data = pd.read_excel(datname + ".xlsx", sheet_name=sheet)
            noms = list(self.data.columns)[::2]
        elif os.path.isfile(datname + ".csv"):
                self.data = pd.DataFrame(np.genfromtxt(datname + ".csv",delimiter = ",",names = True))
        elif os.path.isfile(datname + ".txt"):
                self.data = pd.DataFrame(np.genfromtxt(datname + ".txt",delimiter = ",",names = True))
        else:
            print("Le fichier {} n'as pas été trouvé :(".format(datname))
        #self.renomerCols(list(noms))
        prank = self.giveUnits(units)
        #self.vars = sp.symbols(list(self.data.columns))
        self.formules = {}

        if AutoInsert:
            l = len(self.data.columns)
            lol = True
            while lol:
                lol = False
                for i in range(len(self.data.columns)-1):
                    if not isUncertain(self.data.columns[i]) and not isUncertain(self.data.columns[i+1]):
                        self.data.insert(i+1,delt(self.data.columns[i]),[0 for i in range(len(self.data))])
                        lol = True

            if not isUncertain(self.data.columns[-1]):
                self.data.insert(len(self.data.columns),delt(str(self.data.columns[-1])),[0 for i in range(len(self.data))])

    def __repr__(self):
        return self.data.to_markdown()

    def __getitem__(self, x):
        return np.array(self.data[x])

    def __iter__(self):
        yield from self.data.columns

    def nouvCol(self, nom, rexpression, extra=None, pos=None, NoIncert = False):
        if pos == None:
            pos = len(self.data.columns)
        expression = sp.sympify(preSymp(rexpression))
        self.formules[nom] = expression
        out = np.empty((len(self.data), 2))
        lamb = sp.lambdify(expression.free_symbols,expression)
        # Le cacul est fait ligne par ligne, probablement très optimisable
        for i in range(len(self.data)):
            expr = expression
            formules = formule_incertitude(rexpression)
            for j in self.data.iloc[[i]]:
                # (valeurs) substitue chaque valeurs impliqué une à la fois, à optimiser encore une fois
                expr = expr.subs(j, float(self.data.iloc[[i]][j]))
                formules = formules.subs(
                    j, float(self.data.iloc[[i]][j]))  # (incertitude)
            for var in list(expr.free_symbols):
                expr = expr.subs(var, exval[str(var)])
            for var in list(formules.free_symbols):
                formules = formules.subs(var, exval[str(var)])
            if NoIncert:
                out[i] = [expr.evalf(), 0]
            else:
                out[i] = [expr.evalf(), formules.evalf()]
        # Nomme les colonnes et les ajoute au tableau
        info = {}
        for i in exUnis.keys():
            if i in map(str, expression.free_symbols):
                info[i] = exUnis[i]

        func = sp.lambdify(expression.free_symbols, expression, "math")
        self.units[nom] =func(*(self.units[str(i)] if str(i) in self.units else exUnis[str(i)]\
            for i in self.formules[nom].free_symbols))
        self.units[delt(nom)] = self.units[nom]
        self.data.insert(pos, nom, out[:, 0])
        self.data.insert(pos + 1, delt(nom), out[:, 1])

    def delCol(self,nom):
        self.data = self.data.drop(columns=[nom,delt(nom)])

    def giveUnits(self, units):
        for col in self.data.columns:
            try:
                if type(units[col]) == type(""):
                    units[col] = unit(units[col])
                self.units[col] = units[col]
            except:
                try:
                    self.units[col] = self.units[col]
                except:
                    self.units[col] = unit("1")

    def changeUnits(self, units):
        for col in units:
            fact = self.units[col].to(units[col])
            self.data[col] *= fact
            self.data[delt(col)] *= fact

    def renomerCols(self, noms):
        vars = sp.symbols(noms)
        try:
            self.data.columns = "".join([f"{i} {sp.symbols(delt(str(i)))} " for i in vars])[:-1].split(" ")
        except:
            vars = [vars]
            self.data.columns = "".join([f"{i} {sp.symbols(delt(str(i)))} " for i in vars])[:-1].split(" ")

    def squish(self):  # Combine les colonnes de valeurs et d'incertitude et applique de la mise en forme des données
        out = pd.DataFrame([])
        for col in [i for i in self.data if "Delta" not in i]:
            outCol = []
            # prendre toutes les valeurs + leurs delta individuellement
            for x in range(len(self.data[col])):
                # calculer le nombre de nombre de chiffres significatif à accorder à ceux ci
                val = self.data[col][x]
                d = self.data[delt(col)][x]
                # formater la colonne en fonction des résultats trouvée précédemment
                if d == 0:
                    s = str(val)
                    while '.' in s and s.endswith('0') or s.endswith('.'):
                        s = s[:-1]
                    outCol.append("$" + s + "$")
                elif val == 0:
                    rd = roundUp(d)
                    srd = str(rd)
                    if srd[0:1] == "0.":
                        l = len(srd[2:])
                    else:
                        l=0
                    outCol.append("$0" +"."*(l>0) +"0"*l + r"\pm "+ srd + "$")
                else:
                    miam = "{{:.{}g}}".format(-math.ceil(-sp.log(abs(val), 10)) +
                                              math.ceil(-sp.log(d, 10)) + 1).format(val)
                    outCol.append("$" + miam + " \\pm " +
                              "{:.1g}".format(roundUp(self.data[delt(col)][x])) + "$")
            out["$" + col + "$"] = outCol
        return(out)

    # Export un ficher .tex en faisant toutes les modifications nécessaire pour que le code soit compris par latex
    def makeGoodTable(self,nom, unite=None):
        try:
            os.mkdir("tableaux")
        except:
            pass
        self.fixUnits()
        exp = self.squish()
        noms = []
        for col in self.data.columns[::2]:
            if str(self.units[col].SIval) in ["1", '1.00000000000000']:
                noms.append("${}$".format(col))
            else:
                noms.append("${}$ ({})".format(col, self.units[col].symb))
        exp.columns = noms
        latex = exp.to_latex(index=False)\
            .replace("\\textbackslash ", "\\")\
            .replace("\\_", "_")\
            .replace("\\\\", "\\\\ \\hline")\
            .replace("\\$", "$").replace("e+", "e")\
            .replace("\\toprule", "\\hline")\
            .replace("\\midrule", "")\
            .replace("\\bottomrule", "")\
            .replace('\\textasciicircum', '^')\
            .replace('\{', '{')\
            .replace('\}', '}')\
            .replace('newton',"N")\
            .replace('$tau',r'$\tau')\
            .replace("l" * len(exp.columns), "|" + "c|" * len(exp.columns))\
            .replace("pourcent","\%")\
            .replace("$omega","$\omega")
        with open(f"tableaux\\{nom}.tex", "w+") as final:
            final.write(latex)

    def fixUnits(self):
        for col in self.data.columns:
            const = self.units[col].exctractConstant()
            self.data[col] *= const
            self.units[col] = unit(str(self.units[col].symb / const))

    def errorbar(self,kx,ky,show = True):
        nouvTableau = tableau("temp","rien",data =self.data)
        if kx not in nouvTableau:
            nouvTableau.nouvCol(kx,kx)
        if ky not in nouvTableau:
            nouvTableau.nouvCol(ky,ky)
        x=np.array(list(map(float,nouvTableau[kx])))
        y=np.array(list(map(float,nouvTableau[ky])))
        yerr = np.array(list(map(float,nouvTableau[delt(ky)])))
        xerr = np.array(list(map(float,nouvTableau[delt(kx)])))
        if show:
            plt.errorbar(x,y,yerr = yerr,xerr=xerr,fmt=".",color = "darkorchid",label= "données expérimentales")
            plt.xlabel("$"+kx+"$ "+"( "+str(self.units[kx])+" )")
            plt.ylabel("$"+ky+"$ "+"( "+str(self.units[ky])+" )")
        return x,y,yerr,xerr

    def linFit(self,kx,ky,show = True):
        x,y,yerr,xerr = self.errorbar(kx,ky,show=show)
        popt,pcov = curve_fit(lin,x,y)
        if show:
            xf = np.linspace(min(x),max(x),2)
            plt.plot(x,lin(x,*popt),color = 'darkturquoise',label= "fit linéaire")
        return popt,pcov

    def importCol(self,name,df,index = None):
        if index == None:
            index = (len(self.data)-1)/2
        self.data.insert(index*2,name,df[name],False)
        self.data.insert(index*2+1,delt(name),df[delt(name)],False)

    def __add__(self,other,nom = None):
        if nom == None:
            nom = self.nom
        if type(other) == tableau:
            return tableau(nom,"nothing_cus_this_is_dumb",data = pd.concat([self.data,other.data]))
        else:
            print("oopsie doopsie, that did nothing")
            return self

    def sort(self,col, ascending=True):
        self.data = self.data.sort_values(by=col,ascending=ascending)


# %%

lin = lambda x,a,b: a*x +b

def preSymp(expr):
    filter = r"\\?\b(?:(?!sin|cos|tan|ln)[a-zA-Z]+[a-zA-Z0-9_]*)'?"
    new  = tuple(f'symbols("{i.group()}")' for i in re.finditer(filter,expr))
    for i in new:
        expr = re.sub(filter,r"{:}",expr)
    expr = expr.format(*new)
    return expr

def roundUp(num, n=1):  # arrondis vers le haut à la bonne décimale
    pow = math.floor(sp.log(num, 10).evalf()) + n-1  # Position de la décimla
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
    eq = sp.sympify(preSymp(eq))
    variables = list(eq.free_symbols)  # liste de tout les variables
    # liste de touts les incertitudes asssociées au variables
    uncertain = [sp.symbols(delt(x)) for x in variables]
    fIncert = sp.sqrt(sum(
        [(sp.diff(eq, variables[i]) * uncertain[i])**2 for i in range(len(variables))]))
    return sp.simplify(fIncert)


# permet de reonommer les colonnes pour pouvoir travailler avec des variables plus simples

def defVal(vals):
    for i in vals:
        a = vals[i].split(" ")
        exval[i] = a[0]
        exval[delt(i)] = a[1]
        exUnis[i] = unit(a[2])

def exctractCsv(rep):
    for i in os.listdir(rep):
        if os.path.isdir(rep+"\\"+i):
            exctractCsv(rep+"\\"+i)
        elif ".csv" in i.lower():
            add = ""
            _Data_[rep+"\\"+i] = np.genfromtxt(rep+"\\"+i,delimiter = ",")

def exctractData(rep):
    for i in os.listdir(rep):
        com = ""
        i = i.lower()
        if os.path.isdir(rep+"\\"+i):
            exctractData(rep+"\\"+i)
        elif ".csv" in i:
            com = "_Data_['{0}'] = np.genfromtxt('{1}',delimiter = ',')"\
            .format(i.replace(".csv","").replace(".",""),rep+"\\"+i)
            print("found: {}!".format(i))
        elif ".xlsx" in i:
            com = "_Data_['{0}'] = pd.read_excel('{1}')"\
            .format(i.replace(".xlsx","").replace(".",""),rep+"\\"+i)
            print("found: {}!".format(i))
        com = com.replace("\\","\\\\")
        exec(com)

def clearfile(file):
    for i in os.listdir(file):
        if os.path.isdir(file+"/"+i):
            clearfile(file+"/"+i)
            os.rmdir(file+"/"+i)
        else:
            os.remove(file+"/"+i)

def extrem(l):
    return float(min(l)),float(max(l))

def isUncertain(string):
    marqeurs = ['delta',
                'Delta',
                'Δ',
                'Î”']
    return any([marqeur in str(string) for marqeur in marqeurs])

def getTexIncert(x,special = False):
    out = sp.latex(formule_incertitude(x)).replace("\!"," ")
    if special == True:
        filter = r"(\\left)(\()\2"
        new  = tuple(f'a' for i in re.finditer(filter,expr))
        for i in new:
            expr = re.sub(filter ,r"{:}\2",expr)
        out = out.replace("\left(","\qty")
    print()

def fastPlot(dir, save = False,xlabel = None,ylabel =None):
    bip = dir.split("\\")
    if bip == []: bip = [""]
    sumStr = lambda l: l[0] + sumStr(l[1:]) if len(l) > 0 else ""
    bloop = sumStr(["\\" + b for b in bip[:-1]])
    if save and not os.path.isdir("graphiques" + bloop):
        os.mkdir("graphiques" + bloop)
    if(os.path.isdir(dir)):
        for i in os.listdir(dir):
            print(i)
            fastPlot(dir + "\\" + i,save,xlabel,ylabel)
    elif ".csv" in dir or ".txt" in dir:
        array = np.genfromtxt(dir,skip_header=True,delimiter = "\t")
        plt.plot(array[:,0],array[:,1],color = "darkturquoise")
        if xlabel:plt.xlabel(xlabel)
        if ylabel:plt.ylabel(ylabel)
        if save: plt.savefig("graphiques\\" + dir.replace(".csv",".png").replace(".txt",".png"))
        plt.show()
