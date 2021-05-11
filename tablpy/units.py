import sympy as sp
from . import extra_funcs as ef

prefixs = {"y": 1e-24, "z": 1e-21, "a": 1e-18, "f": 1e-15, "p": 1e-12, "n": 1e-9,
           "\\mu": 1e-6, "μ": 1e-6, "u": 1e-6, "m": 1e-3, "c": 1e-2, "d": 1e-1,
           "": 1, "h": 1e2, "k": 1e3, "M": 1e6, "G": 1e9, "T": 1e12, "P": 1e15,
           "B": 1e15, "E": 1e18, "Z": 1e21, "Y": 1e24}
SI = sp.symbols("g s m K A mol cd")
SIB = sp.symbols("Hz newton Pa J W C V F ohms S H Wb T H °C lm lx Bq Gy Sv kat"
                 "percent")
defs = ["1/s", "1000*g*(m/s**2)", "1000*g/m/s**2", "1000*g*m**2/s**2",
        "1000*g*m**2/s**3", "s*A", "1000*g*m**2/s**3/A",
        "A**2*s**4/(1000*g)/m**2", "1000*g*m**2/s**3/A**2",
        "A**2/(1000*g)/m**2*s**3", "m**2*(1000*g)/s**2/A**2"]

SIBT = dict(zip(SIB, sp.sympify(defs)))


class unit:
    def __init__(self, s):
        self.str = s
        s = s.replace("N", "newton").replace("Ω", "ohms")\
             .replace(r"\Omega", "ohms").replace("%", "percent")\
             .replace('°', 'deg').replace('^\circ', 'deg')
        self.symb = sp.sympify(ef.preSymp(s))
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
        raise(Exception(f"Cannot add units {self} and {other}"))

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
        return sp.latex(self.symb).replace("newton", 'N')\
             .replace("ohms", r"\Omega").replace("percent", "%")\
             .replace('deg', '^\circ')

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
        if isinstance(nunit, str):
            nunit = unit(nunit)
        factor = nunit.SIval / self.SIval
        if not len(factor.free_symbols):
            self.SIval = nunit.SIval
            self.str = nunit.str
            self.symb = nunit.symb
        else:
            print(f"La converstion d'unité à échoué car les {self.str} et les"
                  f"{nunit.str} sont incompatibles")
            factor = 1
        return 1 / factor

    def exctractConstant(self):
        const = self.symb
        for var in const.free_symbols:
            const = const.subs(var, 1)
        return const
