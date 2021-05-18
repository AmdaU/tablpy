from . import extra_funcs as ef
#import extra_funcs as ef
import sympy as sp

prefixs = {"y": 1e-24, "z": 1e-21, "a": 1e-18, "f": 1e-15, "p": 1e-12,
           "n": 1e-9, "\\mu": 1e-6, "μ": 1e-6, "u": 1e-6, "m": 1e-3, "c": 1e-2,
           "d": 1e-1, "": 1, "h": 1e2, "k": 1e3, "M": 1e6, "G": 1e9, "T": 1e12,
           "P": 1e15, "B": 1e15, "E": 1e18, "Z": 1e21, "Y": 1e24}
SI = sp.symbols("g s m K A mol cd")

SIBT = {'Hz':       '1/s',
        'newton':   '1000*g*(m/s**2)',
        'Pa':       '1000*g/m/s**2',
        'J':        '1000*g*m**2/s**2',
        'W':        '1000*g*m**2/s**3',
        'C':        's*A',
        'V':        '1000*g*m**2/s**3/A',
        'F':        'A**2*s**4/(1000*g)/m**2',
        'ohms':     '1000*g*m**2/s**3/A**2',
        'S':        'A**2/(1000*g)/m**2*s**3',
        'H':        'm**2*(1000*g)/s**2/A**2',
       r'\AA':      '0.0000000001*m',
        'Wb':       '1000*m**2*g/s**2/A',
        'T':        '1000*g/A/s**2',
        'lm':       'cd',
        'lx':       'cd/m**2',
        'Bq':       '1/s',
        'Gy':       'm**2/s**2',
        'G':        '0.1*g/(A*s**2)',
        'lb':       '453.59237*g',
        'oz':       '28.3495231*g',
        'L':        'm**3',
        'percent':  '0.01',
        'deg':      '2 * 3.1415926535/360'
       }
SIBT = dict(zip(sp.symbols(list(SIBT.keys())), SIBT.values()))
# SIBT = dict(zip(SIB, sp.sympify(defs)))


class unit:
    def __init__(self, s):
        self.str = s
        s = s.replace("N", "newton").replace("Ω", "ohms")\
             .replace(r"\Omega", "ohms").replace("%", "percent")\
             .replace('°', 'deg').replace(r'^\circ', 'deg')
        self.symb = sp.sympify(ef.preSymp(s))
        self.SIval = self.symb
        self.prefix = None
        for var in self.SIval.free_symbols:
            if len(str(var)[1:]):
                if sp.symbols(str(var)[1:]) in list(SIBT.keys()):
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
             .replace('deg', r'^\circ')

    def __repr__(self):
        return str(self.symb)

    def __pow__(self, other):
        return unit(str(self.symb**other))

    def to(self, nunit):
        """\n

        Convertes units to others and outputs de ratio

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
            print(f"Conversion failed due to {self.str} and "
                  f"{nunit.str} being incompatible")
            factor = 1
        return 1 / factor

    def exctractConstant(self):
        const = self.symb
        for var in const.free_symbols:
            const = const.subs(var, 1)
        return const
