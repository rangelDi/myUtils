import sympy
from sympy import latex
from IPython.display import display, Math

class env:

    def __init__(self, symbols=None):
        """
        Initializes the environment with a list of symbols.

        Parameters:
            symbols (list or None): List of tuples containing sympy symbols and their numeric values.
                                    If None, the environment starts empty.
        """
        self.symbols = {}
        self.expressions = {}
        if symbols is not None:
            for symbol, value in symbols:
                self.s(symbol, value)

    def s(self, name, value=None):
        """
        Adds a symbol to the environment. When using LaTex syntax, pass the name as a raw string.
    
        Parameters:
            name (str): Name of the symbol.
            value (float or None): Numeric value of the symbol. Default is None.
                                If None, the value is set to None.
    
        Returns:
            sympy.Symbol: The created symbol.
        """
        symbol = sympy.Symbol(name)
    
        # Check if an expression with the same name exists
        if name in self.expressions:
            del self.expressions[name]
    
        self.symbols[symbol] = value
        return symbol
    
    def e(self, name, expr):
        """
        Adds an expression to the environment.
    
        Parameters:
            name (str): Name of the expression.
            expr (sympy.Expr): The symbolic expression.
    
        Returns:
            sympy.Expr: The added expression.
        """
        expr_symbols = expr.free_symbols
    
        # Check if a symbol with the same name exists and remove it
        for symbol, value in self.symbols.items():
            if str(symbol) == name:
                del self.symbols[symbol]
                break
    
        # Store the expression
        self.expressions[name] = expr
    
        # Add symbols from the expression
        for symbol in expr_symbols:
            if symbol not in self.symbols:
                self.symbols[symbol] = None
    
        return expr

    def a(self, name, value_or_expr):
        """
        Adds a single symbol or expression to the environment.

        Parameters:
            name (str): Name of the symbol or expression.
            value_or_expr (float, int, sympy.Expr): Numeric value for a symbol, or a sympy expression.

        Returns:
            sympy.Symbol or sympy.Expr: The added symbol or expression.
        """
        if isinstance(value_or_expr, (float, int)) or value_or_expr is None:
            return self.s(name, value_or_expr)
        elif isinstance(value_or_expr, sympy.Expr):
            return self.e(name, value_or_expr)
        else:
            raise ValueError("Invalid value_or_expr format.")

    def ad(self, items):
        """
        Adds symbols and expressions to the environment.

        Parameters:
            items (list): List containing tuples of the form (str, float) or (str, sympy expression).

        Returns:
            tuple: A tuple containing the sympy objects created from the items.
        """
        sympy_objects = []
        for item in items:
            if isinstance(item, tuple) and len(item) == 2 and isinstance(item[0], str):
                if isinstance(item[1], (float, int)) or item[1] is None:
                    sympy_objects.append(self.s(item[0], item[1]))
                elif isinstance(item[1], sympy.Expr):
                    sympy_objects.append(self.e(item[0], item[1]))
                else:
                    raise ValueError("Invalid item format.")
            else:
                raise ValueError("Invalid item format.")
        # Retain the original functionality by not returning anything if no objects are added
        if sympy_objects:
            return tuple(sympy_objects)

    def eval(self, expression=None):
        """
        Evaluates expressions in the environment numerically with the given substitutions.

        Parameters:
            subs (dict or None): Dictionary of substitutions for symbols in the expressions.
                                 If provided, updates the values of corresponding symbols.
            name (str or None): Optional name of the expression to evaluate.
                                 If provided, only evaluates the expression with that name.
                                 If not provided, evaluates all expressions.

        Returns:
            dict: Dictionary containing the names of expressions as keys and their numerical values as values.
        """
        evaluated_expressions = {}
        subs = {}

        for symbol in self.symbols:
            if self.symbols[symbol] is not None:
                subs[symbol] = self.symbols[symbol]

        if expression is not None:
            if expression in self.expressions:
                evaluated_expressions[expression] = self.expressions[expression].evalf(subs=subs)
            else:
                raise ValueError(f"No expression found with name '{name}'.")
        else:
            for name, expr in self.expressions.items():
                evaluated_expressions[name] = expr.evalf(subs=subs)

        for expression, value in evaluated_expressions.items():
            if isinstance(value,sympy.core.numbers.Float):
                self.a(expression,float(value))
            else:
                self.a(expression,value)

        for symbol in self.symbols:
            aVal = self.symbols[symbol] 
            if aVal is not None:
                evaluated_expressions[str(symbol)] = float(aVal)

        return evaluated_expressions


    def p(self,name = None):
        """
        Prints all, some, or one symbol and/or expression of the environment.
        """

        if name is None:
            self.ps()
            self.pe()
            return
        elif isinstance(name,str):
            self.ps(name)
            self.pe(name)
            return
        elif isinstance(name,list) and len(name) > 0:
            print("Symbols:")
            for n in name:
                self.ps(n)
            print("\nExpressions:")
            for n in name:
                self.pe(n)

    def ps(self,name = None):
        """
        Prints all, some, or one symbol of the environment.
        """
        if name is None:
            for symbol in self.symbols:
                print(f"{symbol} = {self.symbols[symbol]}")
        elif isinstance(name,str):
            for symbol in self.symbols:
                if str(symbol) == name:
                    print(f"{symbol} = {self.symbols[symbol]}")
        elif isinstance(name,list) and len(name) > 0: 
            for n in name:
                for symbol in self.symbols:
                    if str(symbol) == n:
                        print(f"{symbol} = {self.symbols[symbol]}")

    def pe(self, name = None):
        """
        Prints all, some, or one expression of the environment.
        """
        if name is None:
            for n, expr in self.expressions.items():
                print(f"{n} = {expr}")
        elif isinstance(name,str):
            if name in self.expressions:
                expr = self.expressions[name]
                print(f"{name} = {expr}")
        elif isinstance(name,list) and len(name) > 0:
            for na in name:
                if na in self.expressions:
                    expr = self.expressions[na]
                    print(f"{na} = {expr}")

    def g(self):
        """
        Gets all symbols and expressions in the environment.

        Returns:
            tuple: A tuple containing two tuples.
                   The first tuple contains the symbols in the environment,
                   and the second tuple contains the expressions in the environment.
        """
        return tuple(self.symbols.keys()), tuple(self.expressions.items())

    def gs(self):
        """
        Gets all symbols in the environment.

        Returns:
            tuple: A tuple containing all symbols in the environment.
        """
        return tuple(self.symbols.keys())

    def val(self, name):

        for symbol in self.symbols:
            if str(symbol) == name:
                return self.symbols[symbol]


    def ge(self, name):
        """
        Gets the expression with the specified name from the environment.

        Parameters:
            name (str): Name of the expression to retrieve.

        Returns:
            sympy.Expr or None: The expression with the specified name, or None if not found.
        """
        return self.expressions.get(name)

    def d(self, name = None):
        """
        Displays all, some, or one symbol and/or expression of the environment.
        """
        if name is None:
            self.ds()
            self.de()
            return
        elif isinstance(name,str):
            self.ds(name)
            self.de(name)
            return
        elif isinstance(name,list) and len(name) > 0:
            for n in name:
                self.ds(n)
            for n in name:
                self.de(n)

    def ds(self, name = None):
        """
        Displays all, some, or one symbol of the environment.
        """
        if name is None:
            for symbol, value in self.symbols.items():
                if value is None:
                    value_str = "?"
                else:
                    value_str = value
                display(Math(fr"{symbol} = {value_str}"))
        elif isinstance(name,str):
            for symbol in self.symbols:
                if str(symbol) == name:
                    value = self.symbols[symbol]
                    if value is None:
                        value_str = "?"
                    else:
                        value_str = value
                    display(Math(fr"{symbol} = {value_str}"))
        elif isinstance(name,list) and len(name) > 0: 
            for n in name:
                for symbol in self.symbols:
                    if str(symbol) == name:
                        value = self.symbols[symbol]
                        if value is None:
                            value_str = "?"
                        else:
                            value_str = value
                        display(Math(fr"{symbol} = {value_str}"))

    def de(self, name = None):
        """
        Displays all, some, or one expression of the environment.
        """

        if name is None:
            for n, expr in self.expressions.items():
                variables = ', '.join(str(var) for var in expr.free_symbols)
                latex_expr = latex(expr, mode='plain')
                display(Math(fr"{n}({variables}) = {latex_expr}"))
            return
        elif isinstance(name,str):
            if name in self.expressions:
                expr = self.expressions[name]
                variables = ', '.join(str(var) for var in expr.free_symbols)
                latex_expr = latex(expr, mode='plain')
                display(Math(fr"{name}({variables}) = {latex_expr}"))
        elif isinstance(name,list) and len(name) > 0:
            for na in name:
                if na in self.expressions:
                    expr = self.expressions[na]
                    variables = ', '.join(str(var) for var in expr.free_symbols)
                    latex_expr = latex(expr, mode='plain')
                    display(Math(fr"{na}({variables}) = {latex_expr}"))

    def r(self, name):
        """
        Removes a symbol or expression from the environment.
    
        Parameters:
            name (str): Name of the symbol or expression to remove.
        """
        # Check if the name corresponds to a symbol
        for symbol, value in self.symbols.items():
            if str(symbol) == name:
                del self.symbols[symbol]
                return
    
        # Check if the name corresponds to an expression
        if name in self.expressions:
            del self.expressions[name]
        else:
            print(f"No symbol or expression found with the name '{name}'.")
