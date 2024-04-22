
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

    def ge(self, name):
        """
        Gets the expression with the specified name from the environment.

        Parameters:
            name (str): Name of the expression to retrieve.

        Returns:
            sympy.Expr or None: The expression with the specified name, or None if not found.
        """
        return self.expressions.get(name)

    def d(self):
        """
        Displays all symbols and expressions in the environment.
        """
        self.ds()
        self.de()

    def ds(self):
        """
        Displays all symbols in the environment.
        """
        for symbol, value in self.symbols.items():
            if value is None:
                value_str = "?"
            else:
                value_str = value
            display(Math(fr"{symbol} = {value_str}"))

    def de(self):
        """
        Displays all expressions in the environment.
        """
        for name, expr in self.expressions.items():
            variables = ', '.join(str(var) for var in expr.free_symbols)
            latex_expr = latex(expr, mode='plain')
            display(Math(fr"{name}({variables}) = {latex_expr}"))
