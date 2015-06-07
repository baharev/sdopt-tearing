import sympy as sy

def expr_to_dotgraph(name, expr):
    from sympy.printing.dot import dotprint
    with open(name + ".dot", "w") as dot_file:
        dot_file.write(dotprint(expr))
        
if __name__ == '__main__':
    x, y, alpha = sy.Symbol('x'), sy.Symbol('y'), sy.Symbol('alpha')
    expr = sy.Eq(y, alpha*x/(1+(alpha-1)*x))
    expr_to_dotgraph('sympy_tree', expr)
    # dot -Tsvg -o sympy_tree.svg sympy_tree.dot
