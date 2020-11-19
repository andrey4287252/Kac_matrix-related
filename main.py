import sympy as sp
import sys

debug = False

def Kac_Matrix_element(m, l, x, y, Liouville_parametrization=False):
    '''
    Agregates given partitions m and l to one list and computes corresponding element of Kac-Shapovalov matrix
    :param m: first partition
    :param l: second partition
    :param x: symbol corresponding to $\alpha$/energy with respect to parametrization
    :param y: symbol corresponding to $\beta$/charge with respect to parametrization
    :param Liouville_parametrization: if true expresses h and c in terms of a and b (with respect to known parametrization)
    '''
    L = list(map(int, l)) + list(map(lambda x: -int(x),  m))[::-1]
    if Liouville_parametrization:
        return Auxilary_Kac_Matrix_element(L, x*((y + 1/y)-x), 1+6*(y + 1/y)**2)
    else:
        return Auxilary_Kac_Matrix_element(L, x, y)

def Auxilary_Kac_Matrix_element(L, h, c):
    '''
    Returns element of Kac-Shapovalov matrix indexed by partitions agregated in L
    :param L: list concatenated from pair of indexing partitions, where second inversed in sign and order.
    :param h: symbol corresponding to energy 
    :param c: symbol corresponding to charge
    :return: symbolical expression for element of form indexed by pair of partitions
    '''
    if debug:
        print("Partitons: ", L)
    if len(L) == 0:
       return 1
    if L[len(L) - 1] > 0:
        return 0
    elif L[0] < 0:
        return 0
    elif L[len(L) - 1] == 0:
        return Auxilary_Kac_Matrix_element(L[:len(L) - 1], h, c)*h
    elif L[0] == 0:
        return Auxilary_Kac_Matrix_element(L[1:], h, c)*h
    else:
        for i in range(len(L) - 1):
            if L[i] >= 0 and L[i + 1] < 0:
                if debug:
                    print("Commutes: ", i, " and ", i + 1)
                    print(Auxilary_Kac_Matrix_element(L[:i] + [L[i] + L[i + 1]] + L[i + 2:]), "appears")
                answ = (L[i] - L[i + 1]) * Auxilary_Kac_Matrix_element(L[:i] + [L[i] + L[i + 1]] + L[i + 2:], h, c) + Auxilary_Kac_Matrix_element(L[:i] + [L[i + 1], L[i]] + L[i + 2:], h, c)
                if L[i] == -L[i + 1]:
                    answ += (L[i]**3 - L[i])*c*Auxilary_Kac_Matrix_element(L[:i] + L[i + 2:], h, c)/12
                return answ

if '-d' in sys.argv[1:]:
    debug = True

def Kac_matrix(N, x, y, Liouville_parametrization=False):
    '''
    Gives Kac-Shapovalov matrix on N-th level.
    :param N: size of the partitions enumerated in matrix
    :param x: symbol corresponding to $\alpha$/energy with respect to parametrization
    :param y: symbol corresponding to $\beta$/charge with respect to parametrization
    :param Liouville_parametrization: if true expresses h and c in terms of a and b (with respect to known parametrization)
    :return: matrix of Kac-Shapovalov form on n-th level
    '''
    partitions = list(sp.utilities.iterables.ordered_partitions(N))[::-1]
    if debug:
        print("partitions: ", partitions)
    return sp.Matrix(len(partitions), len(partitions), lambda i,j: sp.factor(Kac_Matrix_element(partitions[i], partitions[j], x, y, Liouville_parametrization)))

Liouville_parametrization = ('-L' in sys.argv[1:])
if Liouville_parametrization:
    print("Enter a")
    a = input()
    if not a.isnumeric():
        a = sp.symbols(a)
    else:
        a = float(a)
    print("Enter b")
    b = input()
    if not b.isnumeric():
        b = sp.symbols(b)
    else:
        b = float(b)

    print("Enter level")
    print(Kac_matrix(int(input()), a, b, True))
else:
    print("Enter symbolycal expression for weight (energy)")
    h = input()
    if not h.isnumeric():
        h = sp.symbols(h)
    else:
        h = float(h)
    print("Enter symbolical expression for charge")
    c = input()
    if not c.isnumeric():
        c = sp.symbols(c)
    else:
        c = float(c)
    print("Enter level")
    print(Kac_matrix(int(input()), h, c, False))
