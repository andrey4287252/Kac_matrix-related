import sympy as sp
import sys

debug = False


def Kac_Matrix_element(m, l, x, y):
    '''
    Agregates given partitions m and l to one list and computes corresponding element of Kac-Shapovalov matrix
    :param m: first partition
    :param l: second partition
    :param x: symbol corresponding energy
    :param y: symbol corresponding charge
    :return: symbolycal expression for Kac-Shapovalov matrix element corresponding partitions given
    '''
    L = list(map(int, l)) + list(map(lambda x: -int(x),  m))[::-1]
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

def Kac_matrix(N, x, y):
    '''
    Gives Kac-Shapovalov matrix on N-th level.
    :param N: size of the partitions enumerated in matrix
    :param x: symbol corresponding to energy
    :param y: symbol corresponding to charge
    :return: matrix of Kac-Shapovalov form on n-th level
    '''
    partitions = sorted([p[::-1] for p in list(sp.utilities.iterables.ordered_partitions(N))], reverse=True)
    if debug:
        print("partitions: ", partitions)
    return sp.Matrix(len(partitions), len(partitions), lambda i,j: sp.simplify(Kac_Matrix_element(partitions[i], partitions[j], x, y)))


def Liouville_parametrization_inp():
    '''
    Manage input if Liouville parametrization required
    :return: N - number of level considering, h - 
    '''
    b = sp.symbols('b')
    if '-s' in sys.argv[1:]:
        print("Enter m")
        m = int(input())
        print("Enter n")
        n = int(input())
        N = m*n
        a = -((m - 1)*b/2 + (n - 1)/(2*b))
    else:
        a = sp.symbols('a')
        print("Enter level")
        N = int(input())
    h = a*((b + 1/b)-a)
    c = 1+6*(b + 1/b)**2
    return N, h, c


def main():
    if '-d' in sys.argv[1:]:
        debug = True
    if '-L' in sys.argv[1:]:
        N, h, c = Liouville_parametrization_inp()
    else:
        h = sp.symbols('h')
        c = sp.symbols('c')
        print("Enter level")
        N = int(input())
    M = Kac_matrix(N, h, c)
    print("Kac-Shapovalov form on level ", N, ":\n", M)
    if ('-L' in sys.argv[1:]) and ('-s' in sys.argv[1:]):
        print("Singular vectors:\n", M.nullspace())

if __name__ == "__main__":
    main()
