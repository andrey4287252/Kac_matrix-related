import sympy as sp
import sys

debug = False

def Kac_Matrix_element(m, l):
    '''
    Agregates given partitions m and l to one list and computes corresponding element of Kac-Shapovalov matrix
    '''
    L = list(map(int, l)) + list(map(lambda x: -int(x),  m))[::-1]
    return Auxilary_Kac_Matrix_element(L)

def Auxilary_Kac_Matrix_element(L):
    '''
    Returns element of Kac-Shapovalov matrix indexed by partitions agregated in L
    '''
    if debug:
        print("Partitons: ", L)
    h, c = sp.symbols("h, c")
    if len(L) == 0:
       return 1
    if L[len(L) - 1] > 0:
        return 0
    elif L[0] < 0:
        return 0
    elif L[len(L) - 1] == 0:
        return Auxilary_Kac_Matrix_element(L[:len(L) - 1])*h
    elif L[0] == 0:
        return Auxilary_Kac_Matrix_element(L[1:])*h
    else:
        for i in range(len(L) - 1):
            if L[i] >= 0 and L[i + 1] < 0:
                if debug:
                    print("Commutes: ", i, " and ", i + 1)
                    print(Auxilary_Kac_Matrix_element(L[:i] + [L[i] + L[i + 1]] + L[i + 2:]), "appears")
                answ = (L[i] - L[i + 1]) * Auxilary_Kac_Matrix_element(L[:i] + [L[i] + L[i + 1]] + L[i + 2:]) + Auxilary_Kac_Matrix_element(L[:i] + [L[i + 1], L[i]] + L[i + 2:])
                if L[i] == -L[i + 1]:
                    answ += (L[i]**3 - L[i])*c*Auxilary_Kac_Matrix_element(L[:i] + L[i + 2:])/12
                return answ

if '-d' in sys.argv[1:]:
    debug = True

print("Enter first partition")
a = list(map(int, list(map(float, input().split()))))
print("Enter second partition")
b = list(map(int, list(map(float, input().split()))))
print(Kac_Matrix_element(a, b))
