import functools

from numpy import matrix, gcd, absolute
from sympy import Matrix, Rational

def equationBalancer(equation_matrix):
    matrix_a = matrix(equation_matrix)

    rref_matrix = (Matrix(matrix_a).rref())[0]
    solution_set = list(rref_matrix[:,len(rref_matrix[:,0])])
    solution_set.append(Rational(1))
    solution_set = absolute(solution_set)

    denoms = []
    for num in solution_set:
        denoms.append(num.q)

    gcd_num = functools.reduce(lambda a,b: a*b // gcd(a,b),denoms)

    solution_set = solution_set * gcd_num
    print(solution_set)

def main():
    equationBalancer('1,0,1,0; 4,0,0,2; 1,2,2,1')
    equationBalancer('1,0,1,0,0,0; 1,0,0,1,0,0; 4,0,0,0,1,0;0,1,1,2,0,2;0,1,0,0,2,0')

if __name__ == '__main__':
    main()


