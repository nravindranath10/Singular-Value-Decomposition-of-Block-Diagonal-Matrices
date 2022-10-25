# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 19:45:12 2022

@author: Ravindranath Nemani
"""


import numpy as np
from scipy.linalg import svd
from scipy.stats import unitary_group
from scipy.linalg import block_diag


'''Part 1: Proof that the preliminary algorithm outlined in the 
math.stackexchange solution indeed works'''


def permutation_matrix_to_make_block_diagonal_by_permuting_rows(A):
    A_list = A.tolist()
    zero_positions = []
    nonzero_positions = []
    for i in range(len(A_list)):
        if A_list[i][0] == 0:
            zero_positions.append(i)
        else:
            nonzero_positions.append(i)
    
    new_ordering_of_rows = nonzero_positions + zero_positions
    print(new_ordering_of_rows)
    input("new ordering")
    
    number_of_rows = len(A_list)
    number_of_columns = len(A_list[0])
    print(number_of_rows)
    print(number_of_columns)
    P = np.zeros((number_of_rows, number_of_columns))
    
    for j in range(len(P)):
        print(j)
        print(new_ordering_of_rows[j])
        P[j][new_ordering_of_rows[j]] = 1

    return P


def return_correct_sigma_matrix(singular, number_of_columns_of_U, number_of_rows_of_V_transpose):
    print(singular)
    Sigma = np.diag(singular)
    Sigma_list = Sigma.tolist()
    print(Sigma_list)
    Sigma_final_list = []
    
    if number_of_rows_of_V_transpose > number_of_columns_of_U:
        print("case1")
        diff = number_of_rows_of_V_transpose - number_of_columns_of_U    
        for k in range(diff):
            for l in Sigma_list:
                l1 = l
                l1.append(0.0)
                Sigma_final_list.append(l1)
    elif number_of_columns_of_U > number_of_rows_of_V_transpose:
        print("case2")
        diff = number_of_columns_of_U - number_of_rows_of_V_transpose   
        for k in range(diff):
            new_zero_row = np.zeros((number_of_rows_of_V_transpose)).tolist()
            Sigma_list.append(new_zero_row)
        Sigma_final_list = Sigma_list
    else:
        print("case3")
        Sigma_final_list = Sigma    
    
    print(Sigma_final_list)
    input("final")
    return Sigma_final_list

M = np.matrix([[1,1,2,0,0], [2,1,3,0,0], [0,0,0,1,2], [0,0,0,2,1]])
#M = M.T


U, singular, V_transpose = svd(M)


number_of_columns_of_U = len(U[0])
number_of_rows_of_V_transpose = len(V_transpose)


Sigma_final_list = return_correct_sigma_matrix(singular, number_of_columns_of_U, number_of_rows_of_V_transpose)
print(Sigma_final_list)
Sigma = np.asarray(Sigma_final_list, dtype=object)
Sigma_original = Sigma
input("Sigma")


print(U.T)
input("U.T")
P1_transpose = permutation_matrix_to_make_block_diagonal_by_permuting_rows(U.T)
P1 = P1_transpose.T
print(P1)
input("P1")


P1_inv = np.linalg.inv(P1)
print(P1_inv)
input("P1_inv")


P2 = permutation_matrix_to_make_block_diagonal_by_permuting_rows(V_transpose)
print(P2)
input("P2")


P2_inv = np.linalg.inv(P2)
print(P2_inv)
input("P2_inv")


print(Sigma_original)
input("Sigma_original")


Sigma_permuted = np.dot(P1_inv, Sigma_original)
Sigma_permuted = np.dot(Sigma_permuted, P2_inv)
print(Sigma_permuted)
input("Sigma_permuted")


#Singular vlaues of individual blocks

A = np.matrix([[1,1,2], [2,1,3]])
#A = A.T


UA, singularA, VA_transpose = svd(A)
print(singularA)
input("Singular values of A block")


B = np.matrix([[1,2], [2,1]])
#B = B.T


UB, singularB, VB_transpose = svd(B)
print(singularB)
input("Singular values of B block")


'''Part 2: Verification that after pre-multiplying M by a unitary block
diagonal matrix U and post-multiplying M with a unitary block diagonal
matrix V, the singular values of UMV remain the same as those of M'''


M = np.matrix([[1,1,2,0,0], [2,1,3,0,0], [0,0,0,1,2], [0,0,0,2,1]])

A = [[1,1,2], [2,1,3]]
B = [[1,2], [2,1]]


U, singular_of_M, V_transpose = svd(M)


UA = unitary_group.rvs(dim=2)
UB = unitary_group.rvs(dim=2)

U = block_diag(UA, UB)

VA = unitary_group.rvs(dim=3)
VB = unitary_group.rvs(dim=2)

V = block_diag(VA, VB)

M2 = np.dot(U, M)
M2 = np.dot(M2, V)

U, singular_of_UMV, V_transpose = svd(M2)

input("next we see comparison of singular values of M and UMV...")
print(singular_of_M)

print(singular_of_UMV)
