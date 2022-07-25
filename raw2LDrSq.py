import numpy as np
# import pandas as pd
# matrix = Trial_opttxt[1:,1:]


def rm_hetero(matrix):
    "This function takes a matrix of two rows of SNPs\
    and return the same matrix without heterozygous SNPs\
    or error SNPs like NN"
    errors = []
    for i in range(len(matrix[1])):
        if matrix[0, i] == "NN" or matrix[1, i]=="NN" or matrix[0,i][0]!=matrix[0,i][1] or matrix[1,i][0]!=matrix[1,i][1]:
            errors.append(i)
    matrix = np.delete(matrix, errors, 1)
    return matrix


# ex = matrix[0:2,]
# ex[1,1] = "NN"
# ex[1,9] = "NN"
# ex[0,8] = "NN"
# ex[0,0] = "CG"
# print(ex)
# print(rm_hetero(ex))

# unit = rm_hetero(ex)

def unique_haplotypes_counts(matrix):
    "This function transform a 2-row matrix into unique\
    haplotype dictionary of related frequencies"
    unique_list_1 = []
    for x in matrix[0,]:
        if x not in unique_list_1:
            unique_list_1.append(x)
    unique_list_2 = []
    for x in matrix[1,]:
        if x not in unique_list_2:
            unique_list_2.append(x)
    lst = unique_list_1+unique_list_2
    haplotypes = {}
    haplotypes[lst[0]+lst[2]] = 0
    haplotypes[lst[0]+lst[3]] = 0
    haplotypes[lst[1]+lst[2]] = 0
    haplotypes[lst[1]+lst[3]] = 0
    new_lst = []
    for i in range(len(matrix[0])):
        new_lst.append(matrix[0,i]+matrix[1,i])
    for i in haplotypes.keys():
        haplotypes[i] = new_lst.count(i)
    return haplotypes

# x = unique_haplotypes_counts(unit)
# print(sum(x.values()))
def LD_r_sq(freq_dict):
    total = sum(freq_dict.values())
    haps = list(freq_dict.values())
    D = ((haps[0]/total)*(haps[3]/total))-((haps[1]/total)*(haps[2]/total))
    r_sq = (D**2)/((haps[0]+haps[2])*(haps[1]+haps[3])*(haps[0]+haps[1])*(haps[2]+haps[3]))
    return r_sq

#print(LD_r_sq(x))


def LD_r_sq_matrix(SNP_matrix):
    matrix_length = len(SNP_matrix)
    print(matrix_length)
    result = np.zeros((matrix_length,matrix_length))
    for i in range(matrix_length):
        print(i)
        for n in range(i):
            row_idx = np.array([i,n])
            pair = rm_hetero(SNP_matrix[row_idx,:])
            haplo = unique_haplotypes_counts(pair)
            result[i,n] = LD_r_sq(haplo)

    return result

  
##For testing
data = np.genfromtxt(fname = "D:\Sarhan_Research_&_Work\Publications\Feenugreek_SNPs\_trial\soda.txt", dtype= 'str')
print(LD_r_sq_matrix(data[1:,1:]))

