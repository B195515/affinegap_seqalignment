# Perform Smith Waterman Alignment in Python (from first principles)
# By exam number B195515
##############
# Code modified from original by Simon Tomlinson, Bioinformatics Algorithms 2022
# Code adapted from https://github.com/dnase/affine-gap-sequence-alignment
    # for affine gap implementation

def read_fasta(filename):
    seq = ""
    with open(filename, 'r') as filehandle:
        for line in filehandle:
            #print(line)
            if(line[0]==">"):
                name = line.replace(">","").rstrip("\n")
                continue
            seq = seq+line
        import re
        pattern = re.compile(r'\s+')
        seq = re.sub(pattern, '', seq)
        seq = seq.replace('\"', '')
        return name, seq

#return match or mismatch score
def _match(seq1, seq2, i, j):
    if seq2[i-1] == seq1[j-1]:
        return seqmatch
    else:
        return seqmismatch

#initializers for matrices
def _init_x(i, j):
    if i > 0 and j == 0:
        return MIN
    else:
        if j > 0:
            return S + (E * j)
        else:
            return 0

def _init_y(i, j):
    # y(0,j) = h + (g*j)
    if j > 0 and i == 0:
        return MIN
    else:
        if i > 0:
            return S + (E * i)
        else:
            return 0

def _init_m(i, j):
    # M(0,0) = 0
    if j == 0 and i == 0:
        return 0
    else:
        if j == 0 or i == 0:
            # other cells in top row and leftmost column = -inf
            return MIN
        else:
            return 0

def distance_matrix(seq1, seq2):
    rows = len(seq2) + 1
    cols = len(seq1) + 1
    # Create the three matrices and score with affine gap penalties
    X = [[_init_x(i, j) for j in range(0, cols)] for i in range(0, rows)]
    Y = [[_init_y(i, j) for j in range(0, cols)] for i in range(0, rows)]
    M = [[_init_m(i, j) for j in range(0, cols)] for i in range(0, rows)]
    for j in range(1, cols):
        for i in range(1, rows):
            # Matrix X: best score given that x[i] aligns to a gap
            X[i][j] = max((S + E + M[i][j-1]), (E + X[i][j-1]), (S + E + Y[i][j-1]))
            # Matrix Y: best score given that y[j] aligns to a gap
            Y[i][j] = max((S + E + M[i-1][j]), (S + E + X[i-1][j]), (E + Y[i-1][j]))
            # Matrix M: best score given than x[i] aligns to y[j]
            M[i][j] = max(_match(seq1, seq2, i, j) + M[i-1][j-1], X[i][j], Y[i][j])
    return [X, Y, M]

def backtrace(seq1, seq2, X, Y, M):
    # Initialize alignment strings
    topstring = ''
    bottomstring = ''
    midstring = ''
    i = len(seq2)
    j = len(seq1)
    print("Building traceback...")
    while (i>0 or j>0):
        # If match, trace diagonally (-1 on both indexes)
        if (i>0 and j>0 and M[i][j] == M[i-1][j-1] + _match(seq1, seq2, i, j)):
            topstring += seq1[j-1]
            bottomstring += seq2[i-1]
            midstring += "|"
            i -= 1; j -= 1
        # If mismatch, trace horizontally/vertically based on alignment
            # with either matrix X or matrix Y at the same index
        elif (i>0 and M[i][j] == Y[i][j]):
            topstring += '-'
            bottomstring += seq2[i-1]
            midstring += " "
            i -= 1
        elif (j>0 and M[i][j] == X[i][j]):
            topstring += seq1[j-1]
            bottomstring += '-'
            midstring += " "
            j -= 1
    return [topstring, bottomstring, midstring]

def get_max(mymatrix):
    max=mymatrix[0][0]
    mrow = 0
    mcol = 0
    rows = len(mymatrix)
    cols = len(mymatrix[0])

    for i in range(1, rows):
        for j in range(1, cols):
            if mymatrix[i][j]<max:
                # "less than" bcs scores are negative
                max = mymatrix[i][j]
                mrow = i
                mcol = j
    print(f"Max score: {max} at r:{mrow} c:{mcol}")
    # get index of lowest scoring point
    return [mrow,mcol]

def print_matrix(mymatrix):
    rows = len(mymatrix)
    cols = len(mymatrix[0])
    s1 = "  "+sequence1
    s2 = " "+sequence2
    print(f"Dimensions: Rows= {rows}, Columns= {cols}")
    for base in s1:
        # print columns
        print(base, "\t", end="")
    print("\n", end="")
    for i in range(0, rows):
        # print rows
        print(s2[i], "\t", end="")
        for j in range(0, cols):
            # print score values in best scoring alignment
            print(mymatrix[i][j], "\t", end="")
        print("\n", end="")
    print("##")

#########

import argparse
parser = argparse.ArgumentParser(description='Aligning sequences...')
parser.add_argument('seq1',action="store",help="First sequence")
parser.add_argument('seq2',action="store",help="Second sequence")
margs = parser.parse_args()
#print(margs)

# set vars as global if in function
S = -10.
E = -0.5
seqmatch = 1.
seqmismatch = -4.
MIN = -float("inf")

# Read sequences from cmd line args
n1, sequence1 = read_fasta(margs.seq1)
n2, sequence2 = read_fasta(margs.seq2)
print("##")
print(f"Name: {n1}\nSequence1: {sequence1}\n##")
print(f"Name: {n2}\nSequence2: {sequence2}\n##")

[X, Y, M] = distance_matrix(sequence1, sequence2)
[topstring, bottomstring, midstring] = backtrace(sequence1, sequence2, X, Y, M)
print_matrix(M)
get_max(M)

print("##\nAlignment with affine gap penalties")
print(topstring[::-1])
print(midstring[::-1])
print(bottomstring[::-1])
print("\n")

print(topstring)
print(midstring)
print(bottomstring)
print("\n")