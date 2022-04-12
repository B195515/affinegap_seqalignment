# Perform Smith Waterman Alignment in Python (from first principles)
# By exam number B195515
##############
# Code modified from original by Simon Tomlinson, Bioinformatics Algorithms 2022
# Code adapted from https://github.com/dnase/affine-gap-sequence-alignment
    # for affine gap implementation
# Modified variable names such as i/x to R (rows), and j/y to C (columns) for readability
# Added functions to initialize the three different matrices
# Score calculation is done in one function for simplification
# Reformatted outputs to screen, added functionality to save output to file

def read_fasta(filename):
    seq = ""
    with open(filename, 'r') as filehandle:
        for line in filehandle:
            #print(line)
            if(line[0]==">"):
                name = line.replace(">","").rstrip("\n")
                continue
            # multiseq file will have concatenated sequences
            seq = seq+line
        import re
        pattern = re.compile(r'\s+')
        seq = re.sub(pattern, '', seq)
        seq = seq.replace('\"', '')
        return name, seq

#return match or mismatch score
def _match(seq1, seq2, R, C):
    if seq2[R-1] == seq1[C-1]: return seqmatch
    else: return seqmismatch

#initializers for matrices
def _init_x(R, C):
    # x(R,0) = h + (g*R)
    if R > 0 and C == 0: return MIN
    else: 
        if C > 0: return S + (E * C)
        else: return 0

def _init_y(R, C):
    # y(0,C) = h + (g*C)
    if C > 0 and R == 0: return MIN
    else:
        if R > 0: return S + (E * R)
        else: return 0

def _init_m(R, C):
    # M(0,0) = 0
    # other cells in top row and leftmost column = -inf
    if C == 0 and R == 0: return 0
    else:
        if C == 0 or R == 0: return MIN
        else: return 0

def createBuild_matrix(seq1, seq2):
    # define rows/cols +1 here is faster than in list comprehension
    rows = len(seq2)+1
    cols = len(seq1)+1
    # Create the three matrices and initial scoring
    X = [[_init_x(R, C) for C in range(cols)] for R in range(rows)]
    Y = [[_init_y(R, C) for C in range(cols)] for R in range(rows)]
    M = [[_init_m(R, C) for C in range(cols)] for R in range(rows)]
    # Then find the max score
    for C in range(1, cols):
        for R in range(1, rows):
            # Matrix X: best score given that x[R] aligns to a gap
            X[R][C] = max((S + E + M[R][C-1]), (E + X[R][C-1]), (S + E + Y[R][C-1]))
            # Matrix Y: best score given that y[C] aligns to a gap
            Y[R][C] = max((S + E + M[R-1][C]), (S + E + X[R-1][C]), (E + Y[R-1][C]))
            # Matrix M: best score given than x[R] aligns to y[C]
            M[R][C] = max(_match(seq1, seq2, R, C) + M[R-1][C-1], X[R][C], Y[R][C])
    return [X, Y, M]

def traceback(seq1, seq2, X, Y, M):
    # Initialize alignment strings
    topstring = ''
    bottomstring = ''
    midstring = ''
    R = len(seq2)
    C = len(seq1)
    print("Building traceback...")
    while (R>0 or C>0):
        # If match, trace diagonally (-1 on both indexes)
        if (R>0 and C>0 and M[R][C] == M[R-1][C-1] + _match(seq1, seq2, R, C)):
            topstring += seq1[C-1]
            bottomstring += seq2[R-1]
            midstring += "|"
            R -= 1; C -= 1
        # If mismatch, trace horizontally/vertically based on alignment
        # with either matrix X or matrix Y at the same index
        elif (R>0 and M[R][C] == Y[R][C]):
            topstring += '-'
            bottomstring += seq2[R-1]
            midstring += " "
            R -= 1
        elif (C>0 and M[R][C] == X[R][C]):
            topstring += seq1[C-1]
            bottomstring += '-'
            midstring += " "
            C -= 1
    return [topstring, bottomstring, midstring]

def get_max(mymatrix):
    max=mymatrix[0][0]
    mrow = 0
    mcol = 0
    rows = len(mymatrix)
    cols = len(mymatrix[0])

    for R in range(1, rows):
        for C in range(1, cols):
            if mymatrix[R][C]<max:
                # "less than" bcs scores are negative
                max = mymatrix[R][C]
                mrow = R
                mcol = C
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
    for R in range(0, rows):
        # print rows
        print(s2[R], "\t", end="")
        for C in range(0, cols):
            # print score values in best scoring alignment
            print(mymatrix[R][C], "\t", end="")
        print("\n", end="")
    print("##")

def perform_affine_gaps():
    # set score rules
    global seqmatch; global seqmismatch; global S; global E; global MIN
    S = -10.
    E = -0.5
    seqmatch = 1.
    seqmismatch = -4.
    MIN = -float("inf")

    import argparse
    parser = argparse.ArgumentParser(description='Aligning sequences...')
    parser.add_argument('seq1',action="store",help="First sequence")
    parser.add_argument('seq2',action="store",help="Second sequence")
    margs = parser.parse_args()
    #print(margs)

    # input sequences + sequence name in file
    global sequence1; global sequence2
    n1, sequence1 = read_fasta(margs.seq1)
    n2, sequence2 = read_fasta(margs.seq2)

    # Print to screen and save to file
    # import sys
    # filepath = f"affinegap_{n1}_{n2}.txt"
    # print(f"Output saved to {filepath}")
    # print("\nEnd.")
    # sys.stdout = open(filepath, 'w')
    print("##")
    print(f"Name: {n1}\nSequence1: {sequence1}\n##")
    print(f"Name: {n2}\nSequence2: {sequence2}\n##")

    [X, Y, M] = createBuild_matrix(sequence1, sequence2)
    [topstring, bottomstring, midstring] = traceback(sequence1, sequence2, X, Y, M)
    print_matrix(M)
    get_max(M)

    print("##\nAlignment with affine gap penalties")
    print(f"Rules:\nGap opening ({S}), Gap extension ({E}), Match ({seqmatch}), Mismatch ({seqmismatch})")
    print(topstring[::-1])
    print(midstring[::-1])
    print(bottomstring[::-1])
    print("\n")
    # sys.stdout.close()

## Run program
perform_affine_gaps()