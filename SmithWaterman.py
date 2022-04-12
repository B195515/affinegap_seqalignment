#Simon Tomlinson Bioinformatics Algorithms
#Perform Smith Waterman Alignment in Python (from first principles)
#contains rows lists each of length cols initially set to 0
#index as my_matrix[1][2] my_matrix[R][C]

def create_matrix(rows, cols):
    my_matrix = [[0 for col in range(cols+1)] for row in range(rows+1)]
    return my_matrix

#x is row index, y is column index
# follows[r][c]
def calc_score(matrix, x, y):
    #print("seq1:",sequence1[y- 1]," seq2: "+sequence2[x - 1],"x:",x," y:",y)
    #check if nucleotide matches
     sc = seqmatch if sequence1[y-1]==sequence2[x-1] else seqmismatch
     base_score = matrix[x-1][y-1] + sc
     insert_score = matrix[x-1][y] + seqgap
     delete_score = matrix[x][y-1] + seqgap
     # get the maximum score for the position
     v = max(0, base_score, insert_score, delete_score)
     return v

#makes a single traceback step
def traceback(mymatrix,maxv):
    x = maxv[0]
    y = maxv[-1]
    val = mymatrix[x][y]
    # Do a -2 because its tracing back two steps
    sc = seqmatch if sequence1[y-2]==sequence2[x-2] else seqmismatch
    base_score = mymatrix[x-1][y-1] + sc
    if base_score==val:
        return [x-1,y-1]
    insert_score = mymatrix[x-1][y] + seqgap
    if insert_score==val:
        return [x-1,y]
    else:
        return [x,y-1]

#builds the initial scoring matrix used for traceback
def build_matrix(mymatrix):
    # cols at 0 because the calculation starts from "nothing"?
    rows = len(mymatrix)
    cols = len(mymatrix[0])
    for i in range(1, rows):
         for j in range(1, cols):
          mymatrix[i][j] = calc_score(mymatrix, i, j)
    return mymatrix

#gets the max value from the built matrix
def get_max(mymatrix):
     max=mymatrix[0][0]
     mrow = 0
     mcol = 0
     rows = len(mymatrix)
     cols = len(mymatrix[0])

     for i in range(1, rows):
         for j in range(1, cols):
             if mymatrix[i][j]>max:
                 max = mymatrix[i][j]
                 mrow = i
                 mcol = j
     print(f"max score:{max}")
     # why do this return?
     return [mrow,mcol]

#print out the best scoring path from the SW matrix
def print_matrix(mymatrix):
    rows = len(mymatrix)
    cols = len(mymatrix[0])
    s1 = "  "+sequence1
    s2 = " "+sequence2

    print("Dimensions: r= %2d , c= %2d" % (rows, cols))

    for a in s1:
        print(a, end ="")
        print(" \t", end ="")
    print("\n",end="")

    for i in range(0, rows):
        print(s2[i], end ="")
        print(" \t", end ="")
        for j in range(0, cols):
           print("%02d\t" % (mymatrix[i][j]),end="")
        print("\n",end="")

#print out the traceback of the best scoring alignment
def print_traceback(mymatrix):
    #this will print as expected with internal gaps
    print("Building traceback...")
    maxv = get_max(mymatrix)
    print(maxv)

    #traverse the matrix to find the traceback elements
    #if more than one path just pick one
    topstring = ""
    midstring = ""
    bottomstring = ""
    #pad the sequences so indexes into the sequences match the matrix indexes
    asequence1 = "#" + sequence1
    asequence2 = "#" + sequence2

    # add first element and go to previous
    topstring += asequence1[maxv[-1]]
    bottomstring += asequence2[maxv[0]]
    if asequence1[maxv[-1]] == asequence2[maxv[0]]:
        midstring += "|"
    else:
        midstring += " "

    #Store the old position so we can track if it is an insertion of deletion
    #code assumes highest score cannot be a gap!
    old_maxv=maxv

    #add the rest of the elements
    search=True
    while(search):
        maxv=traceback(mymatrix,maxv)
        if(mymatrix[maxv[0]][maxv[-1]]==0):
            search=False
            continue

        #deal with single gaps or matches
        if(old_maxv[-1] == maxv[-1]):
            topstring += "-"
            bottomstring += asequence2[maxv[0]]
        elif(old_maxv[0] == maxv[0]):
            # insertion or deletion
            topstring += asequence1[maxv[-1]]
            bottomstring += "-"
        else:
            # add the next element and go to previous
            topstring += asequence1[maxv[-1]]
            bottomstring += asequence2[maxv[0]]

        if (asequence1[maxv[-1]] == asequence2[maxv[0]]) \
        & (old_maxv[-1] != maxv[-1]) \
        & (old_maxv[0] != maxv[0]):
            midstring += "|"
        else:
            midstring += " "
        old_maxv = maxv

    print(topstring[::-1])
    print(midstring[::-1])
    print(bottomstring[::-1])

# Use Biopython for real projects that read fasta!!!
def read_fasta_filename2(filename):
    seq = ""
    with open(filename, 'r') as filehandle:
        for line in filehandle:
            if(line[0]== ">"):
                continue
            seq = seq+line
            #OK but better to use a regular expression re.search(pattern, line):
            #eg re.search("^>",line)
            #We need to really search for ">" as the first character
            #there may be whitespaces!!
        #delete any whitespaces & quotes
        import re
        pattern = re.compile(r'\s+')
        seq = re.sub(pattern, '', seq)
        seq = seq.replace('\"', '')
        return seq

#build the SW alignment...
def perform_smith_waterman():
    #values for weights
    global seqmatch
    global seqmismatch
    global seqgap
    global sequence1
    global sequence2

    #note these are not the exact weights used in the original SW paper
    # These can be cmd line parameters but is default here for now
    seqmatch = 1
    seqmismatch = -1
    seqgap = -1

    import argparse
    parser = argparse.ArgumentParser(description='Aligning sequences...')
    parser.add_argument('seq1',action="store",help="First sequence")
    parser.add_argument('seq2',action="store",help="Second sequence")
    # parser.add_argument('anum',action="store",help="A number",type =int)
    margs = parser.parse_args()
    print(margs)

    #input sequences
    sequence1 = read_fasta_filename2(margs.seq1)
    sequence2 = read_fasta_filename2(margs.seq2)

    print("Sequence1: "+sequence1)
    print("Sequence2: "+sequence2)

    mymatrix = create_matrix(len(sequence2), len(sequence1))
    mymatrix = build_matrix(mymatrix)
    print_matrix(mymatrix)
    print_traceback(mymatrix)

##this calls the SW algorithm when the script loads
perform_smith_waterman()
