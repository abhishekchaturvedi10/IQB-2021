'''
    Name : ABHISHEK CHATURVEDI
    Roll No.: 2019401

    ** The written answers are in the write-up given with this file **
'''




'''
    This function recursively finds the optimal alignments for global alignment.
    It takes the DP matrix, no of columns, no of rows, the srings of alignment(x, y) and the scoring values
    (match, mismatch, gap) as parameters.
'''

def generateOptimalGloabalAlignments(DP, j, i, x, y, match, mismatch, gap):

    '''
        I f we reach the 0,0 cell then we know that we have found and optimal alignment and we add that alignment to he list
        of optimal alignments 
    '''
    if i==0 and j==0:
        OptimalAlignments.append([x, y])
        return

        '''
            Check if the diagonally left cell's coordinates corresponding characters of the sequences match and if they do 
            then add the characters to strings of the alignment and move to the diagonally left cell given that the current 
            cell's score is equal to the match score + diagonally left cell score 
        '''
    if i>0 and j>0 and DP[i][j]==DP[i-1][j-1]+match and sequence1[j - 1] == sequence2[i - 1]:
        generateOptimalGloabalAlignments(DP, j-1, i-1, sequence1[j-1] + x, sequence2[i-1] + y, match, mismatch, gap)


        '''
            If there is no match and mismatch then check for the gap score from the upper cell and compare the score of the 
            current cell to sum of score of upper cell and the gap score and if they are equal then add a gap in Y string of 
            the alignment and corresponding character from the sequence1 to the string X of the alignment
        '''
    if i>0 and DP[i-1][j]+gap==DP[i][j]:
        generateOptimalGloabalAlignments(DP, j, i-1, '-' + x, sequence2[i-1] + y, match, mismatch, gap)

        '''
            If there is no match and mismatch then check for the gap score from the left cell and compare the score of the 
            current cell to sum of score of left cell and the gap score and if they are equal then add a gap in X string of 
            the alignment and corresponding character from the sequence1 to the string Y of the alignment
        ''' 
    if j>0 and DP[i][j-1]+gap==DP[i][j]:
        generateOptimalGloabalAlignments(DP, j-1, i, sequence1[j-1] + x, '-' + y, match, mismatch, gap)

        '''
            Check if the diagonally left cell's coordinates corresponding characters of the sequences mismatch and if they 
            do then add the characters to strings of the alignment and move to the diagonally left cell given that the 
            current cell's score is equal to the mismatch score + diagonally left cell score 
        '''
    if i>0 and j>0 and DP[i][j]==DP[i-1][j-1]+mismatch and sequence1[j - 1] != sequence2[i - 1]:
        generateOptimalGloabalAlignments(DP, i-1, j-1, sequence1[j-1] + x, sequence2[i-1] + y, match, mismatch, gap)






'''
    This function generates the DP matrix for global sequence alignment and also calculates the optimal score for alignment.
    It returns the DP matrix and the optimal score for the alignment.
    It takes the two sequences parameters.
'''
def global_align(sequence1, sequence2):

    # No of columns of the DP matrix 
    n=len(sequence1)

    # No of rows of the DP matrix 
    m=len(sequence2)

    # Creating a DP matrix of n+1 columns and m+1 rows with initial score 0
    DP = [[0] * (n + 1) for i in range(m + 1)]

    # Filling the first column with the row no. times gap score 
    for i in range(m + 1):
        DP[i][0] = i * gap
    
    # Filling the first row with the column no. times the gap score
    for i in range(n + 1):
        DP[0][i] = i * gap
    
    # Traversing DP matrix row by row
    for i in range(1, m + 1):
        for j in range(1, n + 1):

            # If the corresponding characters of the sequences match then we take the max of upper cell score + gap score, 
            # left cell score + gap score and diagonally left cell score + match score
            if(sequence2[i-1]==sequence1[j-1]):
                DP[i][j] = max(DP[i][j-1] + gap, DP[i-1][j] + gap, DP[i-1][j-1] + match)

            # If the corresponding characters of the sequences do not match then we take the max of upper cell score + gap 
            # score, left cell score + gap score and diagonally left cell score + mismatch score
            else:
                DP[i][j] = max(DP[i][j-1] + gap, DP[i-1][j] + gap, DP[i-1][j-1] + mismatch)


    # Returning the DP matrix and the optimal alignment score
    return [DP, DP[m][n]]







'''
    This function generates the DP matrix for local sequence alignment. It also calculates the optimal score for alignment
    and also returns one optimal alignment which has optimal score.
    It returns the DP matrix, optimal score and the alignment sequences for local alignment.
    It takes the two sequences parameters.
'''
def local_align(sequence1, sequence2):
   
    # No of columns of the DP matrix 
    n=len(sequence1)

    # No of rows of the DP matrix
    m=len(sequence2)

    # Creating a DP matrix of n+1 columns and m+1 rows with initial score 0
    DP = [[0] * (n + 1) for i in range(m + 1)]
    
    # A variable to store the optimal score
    optimal_score = 0

    # A variable to store the location/coordinates of cell with highest score or optimal score
    optimal_location = (0,0)

    # Traversing DP matrix row by row
    for i in range(1, m + 1):
        for j in range(1, n + 1):

            # If the corresponding characters of the sequences match then we take the max of upper cell score + gap score, 
            # left cell score + gap score and diagonally left cell score + match score
        	if(sequence2[i-1]==sequence1[j-1]):
        		DP[i][j] = max(0, DP[i][j-1] + gap, DP[i-1][j] + gap, DP[i-1][j-1] + match)

            # If the corresponding characters of the sequences match then we take the max of upper cell score + gap score, 
            # left cell score + gap score and diagonally left cell score + mismatch score
        	else:
        		DP[i][j] = max(0, DP[i][j-1] + gap, DP[i-1][j] + gap, DP[i-1][j-1] + mismatch)

            # Check if the score of the current cell is greater than the optimal score or not and if it is then upadte the 
            # optimal score and the set the optimal location to the coordinates of the current cell
        	if DP[i][j] < optimal_score:
        		continue
        	else:
        		optimal_score = DP[i][j]
        		optimal_location = (i,j)


    # Variables to store the strings of alignment during the traceback in the DP matrix 
    alignment_X = ""
    alignment_Y = ""

    # Variables to traceback on the optimal path from the optimal location
    j, i = optimal_location


    # We only traverse from the optimal location till we encounter a cell with score 0 or we reach the topmost left corner cell    
    while (i > 0 or j > 0) and DP[j][i] > 0:
         
        current_score = DP[j][i]

        '''
            Check if the diagonally left cell's coordinates corresponding characters of the sequences match and if they do 
            then add the characters to strings of the alignment and move to the diagonally left cell given that the current 
            cell's score is equal to the match score + diagonally left cell score 
        '''
        if i > 0 and j > 0 and sequence1[i - 1] == sequence2[j - 1] and current_score == DP[j - 1][i - 1] + match:
            alignment_Y = sequence2[j - 1] + alignment_Y
            alignment_X = sequence1[i - 1] + alignment_X
            i = i - 1
            j = j - 1

            '''
                Check if the diagonally left cell's coordinates corresponding characters of the sequences mismatch and if they 
                do then add the characters to strings of the alignment and move to the diagonally left cell given that the 
                current cell's score is equal to the mismatch score + diagonally left cell score 
            '''
        elif i > 0 and j > 0 and sequence1[i - 1] != sequence2[j - 1] and current_score == DP[j - 1][i - 1] + mismatch:
            alignment_Y = sequence2[j - 1] + alignment_Y
            alignment_X = sequence1[i - 1] + alignment_X
            i = i - 1
            j = j - 1

            '''
                If there is no match and mismatch then check for the gap score from the upper cell and compare the score of the 
                current cell to sum of score of upper cell and the gap score and if they are equal then add a gap in Y string of 
                the alignment and corresponding character from the sequence1 to the string X of the alignment
            '''
        elif i > 0 and current_score == DP[j][i - 1] + gap:
            alignment_Y = "-" + alignment_Y
            alignment_X = sequence1[i - 1] + alignment_X
            i = i - 1

            '''
                If there is no match and mismatch then check for the gap score from the left cell and compare the score of the 
                current cell to sum of score of left cell and the gap score and if they are equal then add a gap in X string of 
                the alignment and corresponding character from the sequence1 to the string Y of the alignment
            ''' 
        elif j > 0 and current_score == DP[j - 1][i] + gap:
            alignment_Y = sequence2[j - 1] + alignment_Y
            alignment_X = "-" + alignment_X
            j = j - 1

    # Return the DP matrix, the optimal score and the strings of the alignment 
    return [DP, DP[m][n], alignment_X, alignment_Y]








# This functions the prints the matrix which here are the DP tables/matrices for the global and local alignment.
def printDP(matrix, NoOfRows, NoOfColumns):

	print("", end="       ")
	for i in range(NoOfRows):
		print('{:3}'.format(sequence1[i]), end="")
	print()

	for i in range(NoOfColumns + 1):
		if i>0:
			print(sequence2[i-1], end=" ")
		else:
			print("", end="  ")
		for j in range(NoOfRows + 1):
			print('{:3}'.format(matrix[i][j]), end="")
		print()
	print()




   


# The given sequences in question
sequence1 = 'ATCAGAGTA'
sequence2 = 'TTCAGTA'






# //////////////////////////////////////////////// QUESTION 1 BEGINS HERE /////////////////////////////////////////////////////////


# The score of match, mismatch and gap according to the given scoring function in question
match = 2
mismatch = -1
gap = -1





print("/////////////////////////////////////////////  QUESTION 1 ///////////////////////////////////////////////")
print()
print("/////////////////////////////////////////// GLOBAL  ALIGNMENT ///////////////////////////////////////////")
print()

# Calling the global alignment function and storing the score returned by it in a variable
global_alignment_result = global_align(sequence1, sequence2)

print("   Dyanammic Programming matrix")
print()

# Calling the function to print the DP matrix for global sequence alignment
printDP(global_alignment_result[0], len(sequence1), len(sequence2))

print("   Optimal Score = ", end="")

# Printing the optimal score
print(global_alignment_result[1])
print()


# List to store the alignments which have optimal score i.e., alignments which are optimal
OptimalAlignments = []

generateOptimalGloabalAlignments(global_alignment_result[0], len(sequence1), len(sequence2), "", "", match, mismatch, gap)


print("   Optimal Alignments")
print()

'''
    This part prints the optimal alignments in the desired form.
'''
for Alignment in OptimalAlignments:
    
    print("   ", end="")
    print("Score : ", global_alignment_result[1])
    print()

    print("   ", end="")
    print(Alignment[0])
    
    print("   ", end="")
    for i in range(len(Alignment[0])):
        if Alignment[0][i]=='-' or Alignment[1][i]=='-' or Alignment[0][i] != Alignment[1][i]:
            print("", end=" ")
        else:
            print("|", end="")
    print()
    
    print("   ", end="")
    print(Alignment[1])
    print()
    print()


# ////////////////////////////////////////////// QUESTION 1 ENDS HERE ///////////////////////////////////////////////////////////




# ////////////////////////////////////////////// QUESTION 2 BEGINS HERE /////////////////////////////////////////////////////////


# The score of match, mismatch and gap according to the given scoring function in question
match = 2
mismatch = -1
gap = -1



print("/////////////////////////////////////////////  QUESTION 2 //////////////////////////////////////////////")
print()
print("/////////////////////////////////////////// LOCAL  ALIGNMENT ///////////////////////////////////////////")
print()


# Calling the local alignment function and storing the score returned by it in a variable
local_alignment_result = local_align(sequence1, sequence2)

print("   Dyanammic Programming matrix")
print()

#Calling the function to print the DP matrix for global sequence alignment
printDP(local_alignment_result[0], len(sequence1), len(sequence2))


#Printing the optimal score
print("   Optimal Score = ", end="")
print(local_alignment_result[1])
print()



# This part prints the optimal alignments in the desired form.
print("   Optimal Alignment")
print()

print("   ", end="")
print("Score : ", local_alignment_result[1])
print()

print("   ", end="")
print(local_alignment_result[2])

print("   ", end="")
for i in range(len(local_alignment_result[2])):
	if local_alignment_result[2][i]=='-' or local_alignment_result[3][i]=='-' or local_alignment_result[2][i]!=local_alignment_result[3][i]:
		print("", end=" ")
	else:
		print("|", end="")
print()

print("   ", end="")
print(local_alignment_result[3])

print()
print()

# //////////////////////////////////////////////// QUESTION 2 ENDS HERE ///////////////////////////////////////////////////////////







# //////////////////////////////////////////////// QUESTION 4 BEGINS HERE /////////////////////////////////////////////////////////


# Redifining the score of match, mismatch and gap according to the given scoring function in question
match = 2
mismatch = -1
gap = -2



print("/////////////////////////////////////////////  QUESTION 4 ///////////////////////////////////////////////")
print()
print("/////////////////////////////////////////// GLOBAL  ALIGNMENT ///////////////////////////////////////////")
print()


# Calling the global alignment function and storing the score returned by it in a variable
global_alignment_result = global_align(sequence1, sequence2)

print("   Dyanammic Programming matrix")
print()

# Calling the function to print the DP matrix for global sequence alignment
printDP(global_alignment_result[0], len(sequence1), len(sequence2))

print("   Optimal Score = ", end="")

# Printing the optimal score
print(global_alignment_result[1])
print()


# List to store the alignments which have optimal score i.e., alignments which are optimal
OptimalAlignments = []

generateOptimalGloabalAlignments(global_alignment_result[0], len(sequence1), len(sequence2), "", "", match, mismatch, gap)

print("   Optimal Alignments")
print()

'''
    This part prints the optimal alignments in the desired form.
'''
for Alignment in OptimalAlignments:
    
    print("   ", end="")
    print("Score : ", global_alignment_result[1])
    print()

    print("   ", end="")
    print(Alignment[0])
    
    print("   ", end="")
    for i in range(len(Alignment[0])):
        if Alignment[0][i]=='-' or Alignment[1][i]=='-' or Alignment[0][i] != Alignment[1][i]:
            print("", end=" ")
        else:
            print("|", end="")
    print()
    
    print("   ", end="")
    print(Alignment[1])
    print()
    print()





print("/////////////////////////////////////////////  QUESTION 4 //////////////////////////////////////////////")
print()
print("/////////////////////////////////////////// LOCAL  ALIGNMENT ///////////////////////////////////////////")
print()


# Calling the local alignment function and storing the score returned by it in a variable
local_alignment_result = local_align(sequence1, sequence2)

print("   Dyanammic Programming matrix")
print()

#Calling the function to print the DP matrix for global sequence alignment
printDP(local_alignment_result[0], len(sequence1), len(sequence2))


#Printing the optimal score
print("   Optimal Score = ", end="")
print(local_alignment_result[1])
print()



# This part prints the optimal alignments in the desired form.
print("   Optimal Alignment")
print()

print("   ", end="")
print("Score : ", local_alignment_result[1])
print()

print("   ", end="")
print(local_alignment_result[2])

print("   ", end="")
for i in range(len(local_alignment_result[2])):
    if local_alignment_result[2][i]=='-' or local_alignment_result[3][i]=='-' or local_alignment_result[2][i]!=local_alignment_result[3][i]:
        print("", end=" ")
    else:
        print("|", end="")
print()

print("   ", end="")
print(local_alignment_result[3])

# ////////////////////////////////////////////////// QUESTION 4 ENDS HERE //////////////////////////////////////////////////////////