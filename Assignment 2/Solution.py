from tabulate import tabulate


# Given protein sequence in the question
protein_sequence = 'SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDLLIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF'



# Table/dictionary to store the single letter code for amino acid names
Amino_acid = {}

Amino_acid['A'] = "alanine"
Amino_acid['R'] = "arginine"
Amino_acid['N'] = "asparagine"
Amino_acid['D'] = "aspartic acid"
Amino_acid['C'] = "cysteine"
Amino_acid['E'] = "glutamic acid"
Amino_acid['Q'] = "glutamine"
Amino_acid['G'] = "glycine"
Amino_acid['H'] = "histidine"
Amino_acid['I'] = "isoleucine"
Amino_acid['L'] = "leucine"
Amino_acid['K'] = "lysine"
Amino_acid['M'] = "methionine"
Amino_acid['F'] = "phenylalanine"
Amino_acid['P'] = "proline"
Amino_acid['S'] = "serine"
Amino_acid['T'] = "threonine"
Amino_acid['W'] = "tryptophan"
Amino_acid['Y'] = "tyrosine"
Amino_acid['V'] = "valine"



# Table/dictionary to store the P_alpha values for the amino acids
Pa = {}

Pa["glutamic acid"] = 1.53
Pa["alanine"] = 1.45
Pa["leucine"] = 1.34
Pa["histidine"] = 1.24
Pa["methionine"] = 1.20
Pa["glutamine"] = 1.17
Pa["tryptophan"] = 1.14
Pa["valine"] = 1.14
Pa["phenylalanine"] = 1.12
Pa["lysine"] = 1.07
Pa["isoleucine"] = 1.00
Pa["aspartic acid"] = 0.98
Pa["threonine"] = 0.82
Pa["serine"] = 0.79
Pa["arginine"] = 0.79
Pa["cysteine"] = 0.77
Pa["asparagine"] = 0.73
Pa["tyrosine"] = 0.61
Pa["proline"] = 0.59
Pa["glycine"] = 0.53



# Table/dictionary to store the values of P_beta for the amino acids
Pb = {}

Pb["glutamic acid"] = 0.26
Pb["alanine"] = 0.97
Pb["leucine"] = 1.22
Pb["histidine"] = 0.71
Pb["methionine"] = 1.67
Pb["glutamine"] = 1.23
Pb["tryptophan"] = 1.19
Pb["valine"] = 1.65
Pb["phenylalanine"] = 1.28
Pb["lysine"] = 0.74
Pb["isoleucine"] = 1.60
Pb["aspartic acid"] = 0.80
Pb["threonine"] = 1.20
Pb["serine"] = 0.72
Pb["arginine"] = 0.90
Pb["cysteine"] = 1.30
Pb["asparagine"] = 0.65
Pb["tyrosine"] = 1.29
Pb["proline"] = 0.62
Pb["glycine"] = 0.81



# Boolean array to store the if the amino acid at the some position could be part of alpha helix or not
Helix = [False] * len(protein_sequence)

# Boolean array to store the if the amino acid at the some position could be part of beta strand or not 
Strand = [False] * len(protein_sequence)



# Function to find the nucleation site for aplha helix 
def check_helix():
    
    start = 0

    while start + 6 < len(protein_sequence):
        
        # Variable to count the number of amino acids in the window of 6 contigous amino acids that have P_alpha value >= 1
        numgood = 0
        
        # start to start+6 means window of 6 amino acids for nucleation site for helix
        for i in range(start, start+6):
            if Pa[Amino_acid[protein_sequence[i]]] >= 1:
                numgood += 1

        # If the number of amino acids having P_alpha value >=1 is equal to or greater than 4 then call a function to extend the nucleation site  
        if numgood >= 4:
            extend_helix(start, start+6)

        start = start + 1



# Function to extend the aplha helix nucleation site 
def extend_helix(start, end):


    # Variable to store the average of P_alpha values of window of 4 rightmost amino acids  
    right_extension_avg  = sum([Pa[Amino_acid[x]] for x in protein_sequence[end-3: end+1]]) / 4
    
    # While the right extension average >=1 we keep extending to the right 
    while end < len(protein_sequence) and right_extension_avg >= 1:
        end = end + 1
        
        right_extension_avg  = sum([Pa[Amino_acid[x]] for x in protein_sequence[end-3: end+1]]) / 4


    # Variable to store the average of P_alpha values of window of 4 leftmost amino acids
    left_extension_avg = sum([Pa[Amino_acid[x]] for x in protein_sequence[start-1: start+3]]) / 4

    # While the left extension average >=1 we keep extending to the left
    while  start > 0 and left_extension_avg >= 1:
        start = start - 1
        
        left_extension_avg = sum([Pa[Amino_acid[x]] for x in protein_sequence[start: start+3]]) / 4


    # Marking the amino acids from the start position to the end position i.e., extended region as true which means that they can be alpha helix 
    for i in range(start, end):
        Helix[i] = True



# Function to find the nucleation site for aplha helix 
def check_strand():

    start = 0

    # Variable to count the number of amino acids in the window of 5 contigous amino acids that have P_alpha value >= 1
    while start + 5 < len(protein_sequence):
        
        # Variable to count the number of amino acids in the window of 5 contigous amino acids that have P_alpha value >= 1
        numgood = 0
        
        # start to start+5 means window of 5 amino acids for nucleation site for beta strand
        for i in range(start, start+5):
            if Pb[Amino_acid[protein_sequence[i]]] >= 1:
                numgood += 1
        
        # If the number of amino acids having P_alpha value >=1 is equal to or greater than 3 then call a function to extend the nucleation site
        if numgood >= 3:
            extend_strand(start, start+5)

        start = start + 1



# Function to extend the beta strand nucleation site 
def extend_strand(start, end):


    # Variable to store the average of P_alpha values of window of 4 rightmost amino acids
    right_extension_avg  = sum([Pb[Amino_acid[x]] for x in protein_sequence[end-3: end+1]]) / 4
    
    # While the right extension average >=1 we keep extending to the right
    while end < len(protein_sequence) and right_extension_avg >= 1:
        end = end + 1
        
        right_extension_avg = sum([Pb[Amino_acid[x]] for x in protein_sequence[end-3: end+1]]) / 4


    # Variable to store the average of P_alpha values of window of 4 leftmost amino acids
    left_extension_avg = sum([Pb[Amino_acid[x]] for x in protein_sequence[start-1: start+3]]) / 4

    # While the left extension average >=1 we keep extending to the left
    while start > 0 and left_extension_avg >= 1:
        start = start - 1

        left_extension_avg = sum([Pb[Amino_acid[x]] for x in protein_sequence[start: start+3]]) / 4


    # Marking the amino acids from the start position to the end position i.e., extended region as true which means that they can be beta strand    
    for i in range(start, end):
        Strand[i] = True



# Function to call the functions to find the nucleation sites and extended regions so as to fiind the amino acids which can be part of alpha helix or beta strand
def chou_fasman():
   
    check_helix()
    
    check_strand()



# Function to find and return the secondary structure sequence 
def secondary_structure():

    i = 0
    secondary_sequence = ""

    while i < len(protein_sequence):
    
        # Condition that the amino acid at the ith position cannot be helix or strand 
        if Helix[i] == False and Strand[i] == False:

            j = i

            # Loop forward till the amino acid cannot be helix or strand and add turn to the secondary structure sequence
            while j < len(protein_sequence) and Helix[j] == False and Strand[j] == False:
                secondary_sequence += 'T'
                j += 1

            i = j

        # Condition that the amino acid at the ith position is helix and is not strand    
        elif Helix[i] == True and Strand[i] == False:

            j = i

            # Loop forward till the amino acid is helix and is not strand and add helix to the secondary structure sequence
            while j < len(protein_sequence) and Helix[j] == True and Strand[j] == False:
                secondary_sequence += 'H'
                j += 1

            i = j

        # Conditon to check that the amino acid at the it position is strand and is not helix
        elif Helix[i] == False and Strand[i] == True:

            j = i

            # Loop forward till the amino acid is strand and is not helix and add strand to the secondary structure sequence
            while j < len(protein_sequence) and Helix[j] == False and Strand[j] == True:
                secondary_sequence += 'S'
                j += 1

            i = j

        # Condition to check that the amino acid at the ith position can be helix and strand both
        elif Helix[i] == True and Strand[i] == True:

            j = i
            
            # Variable to store the sum of Pa values of the amino acids
            Helix_sum = 0

            # Varibale to store the sum of Pb values of the amino acids 
            Strand_sum = 0


            # Loop forward till the amino acid can be both helix and strand and add the Pa and Pb values of the amino acid to the Helix_sum and Strand_sum respectively
            while j < len(protein_sequence) and Helix[j] == True and Strand[j] == True:
                Helix_sum += Pa[Amino_acid[protein_sequence[j]]]
                Strand_sum += Pb[Amino_acid[protein_sequence[j]]]
                j += 1

            # If the Helix_sum is greater than or equal to the Strand_sum add helix otherwise strand to the secondary structure sequence for the region 
            # where the amino acid(s) can be both helix and strand 
            for z in range(i, j):
                if Helix_sum >= Strand_sum:
                    secondary_sequence += 'H'
                else:
                    secondary_sequence += 'S'

            i = j
            

    return secondary_sequence



# Function to find and return the differing regions between the secondary structure sequences obtained from chou fasman algorithm and STRIDE webserver
# Function to find the regions of similarity in the sequences
def difference_between_sequences():

    # Array to store the regions which are not same in the two sequences
    regions_differing_in_sec_seqs = []

    differing_regions = ""

    # String to find the similar regions between the sequences
    same_regions = ""

    i = 0

    while i < len(protein_sequence):
    
        # If the amino acids at the ith position are equal then we move one position forward 
        if chou_fasman_sec_seq[i] == STRIDE_sec_seq[i]:
            same_regions += chou_fasman_sec_seq[i]
            i += 1
            differing_regions += '-'
        
        # If the amino acids at the ith position are not equal
        else:
            
            # Array to store the starting and ending positions of the region where the amino acids are not same
            # i+1 is the starting position 
            cur = [i+1]
            j = i + 1
            
            differing_regions += protein_sequence[i]
            same_regions += "-"

            # Loop forward till the amino acids at the position in the two sequences are not equal  
            while j < len(protein_sequence) and chou_fasman_sec_seq[j] != STRIDE_sec_seq[j]:
                differing_regions += protein_sequence[j]
                same_regions += '-'
                j += 1
            
            # j is the ending position 
            cur.append(j)

            # Storing the sequence from chou_fasman secondary sequence from pos i to j
            cur.append(chou_fasman_sec_seq[i:j])

            # Storing the sequence from STRIDE secondary sequence from pos i to j
            cur.append(STRIDE_sec_seq[i:j])
        
            i = j

            # Add the region to array of differing regions 
            regions_differing_in_sec_seqs.append(cur)


    return regions_differing_in_sec_seqs, differing_regions, same_regions


# Calling the chou_fasman function  
chou_fasman()

# Calling secondary_structure and store the returned secondary structure string in chou_fasman_sec_seq 
chou_fasman_sec_seq = secondary_structure()


# Output of Q1

print("Q1")
print()

print("CHOU FASMAN SECONDARY STRUCTURE SEQUENCE")
print()

print(chou_fasman_sec_seq)

print()
print()
print()


# Secondary sequence obtained from STRIDE webserver
STRIDE_sec_seq = "TTTTCCCCCHHHHHHCSSSSSSTTSSSSSSSSTTSSSSSGGGGCCHHHHHCCCHHHHHHHCCGGGCSSSSTTSSSCSSSSSSSTTSSSSSSCCCTTTTCCCCCCCCTTTSSSSSSSSSTTSSSSSSSSSSTTTTCBCCCCTTTTTTTSSC"


# Output of Q2(A) 

print("Q2 (A)")
print()

print("STRIDE WEBSERVER SECONDARY STRUCTURE SEQUENCE")
print()

print(STRIDE_sec_seq)

print()
print()

# Storing the array and strings returned from difference_between_sequences in variables
differing_region_between_sequences, differing_regions_protein, similar_regions = difference_between_sequences()

print("REGIONS DIFFERING IN THE CHOU-FASMAN AND STRIDE-WEBSERVER SEQUENCES")
print()

# Presenting data in a tabular form
print (tabulate(differing_region_between_sequences, headers=["Start position", "End position", "Chou-Fasman subsequence", "Stride-Webserver subsequence"]))

print()
print()

print("REGIONS IN INPUT PROTEIN SEQUENCE DIFFERING IN TWO SEQUENCES (Similar regions marked as -)")
print()

print(differing_regions_protein)

print()
print()

print("REGIONS SAME IN SEQUENCES (Dissimilar regions marked as -)")
print()

print(similar_regions)

print()
print()
print()

# Output of Q2(B)

print("Q2 (B)")
print()

print("STRIDE's method principle makes use of combination of calculated backbone torsional angle and hydrogen bond energy to predict secondary structures. The Chou-Fasman method works on the principle examining of relative frequencies of every amino acid in helices, beta sheets and turns using the known protein structures calculated using crystallography (X-Ray) and using these obtained frequencies the probability values (propensity values) for every amino acid were calculated for the occurrence as helix, beta sheet or turn. The type of secondary structure formed by a sequence is predicted using the probability values or  propensity values. The differences in the sequences are significant since the STRIDE's method of prediction is much more accurate than what we are using for the Chou-Fasman. Our method just calculates and predicts the alpha helix and beta strand and assumes that if none of alpha or beta is predicted then it is turn while the STRIDE's method also calculates and predicts turns, coils, bridges and 310 helixes in addition to alpha helixes and beta strand. The method used by STRIDE is more accurate with respect to crystallography also.")
print()