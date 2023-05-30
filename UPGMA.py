import numpy as np
import pandas as pd

def Updating_Sequences(sequences):
    for i in range(len(sequences)):
        if i == dim[0]:
            sequences.remove(f_Sequence)
        elif i == dim[1]:
            sequences.remove(s_Sequence)
    cluster = f_Sequence + s_Sequence
    sequences.append(cluster)
    sequences.sort()
    return sequences


# Taking Inputs
No_Sequences = int(input("Enter the number of Sequences: "))

print("Enter The Names of Sequences: ")
sequences = []
for i in range(No_Sequences):
    sequences.append(input())


Distance_Matrix = []
print("Enter the Distance Matrix: ")
for i in range(No_Sequences):
    Distance_Matrix.append([int(j) for j in input().split()])

# indices = ['A','B','C','D','E']
# Distance_Matrix = [[0,20,60,100,90],[20,0,50,90,80],[60,50,0,40,50],[100,90,40,0,30],[90,80,50,30,0]]

# Transfer Distance Matrix to Frame
Frame_Matrix = pd.DataFrame(data=Distance_Matrix, index=sequences, columns=sequences, dtype=int)

for i in range(len(sequences) - 1):

    # Finding Smallest Number in Distance Matrix
    Smallest = Frame_Matrix.iat[0, 1]
    dim = [0, 1]
    for i in range(len(Frame_Matrix) + 1):
        for j in range(i + 1, len(Frame_Matrix)):
            if Frame_Matrix.iat[i, j] < Smallest:
                Smallest = Frame_Matrix.iat[i, j]
                dim = [i, j]

    f_Sequence = sequences[dim[0]]
    s_Sequence = sequences[dim[1]]

    print("The Smallest Distance is %s" % Smallest)
    print("The Nearest Two Clusters are [%s, %s]" % (f_Sequence, s_Sequence))

    # Updating The Sequences
    sequences = Updating_Sequences(sequences)
    cluster = f_Sequence + s_Sequence

    # Print New Arrange Of Sequences
    print(sequences, "\n")

    # Initialize New Distance Matrix with Zeros
    New_Distance_Matrix = np.zeros((len(sequences), len(sequences)))
    New_Frame_Matrix = pd.DataFrame(data=New_Distance_Matrix, index=sequences, columns=sequences, dtype=int)

    # Filling The Matrix with Calculate Data
    for i in range(len(sequences)):
        for j in range(len(sequences)):
            if i == j:
                New_Frame_Matrix.iloc[i, j] = 0
            elif sequences[i] == cluster:
                FC = Frame_Matrix.at[sequences[j], f_Sequence]
                SC = Frame_Matrix.at[sequences[j], s_Sequence]
                New_Frame_Matrix.iloc[i, j] = (FC + SC) / 2
            elif sequences[j] == cluster:
                FC = Frame_Matrix.at[sequences[i], f_Sequence]
                SC = Frame_Matrix.at[sequences[i], s_Sequence]
                New_Frame_Matrix.iloc[i, j] = (FC + SC) / 2
            else:
                New_Frame_Matrix.iloc[i, j] = Frame_Matrix.at[sequences[i], sequences[j]]

    # Assigning New Frame Matrix To Old Matrix To Continuing The Loop
    Frame_Matrix = New_Frame_Matrix
    # print(New_Frame_Matrix)

    # 0 20 60 100 90
    # 20 0 50 90 80
    # 60 50 0 40 50
    # 100 90 40 0 30
    # 90 80 50 30 0



