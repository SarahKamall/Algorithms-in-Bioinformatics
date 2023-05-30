firstString = "park"
secondString = "spake"
matrix = [[0 for j in range(len(firstString) + 1)] for i in range(len(secondString) + 1)]
notations = [['' for j in range(len(firstString) + 1)] for i in range(len(secondString) + 1)]

for i in range(len(firstString) + 1):
    matrix[0][i] = i
for i in range(len(secondString) + 1):
    matrix[i][0] = i
for i in range(1, len(secondString) + 1):
    for j in range(1, len(firstString) + 1):
        if firstString[j - 1] == secondString[i - 1]:
            matrix[i][j] = matrix[i - 1][j - 1]
            notations[i][j] = "\\"
        else:
            minimum = min(matrix[i - 1][j], matrix[i][j - 1], matrix[i - 1][j - 1])
            matrix[i][j] = minimum + 1
            if(minimum == matrix[i-1][j]):
                notations[i][j] = '|'
            else:
                notations[i][j] = '-'
fs=""
ss=""
i = len(secondString)
j = len(firstString)
while i > 0 and j > 0:
    if(notations[i][j] == "\\"):
        i = i-1
        j = j-1
        fs = fs + firstString[j]
        ss = ss + secondString[i]
    elif notations[i][j] == '|':
        i = i-1
        fs = fs + '-'
        ss = ss + secondString[i]
    elif notations[i][j] == '-':
        j = j - 1
        fs = fs + firstString[j]
        ss = ss + '-'
while i > 0:
    i = i-1
    fs = fs + '-'
    ss = ss + secondString[i]
while j > 0:
    j = j - 1
    fs = fs + firstString[j]
    ss = ss + '-'
print(matrix)
print(notations)
print(ss[::-1])
print(fs[::-1])
