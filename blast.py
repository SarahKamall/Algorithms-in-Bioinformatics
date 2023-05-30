import re
import warnings
from Bio import BiopythonDeprecationWarning
warnings.simplefilter('ignore', BiopythonDeprecationWarning)
from Bio.SubsMat.MatrixInfo import blosum62

# take query sequence
print("enter the query sequence")
query = input()

# word length
word_len = 3


# step 1 remove low complexity and repeated sequence
def remove_repeats():
    count = 0
    new_seq = ""
    seq1 = query[0]
    seq2 = query[1]
    var1 = ''
    var2 = ''
    # QCEgcgcgcGHI
    for i in range(2, len(query)):
        if query[i] == seq1 and query[i + 1] == seq2:
            count += 1
            var1 = seq1
            var2 = seq2
        seq1 = seq2
        seq2 = query[i]
    if count >= 3:
        for i in range(len(query)):
            if query[i] != var1 and query[i] != var2:
                new_seq += query[i]
        return new_seq
    else:
        return query


print("sequence without repeats " + remove_repeats())

# open file
databasefile = open('database_sequence.txt')

# step 2 splitting sequence into words
new_query = remove_repeats()
word = []
for i in range(0, len(new_query)):
    if i + word_len - 1 < len(new_query):
        word.append(new_query[i:i + word_len])
print("words ", end="")
print(word)


# calculate score
def cal_score(seq1, y, matrix):
    score = 0
    for i in range(len(seq1)):
        seq = (seq1[i], y[i])
        if seq not in blosum62:
            reverse_pair = tuple(reversed(seq))
            score += blosum62[reverse_pair]
        else:
            score += blosum62[seq]
    return score

# calculate score
def calculate_score(x, y, matrix):
    score = 0
    seq = (x, y)
    if seq not in blosum62:
        reverse_pair = tuple(reversed(seq))
        score = blosum62[reverse_pair]
    else:
        score = blosum62[seq]
    return score

# step 3 neighbourhood words
def neighbourwords():
    AminoAcid = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
                 'B', 'Z', 'X']
    str2 = ""
    score = 0

    compination_dic = {}
    for i in range(len(word)):
        for j in range(0, 2):
            for l in range(len(AminoAcid)):
                if j == 0:
                    str2 = word[i]
                    str2 = AminoAcid[l] + str2[1] + str2[2]
                    score = cal_score(word[i], str2, blosum62)
                    compination_dic[str2] = score
                elif j == 1:
                    str2 = word[i]
                    str2 = str2[0] + AminoAcid[l] + str2[2]
                    score = cal_score(word[i], str2, blosum62)
                    compination_dic[str2] = score
                elif j == 2:
                    str2 = word[i]
                    str2 = str2[0] + str2[1] + AminoAcid[l]
                    score = cal_score(word[i], str2, blosum62)
                    compination_dic[str2] = score
    return compination_dic


print("n" , neighbourwords())


# filter seeds for specific threshold
def threshold(t=11):
    seeds = {}
    neighwords = neighbourwords()
    key = list(neighwords.keys())
    value = list(neighwords.values())
    for i in range(len(key)):
        if value[i] >= t:
            seeds[key[i]] = value[i]
    return seeds


print("threshold " ,threshold())


def exactmatch():
    databaseseq = databasefile
    word_hit = {}
    seeds = threshold()
    key = list(seeds.keys())
    for seed in key:
        list_Seq = []
        for seq in databaseseq:
            if seed in seq:
                list_Seq.append(seq[:-2])
        word_hit[seed] = list_Seq
        databaseseq.seek(0)
        if word_hit[seed] == []:
            del word_hit[seed]
    return word_hit


print("extend ", exactmatch())


# extend seeds
def extention(dic):
    database = databasefile.readlines()
    global query
    for seed in dic:
        scores = {}
        for seq in dic[seed]:
            final_query = ''
            score = 0
            database_index = seq.find(seed)
            F_index_database = database_index
            L_index_database = database_index+2
            pattern1 = '.' + seed[1] + seed[2]
            pattern2 = seed[0] + '.' + seed[2]
            pattern3 = seed[0] + seed[1] + '.'
            pattern = pattern1 + '|' + pattern2 + '|' + pattern3
            try:
                F_index_query = (re.search(pattern, query)).start()
            except:
                continue
            L_index_query = F_index_query + 2
            final_query = seed[0] + seed[1] + seed[2]
            score = cal_score(seed, final_query, blosum62)
            extend_from_right = L_index_query + 1
            extend_from_left = F_index_query - 1
            Move_to_right = L_index_database + 1
            Move_to_left = F_index_database - 1
            while Move_to_right< len(seq) and extend_from_right< len(query) and Move_to_left >= 0 and extend_from_left >= 0:
                score = score + calculate_score(query[extend_from_right], seq[Move_to_right], blosum62)
                final_query= final_query + query[extend_from_right]
                extend_from_right = extend_from_right + 1
                Move_to_right = Move_to_right + 1

                score = score + calculate_score(query[extend_from_left], seq[Move_to_left], blosum62)
                final_query = query[extend_from_left] + final_query
                extend_from_left = extend_from_left - 1
                Move_to_left = Move_to_left - 1

            while Move_to_right< len(seq) and extend_from_right< len(query):
                score = score + calculate_score(query[extend_from_right],seq[Move_to_right],blosum62)
                final_query= final_query + query[extend_from_right]
                extend_from_right = extend_from_right + 1
                Move_to_right = Move_to_right + 1
            while Move_to_left >= 0 and extend_from_left >=0 :
                score = score + calculate_score(query[extend_from_left], seq[Move_to_left], blosum62)
                final_query = query[extend_from_left] + final_query
                extend_from_left = extend_from_left - 1
                Move_to_left = Move_to_left - 1
            scores[score] = final_query
        try:
            maxKey = max(scores, key=scores.get)
        except:
            continue
        print(scores[maxKey])
        right = extension_threshold_right(F_index_database,seq,seed,scores[maxKey])
        left = extension_threshold_left(L_index_database,seq,seed,scores[maxKey])
        final = ''
        for k in range(left, right, 1):
            final = final + (scores[maxKey])[k]
        print(seed ," ", final)


def extension_threshold_right(F_index_database,seq,seed,final_query):
    array = []
    pattern1 = '.' + seed[1] + seed[2]
    pattern2 = seed[0] + '.' + seed[2]
    pattern3 = seed[0] + seed[1] + '.'
    F_index_query = re.search(pattern1 or pattern2 or pattern3, final_query).start()
    # L_index_query = F_index_query + 2
    accum = 0
    for i in range(F_index_query,len(final_query)):
        print(calculate_score(final_query[i],seq[F_index_database],blosum62))
        accum = accum + calculate_score(final_query[i],seq[F_index_database],blosum62)
        array.append(accum)
        F_index_database = F_index_database + 1
    max = 0
    index = 0
    for i in range(len(array)):
        if array[i] >= max:
            max = array[i]
            index = i
    return index + F_index_query + 1

def extension_threshold_left(L_index_database, seq, seed, final_query):
    array = []
    pattern1 = '.' + seed[1] + seed[2]
    pattern2 = seed[0] + '.' + seed[2]
    pattern3 = seed[0] + seed[1] + '.'
    F_index_query = re.search(pattern1 or pattern2 or pattern3, final_query).start()
    L_index_query = F_index_query + 2
    accum = 0
    for i in range(L_index_query, 0, -1):
        accum = accum + calculate_score(final_query[i], seq[L_index_database], blosum62)
        array.append(accum)
        L_index_database = L_index_database - 1
    max = 0
    index = 0
    for i in range(len(array)):
        if array[i] >= max:
            max = array[i]
            index = i
    return L_index_query - index


extention(exactmatch())
# close file
databasefile.close()
