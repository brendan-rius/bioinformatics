import pickle
from frequent_words import frequent_words_mismatch
from skew import find_minimum_skew

filename = 'salmonella_genome.txt'
with open(filename, 'r') as file:
    genome = file.read()

k = 9
d = 1
min_skews = find_minimum_skew(genome)
print(' '.join([str(v) for v in min_skews]))
skew = min_skews[0]

window1 = genome[skew:skew + 500]
res1 = frequent_words_mismatch(window1, k, d, True)
with open('w1_' + filename, 'wb') as file1:
    pickle.dump((window1, res1), file1)
print(res1)

window2 = genome[skew - 500:skew]
res2 = frequent_words_mismatch(window2, k, d, True)
with open('w2_' + filename, 'wb') as file2:
    pickle.dump((window1, res2), file2)
print(res2)

window3 = genome[skew - 250:skew + 250]
res3 = frequent_words_mismatch(window3, k, d, True)
with open('w3_' + filename, 'wb') as file3:
    pickle.dump((window1, res3), file3)
print(res3)
