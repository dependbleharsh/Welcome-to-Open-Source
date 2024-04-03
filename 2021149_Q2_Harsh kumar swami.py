import numpy as np
import pandas as pd

def local_alignment(x, y, match=2, mismatch=-1, gap=-3):
    nx = len(x)
    ny = len(y)
    
    # Initialization of the matrix with zeros.
    F = np.zeros((nx + 1, ny + 1))
    P = np.zeros((nx + 1, ny + 1), dtype=int)

    # Matrix filling.
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            scores = [F[i-1, j-1] + (match if x[i-1] == y[j-1] else mismatch),
                      F[i-1, j] + gap,
                      F[i, j-1] + gap,
                      0]  # Introducing local alignment condition
            F[i, j] = max(scores)
            P[i, j] = scores.index(F[i, j])
            if F[i, j] < 0:
                F[i, j] = 0  # Resetting negative values to 0 for local alignment

    # Finding the maximum score and its position
    max_score = np.max(F)
    max_pos = np.unravel_index(np.argmax(F), F.shape)

    # Traceback from the position of the maximum score to construct the alignment
    alignment_x = ''
    alignment_y = ''
    i, j = max_pos
    while i > 0 and j > 0 and F[i, j] != 0:
        if P[i, j] == 0:
            alignment_x = x[i - 1] + alignment_x
            alignment_y = y[j - 1] + alignment_y
            i -= 1
            j -= 1
        elif P[i, j] == 1:
            alignment_x = x[i - 1] + alignment_x
            alignment_y = '-' + alignment_y
            i -= 1
        elif P[i, j] == 2:
            alignment_x = '-' + alignment_x
            alignment_y = y[j - 1] + alignment_y
            j -= 1
        else:
            break

    return alignment_x, alignment_y, max_score

seq1 = "GATGCGCAG"
seq2 = "GGCAGTA" 

alignment_x, alignment_y, score = local_alignment(seq1, seq2)

print("Alignment:")
print(alignment_x)
print(alignment_y)
print("Score:", score)
