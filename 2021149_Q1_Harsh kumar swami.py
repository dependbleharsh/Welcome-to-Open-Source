import numpy as np
import pandas as pd

def nw(x, y, match=2, mismatch=-3, gap=-1):
    nx = len(x)
    ny = len(y)
    
    # Initialization of the matrix.
    F = np.zeros((nx + 1, ny + 1))
    F[:, 0] = np.linspace(0, gap * nx, nx + 1) 
    F[0, :] = np.linspace(0, gap * ny, ny + 1)

    # Pointers to trace through all optimal alignments.
    P = np.zeros((nx + 1, ny + 1), dtype=int)
    P[:, 0] = 3
    P[0, :] = 4

    # Matrix filling.
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            scores = [F[i-1, j-1] + (match if x[i-1] == y[j-1] else mismatch),
                      F[i-1, j] + gap,
                      F[i, j-1] + gap]
            F[i, j] = max(scores)
            P[i, j] = scores.index(F[i, j]) + 2
    # Print scoring matrix using pandas DataFrame
    print("Scoring Matrix:")
    df = pd.DataFrame(F, index=['-'] + list(x), columns=['-'] + list(y))
    print(df)

    # Trace through all optimal alignments.
    alignments = []

    def traceback(i, j, rx, ry):
        if i == 0 and j == 0:
            alignments.append((rx[::-1], ry[::-1]))
            return
        if i > 0 and j > 0 and F[i, j] == F[i-1, j-1] + (match if x[i-1] == y[j-1] else mismatch):
            traceback(i - 1, j - 1, rx + x[i - 1], ry + y[j - 1])
        if i > 0 and F[i, j] == F[i-1, j] + gap:
            traceback(i - 1, j, rx + x[i - 1], ry + '-')
        if j > 0 and F[i, j] == F[i, j-1] + gap:
            traceback(i, j - 1, rx + '-', ry + y[j - 1])

    traceback(nx, ny, '', '')
    
    return alignments

seq1 = "GATGCGCAG"
seq2 = "GGCAGTA" 

alignments = nw(seq1, seq2)

print("Number of Alignments:", len(alignments))
print("Alignments and Scores:")
best_alignment = None
best_score = float('-inf')
for i, (alignment1, alignment2) in enumerate(alignments):
    score = sum([2 if a == b else -3 for a, b in zip(alignment1, alignment2)])
    print(f"Alignment {i + 1}:")
    print(alignment1)
    print(alignment2)
    print("Score:", score)
    print()
    if score > best_score:
        best_score = score
        best_alignment = (alignment1, alignment2)

print("Best Alignment (based on max score):")
print(best_alignment[0])
print(best_alignment[1])
print("Score:", best_score)
