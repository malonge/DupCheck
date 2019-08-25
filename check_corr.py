import sys

gc_total = [0 for i in range(101)]  # Sum of all coverage for the 100 possible GC values
gc_count = [0 for i in range(101)]
with open(sys.argv[1], "r") as f:
    for line in f:
        bed_fields = line.rstrip().split("\t")
        gc_total[int(bed_fields[3])] += int(float(bed_fields[-1]))
        gc_count[int(bed_fields[3])] += 1

gc_avg = []
for i, j in zip(gc_total, gc_count):
    if not j:
        gc_avg.append(0)
    else:
        gc_avg.append(i/j)

print (gc_avg)
import matplotlib.pyplot as plt
plt.bar(range(101), gc_avg)
plt.show()