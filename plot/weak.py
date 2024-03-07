import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure
import re

figure(figsize=(5, 4), dpi=160)

t = [1, 2, 4, 8, 16, 32, 64, 128]
r = [1000, 4000, 10000, 40000]

tweak = []

for j in t:
  tweak.append([])
  for i in r:
    with open('../data_weak/t' + str(j) + '-' + str(i * j) + '.out', 'r') as f:
      content = f.read()
      tweak[-1].append(float(re.findall(r'= (.*) seconds', content)[0]))


colors = ['blue', 'purple', 'orange', 'brown', 'green', 'pink', 'red', 'cyan', 'olive']

for i in range(4):
  yval = []
  for j in range(len(t)):
    yval.append(tweak[j][i])
  plt.plot(range(8), yval, marker='o', markersize=4, color = colors[i * 2], label = '#particles/#cores = ' + str(r[i]))

plt.legend()
plt.ylim(0.01, 500)
plt.yscale('log')
plt.xticks(range(8), t)
plt.xlabel('Number of Cores')
plt.ylabel('Time (s)', ha='left', y=1, rotation=0)
plt.tight_layout()
plt.savefig('weak.png')
