import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure
import re

figure(figsize=(5, 4), dpi=160)

n = [1000, 10000, 100000, 1000000, 6000000]
t = [1, 2, 4, 8, 16, 32, 64, 128]

tstrong = []

for j in t:
  tstrong.append([])
  for i in n:
    with open('../data_strong/t' + str(j) + '-' + str(i) + '.out', 'r') as f:
      content = f.read()
      tstrong[-1].append(float(re.findall(r'= (.*) seconds', content)[0]))


colors = ['blue', 'olive', 'orange', 'brown', 'green', 'pink', 'red', 'cyan', 'purple']

for i in range(len(n)):
  yval = []
  for j in range(len(t)):
    yval.append(tstrong[j][i])
  plt.plot(range(8), yval, marker='o', markersize=4, color = colors[i * 2], label = '#particles = ' + str(n[i]))

plt.plot([0, 7], [100, 100 / 128], '--', color = 'black', label = 'Slope = -1')

plt.legend(fontsize = 7)
plt.ylim(0.01, 5000)
plt.yscale('log')
plt.xticks(range(8), t)
plt.xlabel('Number of Cores')
plt.ylabel('Time (s)', ha='left', y=1, rotation=0)
plt.tight_layout()
plt.savefig('strong.png')
