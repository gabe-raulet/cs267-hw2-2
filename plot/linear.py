import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure
import re

figure(figsize=(5, 4), dpi=160)

n = [1000, 6000, 10000, 60000, 100000, 600000, 1000000, 6000000]
t = [1, 4, 16, 64, 128]

tlinear = []

for j in t:
  tlinear.append([])
  for i in n:
    with open('../data/t' + str(j) + '-' + str(i) + '.out', 'r') as f:
      content = f.read()
      tlinear[-1].append(float(re.findall(r'= (.*) seconds', content)[0]))


colors = ['blue', 'olive', 'orange', 'brown', 'green', 'pink', 'red', 'cyan', 'purple']

for i in range(len(t)):
  plt.loglog(n, tlinear[i], marker='o', markersize=4, color = colors[i * 2], label = '#cores = ' + str(t[i]))

plt.loglog([1000, 6000000], [0.1, 600], '--', color = 'black', label = 'Slope = 1')

plt.legend()
plt.xlim(500, 10000000)
plt.ylim(0.01, 2000)
plt.xlabel('Number of Particles')
plt.ylabel('Time (s)', ha='left', y=1, rotation=0)
plt.tight_layout()
plt.savefig('linear.png')
