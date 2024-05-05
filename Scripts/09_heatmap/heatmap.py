# continue from 08_confusion_matrix

import math
matrix = [[0 for i in range(101)] for j in range(121)]
for length in range(0,121,1):
  for identity in range(0,101,1):
    m1 = df[length][identity]
    tp = m1[0]
    tn = m1[1]
    fp = m1[2]
    fn = m1[3]
    if (tp+fp)==0:
      matrix[length][identity] = 0
    else:
      recall = tp/(tp+fn)
      pre = tp/(tp+fp)
      mcc = ((tp*tn)-(fp*fn))/math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
      f1 = tp/(tp+0.5*(fp+fn))
      spe= tn/(tn+fp)
      acc = (tn+tp)/(tn+fp+fn+tp)
      matrix[length][identity] = mcc
matrix = pd.DataFrame.from_dict(matrix)
import matplotlib.pyplot as plt
x_label = []
y_label = []
for x in range(101):
  x_label.append(x)
for y in range(121):
  y_label.append(y)
matrix = pd.DataFrame.from_dict(matrix)
ax = sn.heatmap(matrix, cmap="YlGnBu")
ax.invert_yaxis()
plt.title('MCC’')
plt.xlabel('Identity')
plt.ylabel('Length')
plt.show()
