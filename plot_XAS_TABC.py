import numpy as np
import matplotlib.pyplot as plt
import matplotlib.artist as art



labelSize = 18
titleSize = 20
xstart = 925
xend = 938
eoffset = 9.2

names = [' 0', ' 1', ' 2', ' 3', ' 4', ' 5', ' 6', ' 7', ' 8', ' 9']
for i in range(10,100):
    names.append(str(i))   


fig = plt.figure(figsize=(12,12))
for i in range(100):
    ax = fig.add_subplot(10,10,1+i)
    xas0 = np.loadtxt('XAS-' + names[i] + '.dat')
    ax.plot(xas0[:2001,0]+eoffset, xas0[:2001,1], color = 'k')
    if i == 0: 
        xas = xas0
    else:
        xas = xas + xas0
    ax.set_xlim([xstart, xend])
    
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
ax.plot(xas[:,0]/100+eoffset, xas[:,1]/100)    
ax.set_xlim([xstart, xend])
            
plt.show()
