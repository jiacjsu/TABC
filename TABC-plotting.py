import numpy as np
import matplotlib.pyplot as plt
import matplotlib.artist as art


path = '/Users/chunjing/Documents/Research_Projects/Sqw/TABC-data/'
labelSize = 18
titleSize = 20
xstart = 0
xend = 4
fac = 30

names = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09']
for i in range(10,100):
    names.append(str(i))   

kpoints_8site = ['(0,0)', '(pi,0)', '(pi,pi)', '(pi/2,pi/2)']
fig = plt.figure(figsize=(12,12)) 
#fig2 = plt.figure(figsize=(6,6))
#ax2 = fig2.add_subplot(111)
for j in range(4):   
    for i in range(100):
        ax = fig.add_subplot(10,10,1+i)
        data = np.loadtxt(path + '8_2e/Sqw-' + names[i] + '-' + str(j) + '.dat')
        ax.plot(data[:,0], data[:,1] + j*fac, color = 'k')
        if i == 0: 
            xas = data
        else:
            xas = xas + data
    ax.set_xlim([xstart, xend])

    #ax2.plot(xas[:,0]/100, xas[:,1]/100+j*fac, label = kpoints_8site[j])    
    #ax2.set_xlim([xstart, xend])
plt.legend()
plt.title('U = 8t, t$^\prime$ = -0.30t, 8 site with 2 electrons')            
plt.show()
