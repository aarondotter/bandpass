from pylab import *

x=genfromtxt('Robo-AO_QE_curves_actual.txt',skip_header=2)

f=open('LP600.txt','w')
for i in range(651):
    f.write('{0} {1}\n'.format(10*x[i,0],max(x[i,1],0)))
f.close()

f=open('g.txt','w')
for i in range(651):
    f.write('{0} {1}\n'.format(10*x[i,0],max(x[i,2],0)))
f.close()

f=open('r.txt','w')
for i in range(651):
    f.write('{0} {1}\n'.format(10*x[i,0],max(x[i,3],0)))
f.close()

f=open('i.txt','w')
for i in range(651):
    f.write('{0} {1}\n'.format(10*x[i,0],max(x[i,4],0)))
f.close()

f=open('z.txt','w')
for i in range(651):
    f.write('{0} {1}\n'.format(10*x[i,0],max(x[i,5],0)))
f.close()
