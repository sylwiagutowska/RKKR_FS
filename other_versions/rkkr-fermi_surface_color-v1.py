import numpy as np
import os
from operator import itemgetter
PRECIS=4
nkp=24
try: nkp=int(raw_input('Number of calc. points in one direction: '))
except: nkp=int(input('Number of calc. points in one direction: '))

COMMAND=''
try: yn=(raw_input('Make calculations (y) or analysis only (n)? [n]'))
except: yn=(input('Make calculations (y) or analysis only (n)? [n]'))
if 'y' in yn:   COMMAND='rm fileo*;~/programy/RKKR/rkkr-hop/rkkr004-gf >out'

mag=0
try: yn=(raw_input('Magnetic calculations (y) or not (n)? [n]'))
except: yn=(input('Magnetic calculations (y) or not (n)? [n]'))
if 'y' in yn:   mag=1

liczba=4
if 'rcpa' in COMMAND: liczba=4
elif 'rkkr' in COMMAND: liczba=3
else:
 try: liczba=int(raw_input('Rkkr (3)  or Rcpa (4)? Type 3 or 4 : '))
 except: liczba=int(input('Rkkr (3)  or Rcpa (4)? Type 3 or 4 : '))


def recip_vec_gen(a_vec):
 #      TWOPI = 2*PI
 #      DO 100 I = 1, 3
 #       I1 = 1 + MOD(I,3)
 #       I2 = 1 + MOD(I1,3)
 #        CALL CROSS(BR(1,I1),BR(2,I1),BR(3,I1),
 #     &            BR(1,I2),BR(2,I2),BR(3,I2),
 #     &            BG(1,I),BG(2,I),BG(3,I))
 # 100  CONTINUE
 #      VWS = DABS(BR(1,1)*BG(1,1)+BR(2,1)*BG(2,1)+BR(3,1)*BG(3,1))
 #      VBZ = TWOPI**3/VWS
 #      DO 200 I = 1, 3
 #       DO 150 J = 1, 3
 #        BG(J,I) = TWOPI*BG(J,I)/VWS
 #     SUBROUTINE CROSS(AX,AY,AZ,BX,BY,BZ,CX,CY,CZ)
 #     IMPLICIT NONE
 #     DOUBLE PRECISION AX,AY,AZ,BX,BY,BZ,CX,CY,CZ
 #C
 #     CX = AY*BZ - BY*AZ
 #     CY = BX*AZ - AX*BZ
 #      CZ = AX*BY - BX*AY
 br=np.transpose(a_vec)
 bg=np.zeros((3,3))
 for i in range(3):
  i1=(i+1)%3
  i2=(i1+1)%3
  bg[i]=np.cross( br[i1],br[i2])
 vws=abs(sum([br[0][i]*bg[0][i] for i in range(3)]))
 bg=np.transpose(bg)
 bg=2*np.pi*bg/vws
 return bg
 #should be:
 #   0       r1 =    0.9065      0.0000      0.0000                   0
 #   0       r2 =    0.4532      0.7850     -0.0000                   0
 #    0       r3 =   -0.0000      0.0000      0.8747     

def sorting(allk2):
 Xall=[]
 allk2=sorted(allk2, key=itemgetter(0))
 i=0
 while i<len(allk2): 
  X=[]
  x=allk2[i]
  while i<len(allk2) and x[0]==allk2[i][0]:
   X.append(allk2[i])
   i=i+1
  if len(X)>1: X=sorted(X, key=itemgetter(1))
  Xall.append(X)

 Yall=[]
 for X in Xall:
  x=X[0]
  i=0
  while i<len(X): 
   Y=[]
   x=X[i]
   while i<len(X) and x[1]==X[i][1]:
    Y.append(X[i])
    i=i+1
   Y=sorted(Y, key=itemgetter(2))
   Yall.append(Y)

 allk=[]
 for i in Yall:
  for j in i:
   allk.append(j)
 print (' Sorting - Done!')

 return allk



####MAIN

####make k-grid
table=[" R"+str(i+1) for i in range(3)]
os.system('rm lattice_vectors.dat')
for i in table:
 os.system('grep "'+(i)+'" * >>lattice_vectors.dat')
os.system('grep " a'+' ='+'" * >>lattice_vectors.dat')
os.system('grep " b'+' ='+'" * >>lattice_vectors.dat')
os.system('grep " c'+' ='+'" * >>lattice_vectors.dat')
h=open('lattice_vectors.dat')
xxx=[i for i  in h.readlines() if '0 ' in i]
alat=[float(xxx[-3+m].split()[-2]) for m in range(3)]
a_vec=[ [float(j) for j in i.split()[4:7]] for i in xxx if len(i.split())==8]
print a_vec
b_vec= recip_vec_gen(a_vec)
b_vec0=[b_vec[m]/2./np.pi for m in range(3)]
b_vec=[ b_vec0[m]*alat[m] for m in range(3)]

print b_vec
VEC=[]
for i in range(nkp):
 for j in range(nkp):
  for k in range(2):
   v=np.transpose(np.array([i/float(nkp),j/float(nkp),k-k/float(nkp)]))
   VEC.append([ round(n,3) for n in sum([v[m]*b_vec[m] for m in range(3)])])
h=open('kvec.dat','w')
for i in VEC:
 for j in i:
  h.write(str(j)+' ')
 h.write('\n')
h.close()
####run calc
os.system('cp filei50 filei5')
h=open('filei5','a')
h.write(str(len(VEC)/2)+'\n')
for i,v in enumerate(VEC):
   if i%2==0: h.write(str(nkp)+' ')
   h.write(str(v[0])+' '+str(v[1])+' '+str(v[2])+'\n')
h.close()
os.system(COMMAND)
print('calc done')
#####


#read data
h=open('fileo9')
tmp0=h.readlines()[1:]
tmp=[i.split() for i in tmp0]
h.close()
[n_dir,fermi]=[int(tmp[0][0]), float(tmp[0][1])]
ENE=[]
nbnd=[]
n_dir2=[]
for nl,line in enumerate(tmp[1:]):
 if len(line)==4 and '      ' in tmp0[nl+1]: 
  ENE.append([])
  n_dir2.append(int(line[0]))
 elif (len(ENE[-1])==0) and len(line)==3: continue
 elif len(line)==2:
  if len(ENE[-1])!=0 and len(ENE[-1][-1])!=nbnd[-1]: print('Incorrect number of bands'+str(len(ENE[-1][-1]))+'vs.'+str(nbnd))
  nbnd.append(int(line[1]))
  ENE[-1].append([])
 elif len(line)==liczba:
  if mag==0: ENE[-1][-1].append([float(line[1].replace('D','E'))-fermi, float(line[2].replace('D','E'))])
  else: ENE[-1][-1].append([float(line[1].replace('D','E'))-fermi, float(line[2].replace('D','E'))-fermi, float(line[3].replace('D','E'))])
if len(ENE)!=n_dir: print('Incorrect number of kpoints')

h=open('ene.dat','w')
for i in ENE:
 for j in i:
  for k in j:
   h.write(str(k[0])+' ')
  h.write('\n')
h.close()


print 'from ',max(nbnd),' bands only ',min(nbnd),' are complete. I will approximate the rest.'

line_len=len(ENE[0][0][0])
#detect incomplete bands
nband0=min(nbnd)
nband=max(nbnd)
for i in range(len(ENE)):
 for j in range(len(ENE[i])):
  if len(ENE[i][j])!=nband: 
   ENE[i][j]+=[[k+1 for k in range(line_len)] for k in range(nband-len(ENE[i][j]))]

#rearrange data
#ENE=([ [ j for j in i] for i in ENE])
ENE=[ [ [ENE[i][j][k] for j in range(n_dir2[i]) ] for i in range(n_dir)] for k in range(nband)]


print len(ENE)
for i in ENE: print len(i) 

#check which bands cross EF
ENE2=[]
for i in ENE:
 isfermi=0
 for j in i:
  for nk in range(1,len(j)): 
   if j[nk][0]*j[nk-1][0]<=0 or (mag==1 and j[nk][1]*j[nk-1][1]<=0): 
    isfermi=1
    break
  if isfermi==1: break
 if isfermi==1: ENE2.append(i)

print 'Number of pieces of Fermi surface=',len(ENE2)
###add edges
#for i in range(len(ENE2)):
# for j in range(len(ENE2[i])):
#  ENE2[i][j].append(ENE2[i][j][0])
# ENE[i].append(ENE[i][0])

if mag==0:
 h=open('fermi.frsmf','w')
 for i in range(3): h.write(str(n_dir2[0])+' ')
 h.write('\n1\n'+str(len(ENE2))+'\n')
 for i in range(3):
  for j in range(3): h.write(str(b_vec0[i][j])+' ')
  h.write('\n')
 for i in ENE2:
  for j in i:
   for k in j:
    h.write(str(k[0])+'\n')
 for ni,i in enumerate(ENE2):
  for j in i:
   for k in j:
    if liczba==3: h.write(str(ni)+'\n')
    else: h.write(str(k[1])+'\n')
 #h.write(str(ni)+'\n')
 h.close()

else:
 h=open('fermi_up.frsmf','w')
 for i in range(3): h.write(str(n_dir2[0])+' ')
 h.write('\n1\n'+str(len(ENE2))+'\n')
 for i in range(3):
  for j in range(3): h.write(str(b_vec0[i][j])+' ')
  h.write('\n')
 for i in ENE2:
  for j in i:
   for k in j:
    h.write(str(k[0])+'\n')
 for ni,i in enumerate(ENE2):
  for j in i:
   for k in j:
    if liczba==3: h.write(str(ni)+'\n')
    else: h.write(str(k[2])+'\n')
 #h.write(str(ni)+'\n')
 h.close()
 h=open('fermi_dn.frsmf','w')
 for i in range(3): h.write(str(n_dir2[0])+' ')
 h.write('\n1\n'+str(len(ENE2))+'\n')
 for i in range(3):
  for j in range(3): h.write(str(b_vec0[i][j])+' ')
  h.write('\n')
 for i in ENE2:
  for j in i:
   for k in j:
    h.write(str(k[1])+'\n')
 for ni,i in enumerate(ENE2):
  for j in i:
   for k in j:
    if liczba==3: h.write(str(ni)+'\n')
    else: h.write(str(k[2])+'\n')
 #h.write(str(ni)+'\n')
 h.close()



