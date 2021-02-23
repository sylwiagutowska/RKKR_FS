import numpy as np
import os
from operator import itemgetter
eps1 = 1.0e-4
PRECIS=int(-np.log10(eps1))
print PRECIS
eps2 = 1.0e-5
nkp=24
try: nkp=int(raw_input('Number of calc. points in one direction: '))
except: nkp=int(input('Number of calc. points in one direction: '))

COMMAND=''
try: yn=(raw_input('Make calculations (y) or analysis only (n)? [n]'))
except: yn=(input('Make calculations (y) or analysis only (n)? [n]'))
if 'y' in yn:   COMMAND='rm fileo*;~/programy/RKKR/rkkr-hop/rkkr004-gf >out'

mag=0
#try: yn=(raw_input('Magnetic calculations (y) or not (n)? [n]'))
#except: yn=(input('Magnetic calculations (y) or not (n)? [n]'))
#if 'y' in yn:   mag=1

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


def set_sym_bl(a_vec):
#  ! Provides symmetry operations for all bravais lattices
#  ! Tests first the 24 proper rotations for the cubic lattice;
#  ! then the 8 rotations specific for the hexagonal axis (special axis c);
#  ! then inversion is added
 sin3 = 0.866025403784438597
 cos3 = 0.5
 msin3 =-0.866025403784438597
 mcos3 = -0.5
 # ! s0: the s matrices in cartesian axis
 # ! overlap: inverse overlap matrix between direct lattice
 # ! rat: the rotated of a direct vector ( cartesian )
 # ! rot: the rotated of a direct vector ( crystal axis )
 # ! value: component of the s matrix in axis basis
 # INTEGER :: jpol, kpol, mpol, irot, imat(32)
 # ! counters over the polarizations and the rotations

 S0= [[[ 1.,  0.,  0.],[  0.,  1.,  0.],[  0.,  0.,  1.]], 
          [[-1.,  0.,  0.],[  0., -1.,  0.],[  0.,  0.,  1.]],
          [[-1.,  0.,  0.],[  0.,  1.,  0.],[  0.,  0., -1.]],
          [[ 1.,  0.,  0.],[  0., -1.,  0.],[  0.,  0., -1.]],
          [[ 0.,  1.,  0.],[  1.,  0.,  0.],[  0.,  0., -1.]],
          [[ 0., -1.,  0.],[ -1.,  0.,  0.],[  0.,  0., -1.]],
          [[ 0., -1.,  0.],[  1.,  0.,  0.],[  0.,  0.,  1.]],
          [[ 0.,  1.,  0.],[ -1.,  0.,  0.],[  0.,  0.,  1.]],
          [[ 0.,  0.,  1.],[  0., -1.,  0.],[  1.,  0.,  0.]],
          [[ 0.,  0., -1.],[  0., -1.,  0.],[ -1.,  0.,  0.]],
          [[ 0.,  0., -1.],[  0.,  1.,  0.],[  1.,  0.,  0.]],
          [[ 0.,  0.,  1.],[  0.,  1.,  0.],[ -1.,  0.,  0.]],
          [[-1.,  0.,  0.],[  0.,  0.,  1.],[  0.,  1.,  0.]],
          [[-1.,  0.,  0.],[  0.,  0., -1.],[  0., -1.,  0.]],
          [[ 1.,  0.,  0.],[  0.,  0., -1.],[  0.,  1.,  0.]],
          [[ 1.,  0.,  0.],[  0.,  0.,  1.],[  0., -1.,  0.]],
          [[ 0.,  0.,  1.],[  1.,  0.,  0.],[  0.,  1.,  0.]],
          [[ 0.,  0., -1.],[ -1.,  0.,  0.],[  0.,  1.,  0.]],
          [[ 0.,  0., -1.],[  1.,  0.,  0.],[  0., -1.,  0.]],
          [[ 0.,  0.,  1.],[ -1.,  0.,  0.],[  0., -1.,  0.]],
          [[ 0.,  1.,  0.],[  0.,  0.,  1.],[  1.,  0.,  0.]],
          [[ 0., -1.,  0.],[  0.,  0., -1.],[  1.,  0.,  0.]],
          [[ 0., -1.,  0.],[  0.,  0.,  1.],[ -1.,  0.,  0.]],
          [[ 0.,  1.,  0.],[  0.,  0., -1.],[ -1.,  0.,  0.]],
          [[ cos3,  sin3, 0.],[ msin3,  cos3, 0.],[ 0., 0.,  1.]],
          [[ cos3, msin3, 0.],[  sin3,  cos3, 0.],[ 0., 0.,  1.]],
          [[mcos3,  sin3, 0.],[ msin3, mcos3, 0.],[ 0., 0.,  1.]],
          [[mcos3, msin3, 0.],[  sin3, mcos3, 0.],[ 0., 0.,  1.]],
          [[ cos3, msin3, 0.],[ msin3, mcos3, 0.],[ 0., 0., -1.]],
          [[ cos3,  sin3, 0.],[  sin3, mcos3, 0.],[ 0., 0., -1.]],
          [[mcos3, msin3, 0.],[ msin3,  cos3, 0.],[ 0., 0., -1.]],
          [[mcos3,  sin3, 0.],[  sin3,  cos3, 0.],[ 0., 0., -1.]]]

 S0NAME=['identity                                     ',
                '180 deg rotation - cart. axis [0,0,1]        ',
                '180 deg rotation - cart. axis [0,1,0]        ',
                '180 deg rotation - cart. axis [1,0,0]        ',
                '180 deg rotation - cart. axis [1,1,0]        ',
                '180 deg rotation - cart. axis [1,-1,0]       ',
                ' 90 deg rotation - cart. axis [0,0,-1]       ',
                ' 90 deg rotation - cart. axis [0,0,1]        ',
                '180 deg rotation - cart. axis [1,0,1]        ',
                '180 deg rotation - cart. axis [-1,0,1]       ',
                ' 90 deg rotation - cart. axis [0,1,0]        ',
                ' 90 deg rotation - cart. axis [0,-1,0]       ',
                '180 deg rotation - cart. axis [0,1,1]        ',
                '180 deg rotation - cart. axis [0,1,-1]       ',
                ' 90 deg rotation - cart. axis [-1,0,0]       ',
                ' 90 deg rotation - cart. axis [1,0,0]        ',
                '120 deg rotation - cart. axis [-1,-1,-1]     ',
                '120 deg rotation - cart. axis [-1,1,1]       ',
                '120 deg rotation - cart. axis [1,1,-1]       ',
                '120 deg rotation - cart. axis [1,-1,1]       ',
                '120 deg rotation - cart. axis [1,1,1]        ',
                '120 deg rotation - cart. axis [-1,1,-1]      ',
                '120 deg rotation - cart. axis [1,-1,-1]      ',
                '120 deg rotation - cart. axis [-1,-1,1]      ',
                ' 60 deg rotation - cryst. axis [0,0,1]       ',
                ' 60 deg rotation - cryst. axis [0,0,-1]      ',
                '120 deg rotation - cryst. axis [0,0,1]       ',
                '120 deg rotation - cryst. axis [0,0,-1]      ',
                '180 deg rotation - cryst. axis [1,-1,0]      ',
                '180 deg rotation - cryst. axis [2,1,0]       ',
                '180 deg rotation - cryst. axis [0,1,0]       ',
                '180 deg rotation - cryst. axis [1,1,0]       ',
                'inversion                                    ',
                'inv. 180 deg rotation - cart. axis [0,0,1]   ',
                'inv. 180 deg rotation - cart. axis [0,1,0]   ',
                'inv. 180 deg rotation - cart. axis [1,0,0]   ',
                'inv. 180 deg rotation - cart. axis [1,1,0]   ',
                'inv. 180 deg rotation - cart. axis [1,-1,0]  ',
                'inv.  90 deg rotation - cart. axis [0,0,-1]  ',
                'inv.  90 deg rotation - cart. axis [0,0,1]   ',
                'inv. 180 deg rotation - cart. axis [1,0,1]   ',
                'inv. 180 deg rotation - cart. axis [-1,0,1]  ',
                'inv.  90 deg rotation - cart. axis [0,1,0]   ',
                'inv.  90 deg rotation - cart. axis [0,-1,0]  ',
                'inv. 180 deg rotation - cart. axis [0,1,1]   ',
                'inv. 180 deg rotation - cart. axis [0,1,-1]  ',
                'inv.  90 deg rotation - cart. axis [-1,0,0]  ',
                'inv.  90 deg rotation - cart. axis [1,0,0]   ',
                'inv. 120 deg rotation - cart. axis [-1,-1,-1]',
                'inv. 120 deg rotation - cart. axis [-1,1,1]  ',
                'inv. 120 deg rotation - cart. axis [1,1,-1]  ',
                'inv. 120 deg rotation - cart. axis [1,-1,1]  ',
                'inv. 120 deg rotation - cart. axis [1,1,1]   ',
                'inv. 120 deg rotation - cart. axis [-1,1,-1] ',
                'inv. 120 deg rotation - cart. axis [1,-1,-1] ',
                'inv. 120 deg rotation - cart. axis [-1,-1,1] ',
                'inv.  60 deg rotation - cryst. axis [0,0,1]  ',
                'inv.  60 deg rotation - cryst. axis [0,0,-1] ',
                'inv. 120 deg rotation - cryst. axis [0,0,1]  ',
                'inv. 120 deg rotation - cryst. axis [0,0,-1] ',
                'inv. 180 deg rotation - cryst. axis [1,-1,0] ',
                'inv. 180 deg rotation - cryst. axis [2,1,0]  ',
                'inv. 180 deg rotation - cryst. axis [0,1,0]  ',
                'inv. 180 deg rotation - cryst. axis [1,1,0]  ' ]

#####finding the symmetries 
#    compute the overlap matrix for crystal axis
 ROT=np.array([[sum( [a_vec[kpol][i]*a_vec[jpol][i] for i in range(3)]) for kpol in range(3)] for jpol in range(3)])
 OVERLAP=np.linalg.inv(ROT)
 S=np.zeros((48,3,3))
 SNAME=[ [] for i in range(48)]
 IMAT=[ 0 for i in range(32)]
 nrot=0
 for irot in range(32):
  isign=0
  for jpol in range(3):
   #compute, in cartesian coordinates the rotated vector
   RAT=[ sum([S0[irot][i][mpol]*a_vec[jpol][i] for i in range(3)]) for mpol in range(3)]
   ROT[jpol]=[ sum([a_vec[kpol][i]*RAT[i] for i in range(3)]) for kpol in range(3)]
  # and the inverse of the overlap matrix is applied
  for jpol in range(3):
   for kpol in range(3):
    value=round(sum([OVERLAP[i][jpol]*ROT[kpol][i] for i in range(3)]),PRECIS)
    if abs((int(value))-value)>eps1:
              # if a noninteger is obtained, this implies that this operation
              # is not a symmetry operation for the given lattice
     isign=1
     break
    else: S[nrot][jpol][kpol]=(int(value))
   if isign==1: break
  if isign==1: continue
  SNAME[nrot]=S0NAME[irot]
  IMAT[nrot]=irot
  nrot=nrot+1
#####end of finding the symmetries

 #check number of symmetries
 temp_numbers=[1,2,4,6,8,12,24]
 print nrot
 if nrot not in temp_numbers:
  print("NOTICE: Bravais lattice has wrong number  of symmetries - symmetries are disabled")
  nrot = 0

 #set the inversion symmetry ( Bravais lattices have always inversion
 # !     symmetry )
 for irot in range(nrot):
  SNAME[irot+nrot]=S0NAME[IMAT[irot]+32]
  for kpol in range(3):
   for jpol in range(3):
    S[irot+nrot][jpol][kpol]=-S[irot][jpol][kpol]
 nrot=2*nrot

 print("Found "+str(nrot)+" symmetry operations")
 return  S[:nrot]

def sym_of_group(SYMM_OP):
 SYMM_OP2=[SYMM_OP[0]]
 h=open('filei5','r')
 tmp=[i.split() for i in h.readlines()[2:]]
 h.close()
 for ni,i in enumerate(tmp):
  if len(i)==1: 
   nat=int(i[0])
   break
 ni+=1
 POS=[]
 for i in range(nat):
   npos=int(tmp[ni][0])
   print npos
   ni+=1
   while len(tmp[ni])!=3: ni+=1
   for j in range(npos):
    POS.append([[float(m) for m in tmp[ni]],i])
    ni+=1
 print POS
 for i in SYMM_OP[1:]:
  sign=0
  for nj,j in enumerate(POS):
   pos_new=[ sum([i[k][m]*j[0][m] for m in range(3)]) for k in range (3)]
   for k in POS[nj+1:]:
    if j[1]!=k[1] and sum([ (pos_new[m]-k[0][m])**2])<1e-5: 
     sign=1
     break
   if sign==1: break
  if sign==0:
   SYMM_OP2.append(i)
 print("Found "+str(len(SYMM_OP2))+" symmetry operations")
 return SYMM_OP2

def ruotaijk(s,k):
  return [s[0][0]*k[0]+s[0][1]*k[1]+s[0][2]*k[2],
  s[1][0]*k[0]+s[1][1]*k[1]+s[1][2]*k[2],
  s[2][0]*k[0]+s[2][1]*k[1]+s[2][2]*k[2]]


def check(n,k,kw,ieq,SYM_OP,nmax):
  flag=1
  for s in SYM_OP:
     kr=ruotaijk( s,k[n] )
     for j in range(3):
        while kr[j]>=nmax[j]:
           kr[j]=kr[j]-nmax[j]
        while kr[j]<=-1:
           kr[j]=kr[j]+nmax[j]
     for np in range(n): 
        if  kr[0]==k[np][0] and kr[1]==k[np][1] and kr[2]==k[np][2]:
           kw[n]=0
           naux =np
           while (kw[naux]==0):
              naux=ieq[naux]
           ieq[n]=naux
           kw[naux]=kw[naux]+1
           flag=0
           break
     if flag==0: break
   
def make_k_grid(nkp,b_vec,SYMM_OP,VEC,VEC_all):
 n=0
 KVEC,kw,ieq=[],[1 for i in range(nkp*nkp*nkp)],[0 for i in range(nkp*nkp*nkp)]
 nmax=[nkp,nkp,nkp]
 for i in range(nkp):
  print i
  for j in range(nkp):
   for k in range(nkp):
    KVEC.append(([i,j,k])) 
    check(n,KVEC,kw,ieq,SYMM_OP,nmax)
    n=n+1
 NONEQKW=[]
 NONEQK=[]
 nk=0
 for j in range(n):
    if(kw[j]>0):
        NONEQKW.append(kw[j])
        NONEQK.append([ sum([ KVEC[j][i]*b_vec[l][i]/nmax[i] for i in range(3)]) for l in range(3)])
        print j,KVEC[j],kw[j],ieq[j]
        nk=nk+1
 print nk
  
#### END OF FUNCTIONS


####MAIN


####make k-grid
#1. read data
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
#2.Find symmetry operations using list of atoms from filei5
SYMM_OP=set_sym_bl(a_vec)
#SYMM_OP=sym_of_group(SYMM_OP)

#3. Make k-grid
VEC=[]
VEC_all=[]
make_k_grid(nkp,b_vec,SYMM_OP)

'''
for i in range(nkp):
 print i
 for j in range(nkp):
  for k in range(nkp):
   #v=(([i/float(nkp)-0.5,j/float(nkp)-0.5,k/float(nkp)-0.5])) 
   #v=[round(n,3) for n in sum([v[m]*b_vec[m] for m in range(3)])]
   v=[i,j,k]
   which_k=len(VEC)
   #check symmetry
   sign1=0
   for s in SYMM_OP:
    v_sym=[ sum([s[n][m]*v[m] for m in range(3)]) for n in range (3)]
    for j in range(3):
        while v_sym[j]>=nkp:
           v_sym[j]=v_sym[j]-nkp
        while v_sym[j]<=-1:
           v_sym[j]=v_sym[j]+nkp
    for (nv_old,v_old) in enumerate(VEC):
     if v_sym[0]-v_old[0] == 0 and v_sym[2]-v_old[2] == 0 and v_sym[1]-v_old[1] == 0 :
      which_k=nv_old
      sign1=1
      break
    if sign1==1: break
   if sign1==0: VEC.append(v)
   VEC_all.append([v,which_k])
'''
'''
#3. Make k-grid 
VEC=[]
VEC_all=[]
for i in range(nkp):
 for j in range(nkp):
  for k in range(nkp):
   v=(([i/float(nkp)-0.5,j/float(nkp)-0.5,k/float(nkp)-0.5])) 
   v=[round(n,3) for n in sum([v[m]*b_vec[m] for m in range(3)])]
   VEC_all.append(v)

VEC_all2=[]
nk=len(VEC_all)
for nki,v in enumerate(VEC_all):
   print nki
   VEC.append([v,nki])
   VEC_all2.append([v,nki])
   #check symmetry
   lva=len(VEC_all)
   for nr in range(1,lva-nki-1):
    v2=VEC_all[lva-nr] #od konca # in enumerate(VEC_all[-1:nki]):
    for s in SYMM_OP:
     v_sym=[ sum([s[n][m]*v2[m] for m in range(3)]) for n in range (3)]
     if sum([(v_sym[m]-v[m])**2 for m in range(3)])<eps1:
      VEC_all2.append([v2,len(VEC)])
      del VEC_all[lva-nr]
      break
   if nki==len(VEC_all): break

VEC_all=VEC_all2
'''

print 'No of points=',len(VEC_all), 'no of irreducible points=',len(VEC)

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
   if i%2==0: h.write(str(2)+' ')
   h.write(str(v[0])+' '+str(v[1])+' '+str(v[2])+'\n')
h.close()
print('I am running calculations....')
os.system(COMMAND)
print('...Calc done')
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
 




