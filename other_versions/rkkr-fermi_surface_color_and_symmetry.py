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
if 'y' in yn:   
 #COMMAND='rm fileo*;~/programy/RKKR/rkkr-hop/rkkr004-gf >out'
 COMMAND='rm fileo*; ~/programy/RKKR/SRC_2017/Rcpa005_17/rcpa005 >out'

liczba=4
if 'cpa' in COMMAND: liczba=4
elif 'kkr' in COMMAND: liczba=3
else:
 try: liczba=int(raw_input('Rkkr (3)  or Rcpa (4)? Type 3 or 4 : '))
 except: liczba=int(input('Rkkr (3)  or Rcpa (4)? Type 3 or 4 : '))


def recip_vec_gen(a_vec):
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

def ruotaijk(s,k):
  return [s[0][0]*k[0]+s[0][1]*k[1]+s[0][2]*k[2],
  s[1][0]*k[0]+s[1][1]*k[1]+s[1][2]*k[2],
  s[2][0]*k[0]+s[2][1]*k[1]+s[2][2]*k[2]]


def check(n,k,kw,ieq,SYM_OP,nmax):
  flag=1
  for s in SYM_OP:
     kr=ruotaijk( s,k[n] ) #[s[0][0]*k[n][0]+s[0][1]*k[n][1]+s[0][2]*k[n][2],
#  s[1][0]*k[n][0]+s[1][1]*k[n][1]+s[1][2]*k[n][2],
#  s[2][0]*k[n][0]+s[2][1]*k[n][1]+s[2][2]*k[n][2]] #

     for j in range(3):
        while kr[j]>=nmax:
           kr[j]=kr[j]-nmax
        while kr[j]<=-1:
           kr[j]=kr[j]+nmax
     for npk in range(n): 
        if abs(kr[0]-k[npk][0])<eps1 and abs(kr[1]-k[npk][1])<eps1 and abs(kr[2]-k[npk][2])<eps1:
           kw[n]=-1
           naux =npk
           while (kw[naux]==-1):
              naux=ieq[naux]
           ieq[n]=naux
           kw[naux]=kw[naux]+1
           flag=0
           break
     if flag==0: break
   
def make_k_grid(nkp,b_vec,SYMM_OP):
 n=0
 SYMM_OP=[np.transpose(np.array(s)) for s in SYMM_OP]
 KVEC,kw,ieq=[],[1 for i in range(nkp*nkp*nkp)],[0 for i in range(nkp*nkp*nkp)]
 for i in range(nkp):
  print i
  for j in range(nkp):
   for k in range(nkp):
    #v=(([i/float(nkp),j/float(nkp),k/float(nkp)]))
    KVEC.append([i,j,k])
    #KVEC.append([ round(m2,3) for m2 in [sum([v[m]*b_vec[m3][m] for m in range(3)]) for m3 in range(3)]])
    check(n,KVEC,kw,ieq,SYMM_OP,nkp)
    n=n+1
 NONEQKW=[]
 NONEQK=[]
 ieq_new=[]
 KVEC=[ [ round(sum([v[m]*b_vec[m2][m]/nkp for m in range(3)]),4) for m2 in range(3)] for v in KVEC]
 #KVEC=[  [ kv[i]/nmax[i] for i in range(3)] for kv in KVEC]
 nk=0
 for j in range(n):
    if kw[j]!=-1: 
        NONEQKW.append(kw[j])
        NONEQK.append(KVEC[j])
        KVEC[j].append(nk)
        ieq_new.append(nk)
        nk=nk+1
    else: KVEC[j].append(KVEC[ieq[j]][3])
 return NONEQK,KVEC

def make_k_grid2(nkp,b_vec,SYMM_OP):
# pm=[0,-1,1,-2,2]
# A_vectors=[[round(m2,3) for m2 in [h*b_vec[0][m]+k*b_vec[1][m]+l*b_vec[2][m] for m in range(3)]] for h in pm for k in pm for l in pm]
 n=0
 SYMM_OP=[np.transpose(np.array(s)) for s in SYMM_OP]
 NONEQK,KVEC,kw,which_k=[],[],[1 for i in range(nkp*nkp*nkp)],[]
 for i in range(nkp):
  print i,':'
  for j in range(nkp):
   print j,len(NONEQK)
   for k in range(nkp):
    v=[i/float(nkp),j/float(nkp),k/float(nkp)]
    #v=[i,j,k]
    #v=([ round(m2,PRECIS) for m2 in [sum([v[m]*b_vec[m3][m]/nkp for m in range(3)]) for m3 in range(3)]])
    KVEC.append(v)
    sign=0
    for s in SYMM_OP:
     for nk2,k2 in enumerate(NONEQK):
      k3=ruotaijk( s,k2 ) 
      if int(k3[0]-v[0])==k3[0]-v[0] and int(k3[1]-v[1])==k3[1]-v[1] and int(k3[2]-v[2])==k3[2]-v[2]:
        sign=1
        which_k.append(nk2)
        break
      if sign==1: break
     if sign==1: break
    if sign==0: 
     NONEQK.append(v)
     which_k.append(n)
     n=n+1
 KVEC=[ [ round(sum([v[m]*b_vec[m2][m] for m in range(3)]),PRECIS) for m2 in range(3)] for v in KVEC]
 NONEQK=[ [ round(sum([v[m]*b_vec[m2][m] for m in range(3)]),PRECIS) for m2 in range(3)] for v in NONEQK]
 for i in range(len(which_k)):
  KVEC[i].append(which_k[i])
 return NONEQK,KVEC
  
'''
def make_whole_k_grid(VEC,SYM_OP,e,nkp):
'''
'''
 equiv=[ 0 for i in range (nkp*nkp*nkp)]
 e=np.linalg.inv(np.array(e))
 for ik in range(len(VEC)):
     #xk_frac = ruotaijk(e,VEC[ik])
     for s in SYM_OP:
        kv = [(m*nkp) for m in ruotaijk(s,VEC[ik])]
        ikv=[int(m) for m in kv]
        if abs(kv[0]-ikv[0])>eps1 or abs(kv[1]-ikv[1])>eps1 or abs(kv[2]-ikv[2])>eps1: continue
        ikv = [ m%nkp for m in ikv]
        equiv[ikv[2]*nkp*nkp+ ikv[1]*nkp+ikv[0]] = ik
 return equiv
'''
'''
 #equivalent kpoints
 allk=[]
 mmm=0
 for i in VEC:
  for j in SYM_OP:
   x=[ round(k,PRECIS) for k in ruotaijk(j,i[0:3]) ]
   #if x[0]<maxe[0] and x[1]<maxe[1] and x[2]<maxe[2] and x[0]>-maxe[0] and x[1]>-maxe[1] and x[2]>-maxe[2]:
   sign=0
   for k in allk:
    if x[0]==k[0] and x[1]==k[1] and x[2]==k[2]:
     sign=1
     break
   if sign==0:
    allk.append(np.array([ x[0],x[1],x[2],int(mmm)]))
  mmm=mmm+1
 allk2=sorting(allk)
 return [ int(i[3]) for i in allk2]
'''
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
a_vec=[]
for i in xxx:
 if 'R1' in i and len(a_vec)==0: a_vec.append(i)
 if 'R2' in i and len(a_vec)==1: a_vec.append(i)
 if 'R3' in i and len(a_vec)==2: a_vec.append(i)
a_vec=[ [float(j) for j in i.split()[4:7]] for i in a_vec]
print a_vec
b_vec= recip_vec_gen(a_vec)
b_vec0=[b_vec[m]/2./np.pi for m in range(3)]
b_vec=[ [round(m2) for m2 in b_vec0[m]*alat[m]] for m in range(3)]
#b_vec2=[[    1.0000, 0.0000, 1.0000],
#  [-1.0000, 1.0000, 0.0000],
#  [ 0.0000,-1.0000, 1.0000]]

print b_vec
#2.Find symmetry operations using lattice vectors from filei5
SYMM_OP=set_sym_bl(a_vec)

#3. Make k-grid
VEC,VEC_all=make_k_grid(nkp,b_vec,SYMM_OP)


print 'No of points=',len(VEC_all), 'no of irreducible points=',len(VEC)

h=open('kvec.dat','w')
for ni,i in enumerate(VEC_all):
 for j in i:
  h.write(str(j)+' ')
 h.write('\n')
h.close()

####run calc
os.system('cp filei50 filei5')
h=open('filei5','a')
'''
h.write(str(len(VEC))+'\n')
for i,v in enumerate(VEC):
   h.write(str(1)+' ')
   h.write(str(v[0])+' '+str(v[1])+' '+str(v[2])+'\n')
'''
if len(VEC)%2==0: h.write(str(len(VEC)/2)+'\n')
else:  h.write(str((len(VEC)/2)+1)+'\n')
if len(VEC)%2!=0:
 for i,v in enumerate(VEC[:-1]):
   if i%2==0: h.write(str(2)+' ')
   h.write(str(v[0])+' '+str(v[1])+' '+str(v[2])+'\n')
 h.write(str(2)+' ')
 h.write(str(VEC[-1][0])+' '+str(VEC[-1][1])+' '+str(VEC[-1][2])+'\n')
 h.write(str(VEC[0][0])+' '+str(VEC[0][1])+' '+str(VEC[0][2])+'\n')
else:
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
  ENE[-1][-1].append([float(line[1].replace('D','E'))-fermi, float(line[2].replace('D','E'))])

print 'from ',max(nbnd),' bands only ',min(nbnd),' are complete. I will approximate the rest.'

#from energies of irreducible points to energies of whole grid
ENE2=[]
for i in ENE:
 for j in i:
   ENE2.append(j)
h=open('ene0_irrep.dat','w')
for ni,i in enumerate(ENE2[:len(VEC)]):
  for k in i:
   h.write(str(k[0])+' ')
  h.write('\n')
h.close()
print 'No of irreducible energy points='+str(len(ENE2))+'. Should equal to '+str(len(VEC))
ENE=[]
for k in ((VEC_all)):
    ENE.append(ENE2[k[3]])
print 'No of energy points='+str(len(ENE))+'. Should equal to '+str(nkp*nkp*nkp)
h=open('ene0.dat','w')
for i in ENE:
  for k in i:
   h.write(str(k[0])+' ')
  h.write('\n')
h.close()


#detect incomplete bands
nband0=min(nbnd)
nband=max(nbnd)
'''
line_len=len(ENE[0][0])
for i in range(len(ENE)):
  if len(ENE[i])!=nband: 
   ENE[i]+=[[1000 for k in range(line_len)] for m in range(nband-len(ENE[i]))]
'''

'''
last_b=[ [] for i in range(nband)]
for i in range(len(ENE)):
  last_b[len(ENE[i])-1]=ENE[i][-1]
  if len(last_b[-1])!=0: break
'''

for i in range(len(ENE)):
 if len(ENE[i])==nband:
  ni_max=i
  last_b=ENE[i][:]
  print 'found!'
  break
last_b0=last_b[:]
for i in range(ni_max,len(ENE)):
 # for nj,j in enumerate(ENE[i]):
 #  last_b[nj]=j
  last_b[:len(ENE[i])]=ENE[i][:]
#  last_b[len(ENE[i])-1]=ENE[i][-1]
  if len(ENE[i])!=nband: 
   #ENE[i][j]+=[[k+1 for k in range(line_len)] for m in range(nband-len(ENE[i][j]))]
   ENE[i]+=last_b[len(ENE[i]):nband]
print ENE[-1]
for i in range(ni_max):
  last_b0[:len(ENE[ni_max-i-1])]=ENE[ni_max-i-1][:]
  if len(ENE[ni_max-i-1])!=nband: 
   ENE[ni_max-i-1]+=last_b[len(ENE[ni_max-i-1]):nband]
print ENE[-1]


h=open('ene.dat','w')
for i in ENE:
  for k in i:
   h.write(str(k[0])+' ')
  h.write('\n')
h.close()


#rearrange data
ENE=[  [ENE[i][k][:] for i in range(len(ENE))]  for k in range(nband)]
'''
ENE2=[ [] for k in range(nband)]
for k in range(nband):
 for i in range(len(ENE)):
   ENE2[k].append(ENE[i][k])
ENE=ENE2[:]
'''
print 'No of bands = ',len(ENE),'. Each of them on the grid of',
for i in ENE: print len(i) ,
print 'kpoints'

#check which bands cross EF
ENE2=[]
for i in ENE:
 isfermi=0
 for nj in range(1,len(i)):
  if i[nj][0]*i[nj-1][0]<=0: 
    isfermi=1
    break
 if isfermi==1: ENE2.append(i)

print 'Number of pieces of Fermi surface=',len(ENE2)

h=open('fermi.frsmf','w')
for i in range(3): h.write(str(nkp)+' ')
h.write('\n1\n'+str(len(ENE2))+'\n')
for i in range(3):
 for j in range(3): h.write(str(b_vec[i][j])+' ')
 h.write('\n')
for i in ENE2:
 for j in i:
   h.write(str(j[0])+'\n')
for ni,i in enumerate(ENE2):
 for j in i:
   if liczba==3: h.write(str(ni)+'\n')
   else: h.write(str(j[1])+'\n')
 #h.write(str(ni)+'\n')
h.close()

VEC_all=sorting(VEC_all)
print('Print Fermi surface to file...')
h=open('FS.dat','w')
h.write(str(nkp*nkp*nkp)+' '+str(len(ENE2))+' '+str(0.0)+' ')
for b in range(len(ENE2)):
 h.write('0 ')
h.write('\n')
for i in range(len(VEC_all)):
 h.write(str(VEC_all[i][0])+' '+str(VEC_all[i][1])+' '+str(VEC_all[i][2])+' ')
 for b in range(len(ENE2)):
  h.write(str(ENE2[b][i][0])+' ')
 h.write('\n')
h.close()



