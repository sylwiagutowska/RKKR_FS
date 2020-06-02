import numpy as np
import os
try: from ctypes import *
except: "Ctypes not found"
global eps1
eps1 = 1.0e-4
global PRECIS
PRECIS=int(-np.log10(eps1))

global eps2
eps2 = 1.0e-5




def make_kgrid_C_lib(nkp,SYMM_OP,b_vec,file):
 lib = CDLL(file)
 ll=len(SYMM_OP)*3*3
 lib.FS_make_kgrid.argtypes = c_int,c_int,np.ctypeslib.ndpointer(dtype=np.double)
 lib.FS_make_kgrid.restype =np.ctypeslib.ndpointer(dtype=c_int, shape=(2,nkp*nkp*nkp))
 SYMM_OP2 = []
 for i in SYMM_OP:
  for j in range(3):
   for k in range(3):
    SYMM_OP2.append(i[k][j])
 SYMM_OP2=np.array(SYMM_OP2, dtype=np.double)
 ieq,kw=lib.FS_make_kgrid(c_int(nkp),c_int(len(SYMM_OP)),SYMM_OP2)
 print("Python confirms: C ended")
 KVEC=[]
 for i in range(nkp):
  for j in range(nkp):
   for k in range(nkp):
    KVEC.append([i,j,k])
 NONEQKW=[]
 NONEQK=[]
 ieq_new=[]
 KVEC=[ [ round(sum([v[m]*b_vec[m2][m]/nkp for m in range(3)]),4) for m2 in range(3)] for v in KVEC]
 nk=0
 for j in range(nkp*nkp*nkp):
    if kw[j]!=-1: 
        NONEQKW.append(kw[j])
        NONEQK.append(KVEC[j])
        KVEC[j].append(nk)
        ieq_new.append(nk)
        nk=nk+1
    else: KVEC[j].append(KVEC[ieq[j]][3])
 return NONEQK,KVEC



def ask_for_input(COMMAND):
 nkp=24
 liczba=3
 try: nkp=int(raw_input('Number of calc. points in one direction: '))
 except: nkp=int(input('Number of calc. points in one direction: '))
 try: yn=(raw_input('Make calculations (y) or analysis only (n)? [n]'))
 except: yn=(input('Make calculations (y) or analysis only (n)? [n]'))
 if 'y' not in yn:   
  COMMAND=''
  try: liczba=int(raw_input('Rkkr (3)  or Rcpa (4)? Type 3 or 4 : '))
  except: liczba=int(input('Rkkr (3)  or Rcpa (4)? Type 3 or 4 : '))
 else:
  if 'cpa' in COMMAND: liczba=4
  elif 'kkr' in COMMAND: liczba=3
  else:
   try: liczba=int(raw_input('Rkkr (3)  or Rcpa (4)? Type 3 or 4 : '))
   except: liczba=int(input('Rkkr (3)  or Rcpa (4)? Type 3 or 4 : '))
 return nkp,COMMAND,liczba

def read_data():
 os.system('rm lattice_vectors.dat')
 table=[ 'R'+str(i+1) for i in range(3)]
 for i in table:
  os.system('grep "'+(i)+'" * >>lattice_vectors.dat  2> /dev/null')
 os.system('grep " a'+' ='+'" * >>lattice_vectors.dat  2> /dev/null')
 os.system('grep " b'+' ='+'" * >>lattice_vectors.dat 2> /dev/null')
 os.system('grep " c'+' ='+'" * >>lattice_vectors.dat 2> /dev/null')
 h=open('lattice_vectors.dat')
 xxx=[i for i  in h.readlines() if '0 ' in i]
 alat=[float(xxx[-3+m].split()[-2]) for m in range(3)]
 a_vec=[]
 for i in xxx:
  for nj,j in enumerate(table):
   if j in i and len(a_vec)==nj: a_vec.append(i)
 a_vec=[ [float(j) for j in i.split()[4:7]] for i in a_vec]
 print( a_vec)
 return a_vec,alat

def recip_vec_gen(a_vec,alat):
 br=np.transpose(a_vec)
 bg=np.zeros((3,3))
 for i in range(3):
  i1=(i+1)%3 
  i2=(i1+1)%3
  bg[i]=np.cross( br[i1],br[i2])
 vws=abs(sum([br[0][i]*bg[0][i] for i in range(3)]))
 bg=np.transpose(bg)
 bg=bg/vws
 bg=[ [round(m2) for m2 in bg[m]*alat[m]] for m in range(3)]
 return bg

def run_calc(VEC,alat,COMMAND,file_in, file_in2):
 os.system('cp '+(file_in)+' '+(file_in2))
 h=open(file_in2,'a')
 if len(VEC)%2==0: h.write(str(int(len(VEC)/2))+'\n')
 else:  h.write(str(int(int(len(VEC)/2)+1))+'\n')
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
 print('I am running calculations with command'+COMMAND)
 os.system(COMMAND)
 print('...Calc done')

def read_ene(liczba,file_bands):
#read data
 h=open(file_bands)
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
 print ('from ',max(nbnd),' bands only ',min(nbnd),' are complete. I will approximate the rest.')
 return ENE,nbnd


def ene_from_irreducible_to_whole_grid(ENE,VEC,VEC_all,nkp):
#from energies of irreducible points to energies of whole grid
 ENE2=[]
 for i in ENE:
  for j in i:
    ENE2.append(j)
 print( 'No of irreducible energy points='+str(len(ENE2))+'. Should equal to '+str(len(VEC)))
 ENE=[]
 for k in ((VEC_all)):
     ENE.append(ENE2[k[3]])
 print ('No of energy points='+str(len(ENE))+'. Should equal to '+str(nkp*nkp*nkp))
 h=open('ene0.dat','w')
 for i in ENE:
   for k in i:
    h.write(str(k[0])+' ')
   h.write('\n')
 h.close()
 return ENE

def complete_incomplete_bands(nbnd,ENE):
 #detect incomplete band
 nband0=min(nbnd)
 nband=max(nbnd)
 for i in range(len(ENE)):
  if len(ENE[i])==nband:
   ni_max=i
   last_b=ENE[i][:]
   break
 last_b0=last_b[:]
 for i in range(ni_max,len(ENE)):
  last_b[:len(ENE[i])]=ENE[i][:]
  if len(ENE[i])!=nband: 
   ENE[i]+=last_b[len(ENE[i]):nband]
 for i in range(ni_max):
  last_b0[:len(ENE[ni_max-i-1])]=ENE[ni_max-i-1][:]
  if len(ENE[ni_max-i-1])!=nband: 
   ENE[ni_max-i-1]+=last_b[len(ENE[ni_max-i-1]):nband]
 return nband,ENE

def rearrange_data(nband,ENE):
 ENE=[  [ENE[i][k][:] for i in range(len(ENE))]  for k in range(nband)]
 return ENE

def which_bands_cross_fermi(ENE0):
 ENE=[]
 for i in ENE0:
  isfermi=0
  for nj in range(1,len(i)):
   if i[nj][0]*i[nj-1][0]<=0: 
     isfermi=1
     break
  if isfermi==1: ENE.append(i)
 return ENE

def write_to_file(ENE,VEC,VEC_all,b_vec, alat,
                  nkp,liczba,file_frmsf,file_mayavi):
 h=open(file_frmsf,'w')
 for i in range(3): h.write(str(nkp)+' ')
 h.write('\n1\n'+str(len(ENE))+'\n')
 for i in range(3):
  for j in range(3): h.write(str(b_vec[i][j]/alat[i])+' ')
  h.write('\n')
 for i in ENE:
  for j in i:
    h.write(str(j[0])+'\n')
 for ni,i in enumerate(ENE):
  for j in i:
    if liczba==3: h.write(str(ni)+'\n')
    else: h.write(str(j[1])+'\n')
 #h.write(str(ni)+'\n')
 h.close()

 print('Print Fermi surface to file...')
 h=open(file_mayavi,'w')
 h.write(str(nkp*nkp*nkp)+' '+str(len(ENE))+' '+str(0.0)+' ')
 for b in range(len(ENE)):
  h.write('0 ')
 h.write('\n')
 for i in range(len(VEC_all)):
  h.write(str(VEC_all[i][0])+' '+str(VEC_all[i][1])+' '+str(VEC_all[i][2])+' ')
  for b in range(len(ENE)):
   h.write(str(ENE[b][i][0])+' ')
  h.write('\n')
 h.close()


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
 print (nrot)
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
     kr=ruotaijk( s,k[n] ) 
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
  print( i)
  for j in range(nkp):
   for k in range(nkp):
    KVEC.append([i,j,k])
    check(n,KVEC,kw,ieq,SYMM_OP,nkp)
    n=n+1
 NONEQKW=[]
 NONEQK=[]
 ieq_new=[]
 KVEC=[ [ round(sum([v[m]*b_vec[m2][m]/nkp for m in range(3)]),4) for m2 in range(3)] for v in KVEC]
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

#### END OF FUNCTIONS



