import numpy as np
eps1 = 1.0e-4
eps2 = 1.0e-5
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
    value=round(sum([OVERLAP[i][jpol]*ROT[kpol][i] for i in range(3)]),4)
    print value
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



##########
set_sym_bl([[ 6.9313   ,  -4.0018  ,    0.0000], [0.0000    ,  8.0036  ,    0.0000 ], [0.0000  ,    0.0000    ,  7.1834]])
