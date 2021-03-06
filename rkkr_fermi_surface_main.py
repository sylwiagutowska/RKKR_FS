#! /usr/bin/python
import numpy as np
import os
import sys
import rkkr_fermi_surface_functions as func
##### 0. Be sure, that none of bands is cutted from the bottom (the bottom energy HAS TO be set in the band gap). 
#####    Dont worry if bands are cutted from the top - I can reconstruct it (the top energy HAS NOT TO be set in the band gap)
##### 1. CHANGE COMMAND VARIABLE
##### 2. IN DIRECTORY PUT
#####    a) OUTPUT FILE FROM PREVIOUS CYCLE OF RKKR/RCPA (THE NAME OF FILE DOES NOT MATTER, BECAUSE I MAKE GREP *)
#####    b) copy filei5 or fileo7 to filei50 AND change last number to "999"
##### 3. CALL this program with python2 rkkr-fermi_surface_color_and_symmetry_obiect.py and answer the questions
##### 4. CALL fermisurfer fermi.frsmf

##### Format of file_frmsf:
#first line: nkp x nkp x nkp
#second line: control number
#third line: number of bands
#4-6 lines: reciprocal lattice vectors (in 1/a_B)
#then list of columns:
# 1. E-ef(eV), 2. linewidth (eV) (if rcpa) or band index (if rkkr)
##### Format of file_lifetime_frmsf:
#1-6 lines like in file_frmsf; then
# 1. E-ef(eV), 2. lifetime (ps) (if rcpa) or band index (if rkkr)
##### Format of file_mayavi
#first line: number of kpoints and number of bands; then:
# 1. kx, 2. ky, 3.kz (in 1/a_i units, a_i-lattice constants), 4. E-ef (eV)
class init_calc:
    def __init__(self):
        self.COMMAND_rkkr = 'rm fileo*;~/programy/RKKR/rkkr-hop/rkkr004-gf '
        self.COMMAND_rcpa = 'rm fileo*;/home/syllwinka/programy/RKKR/SRC_2017/Rcpa005_17/rcpa005 '
        self.clib_name="libFS_make_kgrid.so"
        self.clib_dir=os.path.dirname(sys.argv[0])+"/cpp/FS_make_kgrid/"
        self.clib_file="FS_make_kgrid.cpp"
        self.nkp = 0 #number of points in 1 direction
        self.liczba = 0 #3 for rkkr and 4 for rcpa
        self.file_in = 'filei50' #input file (ended with 999)
        self.file_in2 = 'filei5' #=filei50 with list of kpoints added
        self.file_bands = 'fileo9' #rkkr/rcpa output file with bands
        self.file_out = 'FS_out' #my rkkr/rcpa output file with general info
        self.COMMAND=''
        self.file_frmsf = 'fermi.frmsf' #file where data are written (with bandwith if rcpa is used, with the band index in case of rkkr), readable by fermisurfer
        self.file_lifetime_frmsf = 'fermi_lifetime.frmsf' #file where the Fermi surface coloured by lifetime is written, readable by fermisurfer
        self.file_mayavi = 'FS.dat' #file where data are written, readable by my program FS (which uses mayavi library)
        self.imE_to_lifetime_conv=32.91e-17*1e12 #t= conv/Im(E), E in Ry, t in ps
    def ask_for_input(self,file_out):
        self.nkp,self.COMMAND,self.liczba =     \
            func.ask_for_input(self.COMMAND_rkkr,self.COMMAND_rcpa,\
                               self.file_out)


class Obiekt(init_calc):
    def __init__(self):
        init_calc.__init__(self)
        self.a_vec =[] #direct lattice vector
        self.alat=[] #lattice constants
        self.VEC = [] #list of irreducible points
        self.VEC_all = [] #list of all points
        self.ENE=[] #list of energies (from fileo9)
        self.SYMM_OP=[] #symmetry operations
        self.b_vec=[] #reciprocal lattice vectors
        self.nbnd=[] #number of bands in all k-points
        self.nband=0 #=max(nbnd)=number of bands in every k-point. If nbnd[kpoint]<nband, then the lacked band is approximated in this point
    def read_data(self):
        self.a_vec,self.alat = func.read_data()
    def set_sym_bl(self):
        self.SYMM_OP = func.set_sym_bl(self.a_vec)
    def recip_vec_gen(self):
        self.b_vec = func.recip_vec_gen(self.a_vec,self.alat)
    def make_k_grid(self):
        print( "I am making k-grid...")
        #function written in C is much faster, but ctypes lib is needed. If not present, python function is used.
        self.VEC,self.VEC_all=func.choose_and_run_make_kgrid(self.nkp,self.SYMM_OP,\
                                          self.b_vec,  self.clib_dir,self.clib_file)
        print ("...K-grid done.")
    def read_ene(self):
        print ("...I am reading energies...")
        self.ENE,self.nbnd = func.read_ene(self.liczba,self.file_bands)
        print ("...Done.")
    def ene_from_irreducible_to_whole_grid(self):
        print ("...I am expanding energies for whole grid...")
        self.ENE = func.ene_from_irreducible_to_whole_grid(
                 self.ENE,self.VEC,self.VEC_all,self.nkp)
        print ("...Done.")
    def complete_incomplete_bands(self):
        print ("...I am reconstructing incomplete conducting bands...")
        self.nband,self.ENE = func.complete_incomplete_bands(self.nbnd,self.ENE)
        print ("...Done.")
    def rearrange_data(self):
        print ("...I am rearranging data...")
        self.ENE = func.rearrange_data(self.nband,self.ENE)
        print ("...Done.")
    def which_bands_cross_fermi(self):
        self.ENE = func.which_bands_cross_fermi(self.ENE)
    def run_calc(self):
        print ("...I am running calculations...")
        func.run_calc(self.VEC,self.alat,self.COMMAND,self.file_in, self.file_in2)
        print ("...Calculations done.")
    def write_to_file(self):
        print ("...I am writting to file...")
        func.write_to_file(self.ENE,self.VEC,self.VEC_all,self.b_vec, self.alat, \
                           self.nkp,self.liczba,self.imE_to_lifetime_conv,\
                           self.file_frmsf,self.file_lifetime_frmsf,self.file_mayavi)
        print ("...Done")
####MAIN


####make k-grid
calc=init_calc()
obj=Obiekt()
#0. ask for input
obj.ask_for_input(calc.file_out)

#1. read data
obj.read_data()
#1.a.find reciprocal lattice vectors
obj.recip_vec_gen()
print( obj.b_vec)

#2.Find symmetry operations using lattice vectors from filei5
obj.set_sym_bl()

#3. Make k-grid
obj.make_k_grid()
print ('No of points=',len(obj.VEC_all), 'no of irreducible points=',len(obj.VEC))
#3.a.print to file
h=open('kvec.dat','w')
for ni,i in enumerate(obj.VEC_all):
 for j in i:
  h.write(str(j)+' ')
 h.write('\n')
h.close()

#4.run calc
obj.run_calc()
#4.a.write to file
h=open('ene.dat','w')
for i in obj.ENE:
  for k in i:
   h.write(str(k[0])+' ')
  h.write('\n')
h.close()

#5. read energy 
obj.read_ene()
obj.ene_from_irreducible_to_whole_grid()
obj.complete_incomplete_bands()

#6. rearrange data from ENE[k-point][nband] to ENE[nband][k-point]
obj.rearrange_data()
print( 'No of bands = ',len(obj.ENE),'. Each of them on the grid of'),
for i in obj.ENE: print(len(i)),
print( 'kpoints')

#7. check which bands cross EF
obj.which_bands_cross_fermi()
print( 'Number of pieces of Fermi surface=',len(obj.ENE))

#8. write results to file readable by fermisurfer and my program in mayavi
obj.write_to_file()



