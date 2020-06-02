
#include<iostream>
#include <vector>       
double eps1=1e-4;

using namespace std;


extern "C" 
{
void ruotaijk(vector<vector<double> > s,
			std::vector<double>  kp, std::vector<double> &kr)
 {
 for (int m=0;m<3;m++){
   kr[m]=s[m][0]*kp[0]+s[m][1]*kp[1]+s[m][2]*kp[2];
   }
 }

void check(int n,	
	   std::vector<std::vector<double> > &KVEC,
	   std::vector<int> &kw,	   
	   std::vector<int> &ieq,
	   vector<vector<vector<double> > > &SYM_OP, 	
	   int no_of_symm_op, 	
           int nkp)
  {
  int flag=1;
  int naux=0;
  vector<double> kr{0,0,0};
  for (int s=0; s<no_of_symm_op; s++){
     ruotaijk( SYM_OP[s],KVEC[n],kr ); 
     for (int j=0; j<3;j++){
      while (kr[j]>=nkp)   kr[j]=kr[j]-nkp;
      while (kr[j]<=-1)   kr[j]=kr[j]+nkp;
      }
      for (int npk=0;npk<n;npk++){
	 if ( (abs(kr[0]-KVEC[npk][0]))<eps1 and   (abs(kr[1]-KVEC[npk][1]))<eps1 and  (abs(kr[2]-KVEC[npk][2]))<eps1){
//	   if (kr[0]==KVEC[npk][0] and  kr[1]==KVEC[npk][1] and kr[2]==KVEC[npk][2]) {

	         kw[n]=-1;
		 naux=npk;
		 while (kw[naux]==-1)  naux=ieq[naux];
		 ieq[n]=naux;
		 kw[naux]=kw[naux]+1;
		 flag=0;
		 break;
	   }
         if (flag==0) break;
	 }
     if (flag==0) break;
     }
  }

int* FS_make_kgrid(int nkp,
		   int nso,
  		   double* SYMM_OP)
 {
//innitialization
 int no_of_kpoint=nkp*nkp*nkp;
 //double** KVEC=new double*[no_of_kpoint];
 std::vector<vector<double>> KVEC; //STORAGE OF KPOINTS
 KVEC.assign(no_of_kpoint,vector < double > (3, 0)); 
 //int* kw=new int[no_of_kpoint];
 std::vector<int> kw; //STORAGE OF WEIGHTS
 //int* ieq=new int[no_of_kpoint];
 std::vector<int> ieq; //STORAGE OF THE NO OF EQUIVALENT POINT
 for (int i=0;i<no_of_kpoint;i++){
  kw.push_back(1);
  ieq.push_back(0);
  }

 vector<vector<vector<double> > > SYMM_OP2; //STORAGE OF SYMM. OPERATION. SYMM_OP IS 1 dim ARRAY, WE ARE RESHAPING IT TO THE ARRAY[no_of_symmetry_operations][3][3]
 SYMM_OP2.assign(nso, vector < vector < double > >(3, vector < double >(3,0)));
 int n=0;
 for (int i=0; i<nso;i++){
  for (int j=0; j<3;j++){
   for (int k=0; k<3;k++){
    SYMM_OP2[i][j][k]=double(SYMM_OP[n]);
    n=n+1;
    }
   }
  }
//run
 n=0;
 for (int i=0;i<nkp;i++){
  for (int j=0;j<nkp;j++){
   for (int k=0;k<nkp;k++){
    KVEC[n][0]=double(i); KVEC[n][1]=double(j); KVEC[n][2]=double(k);
    check(n,KVEC,kw,ieq,SYMM_OP2,nso,nkp); 
    n=n+1;
    }
   }
  }
 cout<<"C ended\n";

 int* ieq2=new int[2*no_of_kpoint];
 for (int i=0;i<no_of_kpoint;i++){
  ieq2[i]=ieq[i];
  ieq2[no_of_kpoint+i]=kw[i];
  }
 return ieq2;

 }

}



