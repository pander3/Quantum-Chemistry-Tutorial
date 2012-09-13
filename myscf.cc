//http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming
//Projects 1-13
#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <sys/types.h>
#include <string>
#include <sys/resource.h>
#include <time.h>
#include <fstream>
#include "masses.h"

#define INDEX(i,j) (i>j) ? (ioff[i]+j) : (ioff[j]+i) //Define INDEX function
int ioff[10000];

namespace psi
{
  namespace my_scf
  {
    extern "C"
    {
      FILE *infile,*outfile;
      char *psi_file_prefix;
      const char *gprgid(void){ const char *prgid = "MY_SCF"; return (prgid); }
    }
  }
}// namespace psi::mn_scf

using namespace psi::my_scf;

double ****init_4d(int n);
void constructorthogonalization(double ***matrix, double ***eigenvectors,double **eigenvalues,int *sopi,int nirreps);
void orthogonalize(double ***final,double ***orthogonal,double ***initial, int *sopi,int nirreps);
void constructdensity(double ***matrix,double ***AOeigenvects,int *sopi,int irreps,int *filled);
void newfockmatrix(double ***matrix,double ***corehamiltonian,double ***density,double *tei,int *sopi,int nirreps,int *symmo,bool outofcore);
double rms3d(double ***oldmatrix,double ***newmatrix, int *sopi,int nirreps);
double calcSCFenergy(double ***density,double***corehamiltonian,double***fock,int *sopi,int nirreps);
void TEIarrayconverter(double *oned,double ****fourd,int n);
void buildteimo(double *arrAO,double *arrMO,double ***eigenvec,int *sopi,int nirreps,int *symmo,int nao, bool smart);
double calculateMP2(double *arrMO,double **eigenvalues,int *sopi,int nirreps,int *filled,int *symmo);
void buildtaus(double ****tt,double ****t,double **t1,double ****t2,int *m,int n,int *symmo);
void buildF(double **matrix,double ***sfock,double **t1,double ****so,double ****tt,int *sopi,int *filled,int nirreps,int nso,int *symmo);
void buildW(double ****matrix,double ****so,double **t1,double ****t2,double ****tau,int nso,int *filled,int *sopi,int nirreps,int *symmo);
void buildT1(double **matrix,double ***sfock,double **t1,double ****t2,double ****so,double **Fmat,int nso,int *filled,int *sopi,int nirreps,int *symmo);
void buildT2(double ****matrix,double ***sfock,double ****so,double **t1,double ****t2,double ****t,double **Fmat,double ****Wmat,int nso,int *filled,int *sopi,int nirreps,int *symmo);
double calcEcc(double ***sfock,double ****so,double **t1,double ****t2,int *filled,int *sopi,int nirreps,int *symmo);
double calcEt(double ***sfock,double **t1,double ****t2,double ****so,int *filled,int *sopi,int nirreps,int *symmo);
double rms4d(double ****oldmatrix, double ****newmatrix,int nso,bool flag);
void calcErrorSCF(double ***err,double ***fock,double ***density,double ***S,int *sopi,int nirreps);
void load3dsymm(double ****matrix4d,double ***matrix,int maxheight,int height,int *n,int m);
void buildBSCF(double **matrix,double ****errormatrix,int *sopi,int nirreps,int bsize);
double rms2d(double **oldmatrix,double **newmatrix,int n,bool flag);
void calcErrorCC(double **errorT1,double ****errorT2,double **t1,double **newt1, double****t2,double ****newt2,int n);
void buildBCC(double **matrix,double ***errorT1,double*****errorT2,int bs,int n);
void load4d(double *****matrix5d,double ****matrix,int maxheight,int height,int n);
void load2d(double ***matrix3d,double **matrix,int maxheight,int height,int nbu);
void modifyT(double **matrix1,double ****matrix2,double ***it1,double *****it2,int nso, int bsize, double *right);
void bubblesort(double *array,int n);
void checkmp2(double ****t2,double ****so,int m,int *n,int *filled,int *symmo, double Emp2);
void buildinitialT(double **t1,double ****t2,double ****so,double **evals,int nso,int nirreps,int *sopi,int *filled,int *symmo);
void buildsofock(double ***sfock,double **evals,int nirreps,int *sopi);
void buildHCIS(double **matrix,double ***sfock,double ****so,int nirreps,int *sopi,int *filled, int *symmo);
void calcsinglets(double **singlet,double ***sfock,double *tei,int nirreps,int *sopi,int *filled,int *symmo,int occupied,int unoccupied, double *eval,double **evec);
void calctriplets(double **triplet,double ***sfock,double *tei,int nirreps,int *sopi,int *filled,int *symmo,int occupied,int unoccupied, double *eval,double **evec);
void buildA(double **matrix,double **cis,int size);
void buildB(double **matrix,double ****so,int nirreps,int *sopi,int *filled,int *symmo);
void buildTDHF(double **matrix,double **A,double **B,int size);
void setupTDHF(double **A,double **B,double **AplusB,double **AminusB,double **ans,int occupied,int unoccupied);
void comparison(double **ab,double **comp, int initialL,int N,int L,int *location,double *locVal);
void gschmidt(double **v,int rows,int columns);
void proj(double *u,double *v,int n,double *rtrn);
void davidsonlui(int occupied,int unoccupied,double **singlets,int initialL,double *dlvals,int high);
void clockPrint(std::string title);
void findfilled(double *zval,double **e_vals,int natom,int nirreps,int *sopi,int &occupied,int &unoccupied,int *filled);
double dist(double x1, double x2, double y1, double y2, double z1, double z2);
void moleculedata(int atomno,int *zno,double **geom);


extern "C" void dgeev_( char* jobvl,char* jobvr,int* n, double* a,int* lda,double* wr,double* wi,double* vl, int* ldvl,double* vr, int* ldvr, double* work, int* lwork, int* info);


int main(int argc, char *argv[])
{  
  ioff[0]=0;
  for(int i=1;i<10000;i++)
  {
    ioff[i] = ioff[i-1] + i;
  }

  psi_start(&infile,&outfile,&psi_file_prefix,argc-1,argv+1,0);
  ip_cwk_add(":MY_SCF");
  psio_init(); psio_ipv1_config();
  tstart(outfile);

    int nao,nso;
    int natom;
    double *zval;
    int *zv;
    double **geom;
    int *filled;
    double enuc;
    double Escf,Emp2,Ecc,Et;    
    double rms2,rms3,rms4;
    double temp;
    int occupied=0,unoccupied=0;

    double *scratch,***S,***T,***V,***H;
    double *TEI;
    double **e_vals,***e_vecs;
    double ***shalf,***fock,***C,***density,***newfock,***orthofock,***olddensity;
    double *TEI_MO;
    double ****SOarray;        
    double **initialT1,****initialT2,**T1,****T2,**newT1,****newT2;
    double ***sofock;
    double **F,****W,**newF,****newW;
    double ****tau,****tautilde,****newtau,****newtautilde;    
    double ***error,****iterror,****itfock,**Bmat,**B;   
    int count=0,errno=8;    
    int nirreps,*sopi,*symm_offset;

    bool outofcore=false;

    double ***itT1;
    double *****itT2;
    double **T1err;
    double ****T2err;
    double **tempT1,**T1new;
    double ****tempT2,****T2new;
    double ***itT1err,*****itT2err;

    double **HCIS;
    double *CISevals,**CISevecs;
    double **singlets,*evalsinglets,**evecsinglets;
    double**triplets,*evaltriplets,**evectriplets;
    double **TDHF;
    double **A;
    double **AplusB;
    double **AminusB;
    double *val,**vec;
    double **ans;
    double *real,*imag;
    double **workspace;

    int ntri;
    struct iwlbuf InBuf;
    Value *valptr;
    Label *lblptr;
    int lastbuf=0;
    int idx;

    chkpt_init(PSIO_OPEN_OLD);
    nao = chkpt_rd_nao(); //number of atomic orbitals
    nso = 2*nao;  //number of spin orbitals
    enuc = chkpt_rd_enuc();  //nuclear repulsion energy
    nirreps = chkpt_rd_nirreps();  //number of irreducible representations
    sopi = (int*) malloc(nirreps*sizeof(int));  //Number of Symmetric Orbitals per irreducible representation
    sopi = chkpt_rd_sopi();  
    natom = chkpt_rd_natom();  //number of atoms
    zval = (double *) malloc(natom*sizeof(double));  //atom identities
    zv = (int *) malloc(natom*sizeof(int));
    zval = chkpt_rd_zvals();
    for(int i=0;i<natom;i++)
    {
      //Cast the zvals to integers
      zv[i]=zval[i];
    }
    geom=block_matrix(natom,3);  //atomic geometry
    geom = chkpt_rd_geom();
    ntri = nao*(nao+1)/2;  //size of two-electron-integral (tei)
    chkpt_close();	

    //Allocate memory
    filled = (int*) malloc(nirreps*sizeof(int)); //filled orbitals by irrep
    scratch = init_array(ntri);

    H = (double ***) malloc(nirreps*sizeof(double**)); //Core hamiltonian
    shalf = (double ***) malloc(nirreps*sizeof(double**));//S^-1/2 orthogonalization matrix
    fock = (double ***) malloc(nirreps*sizeof(double**)); //Fock matrix
    e_vecs = (double ***) malloc(nirreps*sizeof(double**)); //Eigenvector matrix
    C = (double ***) malloc(nirreps*sizeof(double**)); //Eigenvectors in original non-AO basis
    density = (double ***) malloc(nirreps*sizeof(double**)); //density matrix
    olddensity = (double ***) malloc(nirreps*sizeof(double**)); //Previous iteration's density matrix
    newfock = (double ***) malloc(nirreps*sizeof(double**)); //Iterated Fock matrix
    orthofock = (double ***) malloc(nirreps*sizeof(double**));//Orthogonalized Fock matrix
    e_vals = (double **) malloc(nirreps*sizeof(double*));

    TEI = init_array(ntri*(ntri+1)/2); //Two-Electron integrals	
    TEI_MO = init_array(ntri*(ntri+1)/2); //TEI in Molecular Orbital Basis   
    error = (double ***) malloc(nirreps*sizeof(double**));
    iterror = ( double **** ) malloc(nirreps*sizeof(double***));
    itfock = ( double **** ) malloc(nirreps*sizeof(double***));
    sofock = (double***) malloc(nirreps*sizeof(double**)); //Spin orbital fock matrix

    for(int h=0;h<nirreps;h++)
    {
      H[h] = block_matrix(sopi[h],sopi[h]);
      shalf[h] = block_matrix(sopi[h],sopi[h]);
      fock[h] = block_matrix(sopi[h],sopi[h]);
      e_vecs[h] = block_matrix(sopi[h],sopi[h]);
      C[h] = block_matrix(sopi[h],sopi[h]);
      density[h] = block_matrix(sopi[h],sopi[h]);
      olddensity[h] = block_matrix(sopi[h],sopi[h]);
      newfock[h] = block_matrix(sopi[h],sopi[h]);
      orthofock[h] = block_matrix(sopi[h],sopi[h]);
      e_vals[h] = init_array(sopi[h]);

      error[h] = block_matrix(sopi[h],sopi[h]);
      iterror[h] = (double***) malloc(errno*sizeof(double**));
      itfock[h] = (double***) malloc(errno*sizeof(double**));
      sofock[h] = block_matrix(2*sopi[h],2*sopi[h]);
      for(int j=0;j<errno;j++)
      {
        iterror[h][j] = block_matrix(sopi[h],sopi[h]);
        itfock[h][j] = block_matrix(sopi[h],sopi[h]);
      }
    }       

    SOarray = (double ****) init_4d(nso); //TEI in a 4D array (spin orbit array)
    initialT1 = init_matrix(nso,nso);  //Initial guess of T1 cluster amplitudes
    initialT2 = (double ****) init_4d(nso); //Initial guess of T2 cluster amplitudes
    T1 = init_matrix(nso,nso);  //T1 Cluster amplitudes
    T2 = (double ****) init_4d(nso);  //T2 Cluster amplitudes
    newT1 = init_matrix(nso,nso);
    newT2 = (double ****) init_4d(nso);
    F = init_matrix(nso,nso); //initial One Particle intermediate for CCSD
    newF = init_matrix(nso,nso); //iterative One Particle intermediate for CCSD
    W = (double ****) init_4d(nso); //initial Two Particle intermediate for CCSD
    newW = (double ****) init_4d(nso); //iterative Two Particle intermediate for CCSD
    tau = (double ****) init_4d(nso); 
    tautilde = (double ****) init_4d(nso);
    newtau = (double ****) init_4d(nso);
    newtautilde = (double ****) init_4d(nso);

    itT1 = (double ***) malloc(errno*sizeof(double**));
    itT2 = (double *****) malloc(errno*sizeof(double ****));
    T1err = init_matrix(nso,nso);
    T2err = (double****) init_4d(nso);
    tempT1 = init_matrix(nso,nso);
    tempT2 = (double ****) init_4d(nso);
    T1new = init_matrix(nso,nso);
    T2new = (double ****) init_4d(nso);
    itT1err = (double ***) malloc(errno*sizeof(double**));
    itT2err = (double *****) malloc(errno*sizeof(double****));
    for(int i=0;i<errno;i++)
    {
      itT1[i] = (double **) malloc(nso*sizeof(double*));
      itT2[i] = (double ****) malloc(nso*sizeof(double***));
      itT1err[i] = (double **) malloc(nso*sizeof(double*));
      itT2err[i] = (double ****) malloc(nso*sizeof(double ***));
      for(int j=0;j<nso;j++)
      {
        itT1[i][j]=(double *) malloc(nso*sizeof(double));
        itT2[i][j]=(double ***) malloc(nso*sizeof(double**));
        itT1err[i][j]=(double*) malloc(nso*sizeof(double));
        itT2err[i][j]=(double***) malloc(nso*sizeof(double**));
        for(int k=0;k<nso;k++)
        {
          itT2[i][j][k]=(double**) malloc(nso*sizeof(double*));
          itT2err[i][j][k]=(double**) malloc(nso*sizeof(double*));
          for(int l=0;l<nso;l++)
          {
            itT2[i][j][k][l]=(double*) malloc(nso*sizeof(double));
            itT2err[i][j][k][l]=(double*) malloc(nso*sizeof(double));
          }
        }
      }
    }

    //Obtain One-Electron-Integrals
    symm_offset = init_int_array(nirreps);
    for(int h=1;h<nirreps;h++)
    {
      symm_offset[h] = symm_offset[h-1]+sopi[h-1];
    }
    S = (double ***) malloc(nirreps*sizeof(double**));
    iwl_rdone(PSIF_OEI,PSIF_SO_S,scratch,ntri,0,0,outfile);
    for(int h=0; h<nirreps;h++)
    {
      S[h] = block_matrix(sopi[h],sopi[h]);
      for(int i=0;i<sopi[h];i++)
      {
        for(int j=0;j<=i;j++)
        {
          int ij;
          ij = INDEX(i+symm_offset[h],j+symm_offset[h]);
          S[h][i][j] = S[h][j][i] = scratch[ij];
        }
      }
    }
    T = (double ***) malloc(nirreps*sizeof(double**));
    iwl_rdone(PSIF_OEI,PSIF_SO_T,scratch,ntri,0,0,outfile);
    for(int h=0;h<nirreps;h++)
    {
      T[h] = block_matrix(sopi[h],sopi[h]);
      for(int i=0;i<sopi[h];i++)
      {
        for(int j=0;j<=i;j++)
        {
          int ij;
          ij = INDEX(i+symm_offset[h],j+symm_offset[h]);
          T[h][i][j] = T[h][j][i] = scratch[ij];
        }
      }
    }
    V = (double ***) malloc(nirreps*sizeof(double**));
    iwl_rdone(PSIF_OEI,PSIF_SO_V,scratch,ntri,0,0,outfile);
    for(int h=0;h<nirreps;h++)
    {
      V[h] = block_matrix(sopi[h],sopi[h]);
      for(int i=0;i<sopi[h];i++)
      {
        for(int j=0;j<=i;j++)
        {
          int ij;
          ij = INDEX(i+symm_offset[h],j+symm_offset[h]);
          V[h][i][j] = V[h][j][i] = scratch[ij];
        }
      }
    }	
	
    //Obtain Two-Electron-Integrals
    if(outofcore==false)
    {
      iwl_buf_init(&InBuf,PSIF_SO_TEI,1e-14,1,0);
      do
      {
        iwl_buf_fetch(&InBuf);
        lblptr = InBuf.labels;
        valptr = InBuf.values;
        lastbuf = InBuf.lastbuf;
        for(idx=4*InBuf.idx;InBuf.idx<InBuf.inbuf;InBuf.idx++)
        {
          int i,j,k,l,ij,kl,ijkl;
	  i = abs( (int) lblptr[idx++]);
	  j = (int) lblptr[idx++];
	  k = (int) lblptr[idx++];
	  l = (int) lblptr[idx++];
	  ij = INDEX(i,j);
       	  kl = INDEX(k,l);
	  ijkl = INDEX(ij,kl);
	  TEI[ijkl] = (double) valptr[InBuf.idx];
        }
      }while( !lastbuf );
      iwl_buf_close(&InBuf,1);
    }


    moleculedata(natom,zv,geom);


    //Calculate the core hamiltonian
    for(int h=0;h<nirreps;h++)
    {
      for(int i=0;i<sopi[h];i++)
      {
        for(int j=0;j<sopi[h];j++)
        {
          if(sopi[h] != 0 )
          {
            H[h][i][j]=T[h][i][j]+V[h][i][j];
          }
        }
      }
    }

    //Construct the S^-1/2 Orthogonalization Matrix
    for(int h=0;h<nirreps;h++)
    {
      if(sopi[h] != 0)
      {
        sq_rsp(sopi[h],sopi[h],S[h],e_vals[h],1,e_vecs[h],1e-13);
      }
    }
    constructorthogonalization(shalf,e_vecs,e_vals,sopi,nirreps);

    //Construct the inital guess Fock Matrix and diagonalize
    orthogonalize(fock,shalf,H,sopi,nirreps);
    for(int h=0;h<nirreps;h++)
    {
      if(sopi[h] != 0 )
      {
        sq_rsp(sopi[h],sopi[h],fock[h],e_vals[h],1,e_vecs[h],1e-13);
      }
    }
    findfilled(zval,e_vals,natom,nirreps,sopi,occupied,unoccupied,filled);

    HCIS = block_matrix(occupied*unoccupied,occupied*unoccupied);
    CISevals = init_array(occupied*unoccupied);
    CISevecs = block_matrix(occupied*unoccupied,occupied*unoccupied);
    singlets = block_matrix( (occupied/2)*(unoccupied/2),(occupied/2)*(unoccupied/2) );
    evalsinglets = init_array((occupied/2)*(unoccupied/2) );
    evecsinglets = block_matrix( (occupied/2)*(unoccupied/2),(occupied/2)*(unoccupied/2) );
    triplets=block_matrix( (occupied/2)*(unoccupied/2),(occupied/2)*(unoccupied/2) );
    evaltriplets=init_array( (occupied/2)*(unoccupied/2) );
    evectriplets=block_matrix( (occupied/2)*(unoccupied/2),(occupied/2)*(unoccupied/2));
    TDHF=block_matrix(2*occupied*unoccupied,2*occupied*unoccupied);
    A=block_matrix(occupied*unoccupied,occupied*unoccupied);
    B=block_matrix(occupied*unoccupied,occupied*unoccupied);
    AminusB = block_matrix(occupied*unoccupied,occupied*unoccupied);
    AplusB = block_matrix(occupied*unoccupied,occupied*unoccupied);
    val=init_array(occupied*unoccupied);
    vec=block_matrix(occupied*unoccupied,occupied*unoccupied);
    ans=block_matrix(occupied*unoccupied,occupied*unoccupied);

    //Transform the eigenvectors into the original (non-orthogonal) AO basis 
    //and construct initial density matrix
    for(int h=0;h<nirreps;h++)
    {
      mmult(shalf[h],0,e_vecs[h],0,C[h],0,sopi[h],sopi[h],sopi[h],0);
    } 
    constructdensity(density,C,sopi,nirreps,filled);

    //Compute Self-Consistent-Field electronic energy
    Escf = calcSCFenergy(density,H,H,sopi,nirreps)+enuc;

    int bsize=0; count = 0;
    //Iterate
    do
    {
      temp = Escf;      

      newfockmatrix(newfock,H,density,TEI,sopi,nirreps,symm_offset,outofcore); 

      calcErrorSCF(error,newfock,density,S,sopi,nirreps);
      load3dsymm(iterror,error,errno,count%errno,sopi,nirreps);
      load3dsymm(itfock,newfock,errno,count%errno,sopi,nirreps);

      if(count<errno)
      {
        bsize = count+2;
      }
      else
      {
        bsize = errno+1;
      }

      Bmat = block_matrix(bsize,bsize);
      buildBSCF(Bmat,iterror,sopi,nirreps,bsize);
      int *IPIV;
      IPIV = (int *) malloc(bsize*sizeof(int));
      double *rhs;
      rhs = init_array(bsize);
      rhs[bsize-1] = -1.0;     

      C_DGESV(bsize,1,Bmat[0],bsize,IPIV,rhs,bsize);

      if(count != 0)
      {
        for(int h=0;h<nirreps;h++)
        {
          for(int i=0;i<sopi[h];i++)
          {
            for(int j=0;j<sopi[h];j++)
            {
              newfock[h][i][j]=0.0;
              for(int k=0;k<bsize-1;k++)
              {
                newfock[h][i][j] += itfock[h][k][i][j] * rhs[k]; 
              }
            }
          }
        }     
      }
      free(IPIV);
      free(rhs);      
      free_block(Bmat);
      
      orthogonalize(orthofock,shalf,newfock,sopi,nirreps);
      for(int h=0;h<nirreps;h++)
      {
        if(sopi[h] != 0)
        {
          sq_rsp(sopi[h],sopi[h],orthofock[h],e_vals[h],1,e_vecs[h],1e-13);
        }
      }

      //Back Transform and construct density matrix
      for(int h=0;h<nirreps;h++)
      {
        for(int i=0;i<sopi[h];i++)
        {
          for(int j=0;j<sopi[h];j++)
          {
            C[h][i][j]=0.0;
            density[h][i][j]=0.0;
          }
        }
      }
      for(int h=0;h<nirreps;h++)
      {
        mmult(shalf[h],0,e_vecs[h],0,C[h],0,sopi[h],sopi[h],sopi[h],0);
      }
      constructdensity(density,C,sopi,nirreps,filled);

      Escf = calcSCFenergy(density,H,newfock,sopi,nirreps) + enuc;
      rms3 = rms3d(olddensity,density,sopi,nirreps);      
      count++;
    }while( fabs(temp - Escf)>1.0e-14 || rms3 > 1.0e-12 );
    std::cout << "SCF Energy: "<< std::setprecision(14) <<Escf << std::endl;

    if(outofcore==false)
    { 
      buildsofock(sofock,e_vals,nirreps,sopi);
      buildteimo(TEI,TEI_MO,C,sopi,nirreps,symm_offset,nao,1);
      TEIarrayconverter(TEI_MO,SOarray,nso);

      //Calculate Energy of MP2 Method using spin orbital expression
      Emp2 = calculateMP2(TEI_MO,e_vals,sopi,nirreps,filled,symm_offset);
      std::cout << "MP2 Energy: " << std::setprecision(12) << Emp2 << std::endl;


      //Check to see if the energy of the MP2 methods match
      buildinitialT(initialT1,initialT2,SOarray,e_vals,nso,nirreps,sopi,filled,symm_offset);
      checkmp2(initialT2,SOarray,nirreps,sopi,filled,symm_offset,Emp2);
    
      //Build Initial Matrices for CCSD
      buildtaus(tautilde,tau,initialT1,initialT2,sopi,nirreps,symm_offset);
      buildF(F,sofock,initialT1,SOarray,tautilde,sopi,filled,nirreps,nso,symm_offset);
      buildW(W,SOarray,initialT1,initialT2,tau,nso,filled,sopi,nirreps,symm_offset);
      buildT1(T1,sofock,initialT1,initialT2,SOarray,F,nso,filled,sopi,nirreps,symm_offset);
      buildT2(T2,sofock,SOarray,initialT1,initialT2,tau,F,W,nso,filled,sopi,nirreps,symm_offset);
      Ecc = calcEcc(sofock,SOarray,T1,T2,filled,sopi,nirreps,symm_offset);

      bsize=0;
      count=0;
      do
      {
        temp = Ecc;
        if(count<errno)
        {
          bsize=count+2;
        }
        else
        {
          bsize=errno+1;
        }      
      
        if(count==0)
        { 
          buildtaus(newtautilde,newtau,T1,T2,sopi,nirreps,symm_offset);          
          buildF(newF,sofock,T1,SOarray,newtautilde,sopi,filled,nirreps,nso,symm_offset);
          buildW(newW,SOarray,T1,T2,newtau,nso,filled,sopi,nirreps,symm_offset);              
          buildT1(newT1,sofock,T1,T2,SOarray,newF,nso,filled,sopi,nirreps,symm_offset);     
          buildT2(newT2,sofock,SOarray,T1,T2,newtau,newF,newW,nso,filled,sopi,nirreps,symm_offset);
          calcErrorCC(T1err,T2err,T1,newT1,T2,newT2,nso);
          load2d(itT1,newT1,errno,count,nso);      
          load2d(itT1err,T1err,errno,count,nso);
          load4d(itT2,newT2,errno,count,nso);
          load4d(itT2err,T2err,errno,count,nso);
        }
        else 
        {
          buildtaus(newtautilde,newtau,newT1,newT2,sopi,nirreps,symm_offset);          
          buildF(newF,sofock,newT1,SOarray,newtautilde,sopi,filled,nirreps,nso,symm_offset);
          buildW(newW,SOarray,newT1,newT2,newtau,nso,filled,sopi,nirreps,symm_offset);              
          buildT1(T1new,sofock,newT1,newT2,SOarray,newF,nso,filled,sopi,nirreps,symm_offset);
          buildT2(T2new,sofock,SOarray,newT1,newT2,newtau,newF,newW,nso,filled,sopi,nirreps,symm_offset);
          calcErrorCC(T1err,T2err,T1new,newT1,T2new,newT2,nso);
          load2d(itT1,T1new,errno,count%errno,nso);
          load2d(itT1err,T1err,errno,count%errno,nso);
          load4d(itT2,T2new,errno,count%errno,nso);
          load4d(itT2err,T2err,errno,count%errno,nso);
     
          Bmat=block_matrix(bsize,bsize);
          buildBCC(Bmat,itT1err,itT2err,bsize,nso); 
          int *IPIV;
          IPIV = (int *) malloc(bsize*sizeof(int));
          double *rhs;
          rhs = init_array(bsize);
          rhs[bsize-1] = -1.0;     

          C_DGESV(bsize,1,Bmat[0],bsize,IPIV,rhs,bsize);
    
          modifyT(tempT1,tempT2,itT1,itT2,nso,bsize,rhs); 
          buildtaus(newtautilde,newtau,tempT1,tempT2,sopi,nirreps,symm_offset);
          buildF(newF,sofock,tempT1,SOarray,newtautilde,sopi,filled,nirreps,nso,symm_offset);
          buildW(newW,SOarray,tempT1,tempT2,newtau,nso,filled,sopi,nirreps,symm_offset);
          buildT1(newT1,sofock,tempT1,tempT2,SOarray,newF,nso,filled,sopi,nirreps,symm_offset);
          buildT2(newT2,sofock,SOarray,tempT1,tempT2,newtau,newF,newW,nso,filled,sopi,nirreps,symm_offset);      
        
          free(IPIV);
          free(rhs);      
          free_block(Bmat);
        }
      
        if(count==0)
        {          
          rms2 = rms2d(T1,newT1,nso,0);
          rms4 = rms4d(T2,newT2,nso,0);
          Ecc = calcEcc(sofock,SOarray,newT1,newT2,filled,sopi,nirreps,symm_offset);
        }
        else
        {
          rms2 = rms2d(newT1,T1new,nso,0);
          rms4 = rms4d(newT2,T2new,nso,0);
          Ecc = calcEcc(sofock,SOarray,newT1,newT2,filled,sopi,nirreps,symm_offset);
        }
        count++;
      }while( fabs( temp - Ecc )>1e-12 || rms4 > 1e-12 || rms2 > 1e-12);
      std::cout << "CCSD Energy: " <<std::setprecision(12)<< Ecc << std::endl;

      //Calculate the Triples Correction
      Et = calcEt(sofock,newT1,newT2,SOarray,filled,sopi,nirreps,symm_offset);
      std::cout << "CCSD(T) Energy: " << std::setprecision(12)<<Et << std::endl;
    

      //Calculate Configure Interaction Singles (CIS) Hamiltonian
      zero_mat(HCIS,occupied*unoccupied,occupied*unoccupied);
      zero_mat(CISevecs,occupied*unoccupied,occupied*unoccupied);
      buildHCIS(HCIS,sofock,SOarray,nirreps,sopi,filled,symm_offset);
      sq_rsp(occupied*unoccupied,occupied*unoccupied,HCIS,CISevals,1,CISevecs,1e-13);
 
      calcsinglets(singlets,sofock,TEI_MO,nirreps,sopi,filled,symm_offset,occupied,unoccupied,evalsinglets,evecsinglets);
      calctriplets(triplets,sofock,TEI_MO,nirreps,sopi,filled,symm_offset,occupied,unoccupied,evaltriplets,evectriplets);


      //Construct Time-Dependent Hartree Fock Hamiltonian
      buildA(A,HCIS,occupied*unoccupied);
      buildB(B,SOarray,nirreps,sopi,filled,symm_offset);
      buildTDHF(TDHF,A,B,occupied*unoccupied);


      int sz1=1,sz=occupied*unoccupied;
      int worksize=4*sz;
      char decider='n';
      double lefteigen[0][0],righteigen[0][0];
      int info;
 
      real=init_array(sz);
      imag=init_array(sz);
      workspace=block_matrix(4*sz,4*sz);

      setupTDHF(A,B,AplusB,AminusB,ans,occupied,unoccupied);
      dgeev_(&decider,&decider,&sz,&ans[0][0],&sz,&real[0],&imag[0],&lefteigen[0][0],&sz1,&righteigen[0][0],&sz1,&workspace[0][0],&worksize,&info);
      bubblesort(real,sz);
      for(int i=0;i<sz;i++)
      {
        real[i] = sqrt(real[i]);
      }
      std::cout << std::endl;
      std::cout << "Time Dependent Hartree Fock Calculation Complete"<<std::endl;

      std::cout << "Singlet Excited States" << std::endl;
      for(int i=0;i<(occupied/2)*(unoccupied/2);i++)
      {
        std::cout << std::setprecision(6) << evalsinglets[i] << std::endl;
      }
      std::cout << std::endl;
      std::cout << "Triplet Excited States" << std::endl;
      for(int i=0;i<(occupied/2)*(unoccupied/2);i++)
      {
        std::cout << std::setprecision(6) << evaltriplets[i] << std::endl;
      }
  
      int dlno=1;
      double *dlvals;
      dlvals=init_array(dlno);
      davidsonlui(occupied,unoccupied,singlets,dlno,dlvals,1);
      std::cout << std::endl;
      std::cout << "Davidson-Lui Singlet Eigenvalues:"<<std::endl;
      for(int i=0;i<dlno;i++)
      {
        std::cout << dlvals[i] << std::endl;
      }
      std::cout << "Davidson-Lui Triplet Eigenvalues:"<<std::endl;
      davidsonlui(occupied,unoccupied,triplets,dlno,dlvals,1);
      for(int i=0;i<dlno;i++)
      {
        std::cout << dlvals[i] << std::endl;
      }
    }

    free(scratch);
    free(sopi);
    for(int i=0;i<nirreps;i++)
    {
      free_block(S[i]);
      free_block(V[i]);
      free_block(T[i]);
      free_block(H[i]);
      free_block(shalf[i]);
      free_block(fock[i]);
      free_block(e_vecs[i]);
      free_block(C[i]);
      free_block(density[i]);
      free_block(olddensity[i]);
      free_block(newfock[i]);
      free_block(orthofock[i]);
      free(e_vals[i]);
    }
    free(S);
    free(T);
    free(V);
    free(H);
    free(shalf);
    free(fock);
    free(e_vecs);
    free(C);
    free(density);
    free(olddensity);
    free(newfock);
    free(orthofock);
    free(e_vals);
    free(TEI);
    free(TEI_MO);
    free(symm_offset);
    for(int i=0;i<nirreps;i++)
    {
      for(int j=0;j<errno;j++)
      {
        free_block(iterror[i][j]);
        free_block(itfock[i][j]);
      }
      free_block(error[i]);
      free(iterror[i]);
      free(itfock[i]);
      free_block(sofock[i]);
    }
    free(error);
    free(iterror);
    free(itfock);
    free(sofock);
        
    free_matrix(initialT1,nso);
    free_matrix(T1,nso);
    free_matrix(newT1,nso);
    free_matrix(F,nso);
    free_matrix(newF,nso);
        
    for(int i=0;i<nso;i++)
    {
      for(int j=0;j<nso;j++)
      {
        for(int k=0;k<nso;k++)
        {
          free(SOarray[i][j][k]);
          free(initialT2[i][j][k]);
          free(T2[i][j][k]);
          free(newT2[i][j][k]);
          free(W[i][j][k]);
          free(newW[i][j][k]);
          free(tau[i][j][k]);
          free(newtau[i][j][k]);
          free(tautilde[i][j][k]);
          free(newtautilde[i][j][k]);
        }
        free(SOarray[i][j]);
        free(initialT2[i][j]);
        free(T2[i][j]);
        free(newT2[i][j]);
        free(W[i][j]);
        free(newW[i][j]);
        free(tau[i][j]);
        free(newtau[i][j]);
        free(tautilde[i][j]);
        free(newtautilde[i][j]);
      }
      free(SOarray[i]);
      free(initialT2[i]);
      free(T2[i]);
      free(newT2[i]);
      free(W[i]);
      free(newW[i]);
      free(tau[i]);
      free(newtau[i]);
      free(tautilde[i]);
      free(newtautilde[i]);
    }
    free(SOarray);
    free(initialT2);
    free(T2);
    free(newT2);
    free(W);
    free(newW);
    free(tau);
    free(newtau);
    free(tautilde);
    free(newtautilde);
     
    free_matrix(T1err,nso);
    free_matrix(tempT1,nso);
    free_matrix(T1new,nso);
    for(int i=0;i<nso;i++)
    {
      for(int j=0;j<nso;j++)
      {
        for(int k=0;k<nso;k++)
        {
          free(T2err[i][j][k]);
          free(tempT2[i][j][k]);
          free(T2new[i][j][k]);
        }
        free(T2err[i][j]);
        free(tempT2[i][j]);
        free(T2new[i][j]);
      }
      free(T2err[i]);
      free(tempT2[i]);
      free(T2new[i]);
    }
    free(T2err);
    free(tempT2);
    free(T2new);
    for(int i=0;i<errno;i++)
    {
      for(int j=0;j<nso;j++)
      {
        free(itT1[i][j]);
        free(itT1err[i][j]);
      }
      free(itT1[i]);
      free(itT1err[i]);
    }
    free(itT1);
    free(itT1err);
    for(int i=0;i<errno;i++)
    {
      for(int j=0;j<nso;j++)
      {
        for(int k=0;k<nso;k++)
        {
          for(int l=0;l<nso;l++)
          {
            free(itT2[i][j][k][l]);
            free(itT2err[i][j][k][l]);
          }
          free(itT2[i][j][k]);
          free(itT2err[i][j][k]);
        }
        free(itT2[i][j]);
        free(itT2err[i][j]);
      }
      free(itT2[i]);
      free(itT2err[i]);
    }
    free(itT2);
    free(itT2err);

    free_block(HCIS);
    free(CISevals);
    free_block(CISevecs);
    free_block(singlets);
    free(evalsinglets);
    free_block(evecsinglets);
    free_block(triplets);
    free(evaltriplets);
    free_block(evectriplets);
    free_block(TDHF);
    free_block(A);
    free_block(B);
    free_block(AminusB);
    free_block(AplusB);
    free(val);
    free_block(vec);
    free_block(ans);
    if(outofcore==false)
    {
      free_block(workspace);
      free(real);
      free(imag);
    }
  tstop(outfile);  
  psio_done();
  psi_stop(infile, outfile, psi_file_prefix);
 
  return 0;
}

//Initialize and zero a 4-d array of doubles
double ****init_4d(int n)
{
  double ****array=NULL;
  array = (double****) malloc(n*sizeof(double***));
  for(int i=0;i<n;i++)
  {
    array[i] = (double ***) malloc(n*sizeof(double**));
    for(int j=0;j<n;j++)
    {
      array[i][j] = (double **) malloc(n*sizeof(double*));
      for(int k=0;k<n;k++)
      {
        array[i][j][k] = (double *) malloc(n*sizeof(double));
        for(int l=0;l<n;l++)
        {
          array[i][j][k][l]=0.0;
        }
      }
    }
  }
  return(array);
}

//Construct S^-1/2 Orthogonalization Matrix
void constructorthogonalization(double ***matrix, double ***eigenvectors,double **eigenvalues,int *sopi,int nirreps)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3

  double ***invrooteval;
  double ***t;
  invrooteval=(double ***) malloc(nirreps*sizeof(double**)) ;
  t = (double ***) malloc(nirreps*sizeof(double**));  
  for(int h=0;h<nirreps;h++)
  {
    invrooteval[h] = block_matrix(sopi[h],sopi[h]);
    t[h] = block_matrix(sopi[h],sopi[h]);
    for(int i=0;i<sopi[h];i++)
    {
      for(int j=0;j<sopi[h];j++)
      { 
        invrooteval[h][i][j]=0.0;
        if( i==j )
        {
          if( eigenvalues[h][i] != 0 )
          {
            invrooteval[h][i][j] = sqrt(1/eigenvalues[h][i]);
          }
        }
      }
    }
  }
  for(int h=0;h<nirreps;h++)
  {
    mmult(eigenvectors[h],0,invrooteval[h],0,t[h],0,sopi[h],sopi[h],sopi[h],0);
    mmult(t[h],0,eigenvectors[h],1,matrix[h],0,sopi[h],sopi[h],sopi[h],0);
  }
  for(int h=0;h<nirreps;h++)
  {
    free_block(invrooteval[h]);
    free_block(t[h]);
  }

  free(invrooteval);
  free(t);
}

//Orthogonalize initial using orthogonal to form final
void orthogonalize(double ***final,double ***orthogonal,double ***initial, int *sopi,int nirreps)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3

  double *** t;
  t = (double ***) malloc(nirreps*sizeof(double**));
  for(int h=0;h<nirreps;h++)
  {
    t[h] = block_matrix(sopi[h],sopi[h]);
  } 

  for(int h=0;h<nirreps;h++)
  {
    for(int i=0;i<sopi[h];i++)
    {
      for(int j=0;j<sopi[h];j++)
      {
        final[h][i][j] = 0.0;
      }
    }
  }  
  for(int h=0;h<nirreps;h++)
  {
   mmult(orthogonal[h],1,initial[h],0,t[h],0,sopi[h],sopi[h],sopi[h],0);
   mmult(t[h],0,orthogonal[h],0,final[h],0,sopi[h],sopi[h],sopi[h],0);
  }
  for(int h=0;h<nirreps;h++)
  {
    free_block(t[h]);
  }
  free(t);
}

//Construct density matrix
void constructdensity(double ***matrix,double ***AOeigenvects,int *sopi,int irreps,int *filled)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3

  for(int h=0;h<irreps;h++)
  {
    for(int i=0;i<sopi[h];i++)
    {
      for(int j=0;j<sopi[h];j++)
      {
        matrix[h][i][j] = 0.0;
        for(int k=0;k<filled[h];k++)
        {
          matrix[h][i][j] += AOeigenvects[h][i][k]*AOeigenvects[h][j][k];
        }
      }
    }
  } 	
}
//Calculate Self-Consistent Field Electronic Energy
double calcSCFenergy(double ***density,double***corehamiltonian,double***fock,int *sopi,int nirreps)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3

  double energy = 0.0;
  double ***t;
  t = (double ***) malloc(nirreps*sizeof(double**)); 
  for(int h=0;h<nirreps;h++)
  {
    t[h]=block_matrix(sopi[h],sopi[h]);
  }
  for(int h=0;h<nirreps;h++)
  {
    for(int i=0;i<sopi[h];i++)
    {
      for(int j=0;j<sopi[h];j++)
      {
        t[h][i][j] = corehamiltonian[h][i][j] + fock[h][i][j];
        energy += density[h][i][j] * t[h][i][j];
      }
    }
  }
  for(int h=0;h<nirreps;h++)
  {
    free_block(t[h]);
  }
  free(t);
  return energy;
}

//Build iterative fock matrix
void newfockmatrix(double ***matrix,double ***corehamiltonian,double ***density,double *tei,int *sopi,int nirreps,int *symmo,bool outofcore)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project11

  if(outofcore==false)
  {
    int mu,nu,lambda,sigma,munu,lambdasigma,munulambdasigma,mulambda,nusigma,mulambdanusigma;
    int MU,NU,LAMBDA,SIGMA;
    for(int hmn=0;hmn<nirreps;hmn++)
    {
      for(mu=0;mu<sopi[hmn];mu++)
      {
        MU=mu+symmo[hmn];
        for(nu=0;nu<sopi[hmn];nu++)
        {
          matrix[hmn][mu][nu]=corehamiltonian[hmn][mu][nu];
          NU=nu+symmo[hmn];
          for(int hls=0;hls<nirreps;hls++)
          {
            for(lambda=0;lambda<sopi[hls];lambda++)
            {
              LAMBDA=lambda+symmo[hls];
              for(sigma=0;sigma<sopi[hls];sigma++)
              {
                SIGMA=sigma+symmo[hls];
                munu = INDEX(MU,NU);
                lambdasigma = INDEX(LAMBDA,SIGMA);
                munulambdasigma = INDEX(munu,lambdasigma);
                mulambda = INDEX(MU,LAMBDA);
                nusigma = INDEX(NU,SIGMA);
                mulambdanusigma = INDEX(mulambda,nusigma);
                matrix[hmn][mu][nu] += density[hls][lambda][sigma] * ( (2.0*tei[munulambdasigma]) - tei[mulambdanusigma] );
              }
            }
          }
        }
      }
    }
  }
  else
  {
    struct iwlbuf InBuf;
    Value *valptr;
    Label *lblptr;
    int lastbuf=0;
    int idx;
  
    for(int h=0;h<nirreps;h++)
    {
      for(int i=0;i<sopi[h];i++)
      {
        for(int j=0;j<sopi[h];j++)
        {
          matrix[h][i][j] = corehamiltonian[h][i][j];
        }
      }
    }

    double integral;
    //Obtain Two-Electron-Integrals
    iwl_buf_init(&InBuf,PSIF_SO_TEI,1e-14,1,0);
    do
    {
      iwl_buf_fetch(&InBuf);
      lblptr = InBuf.labels;
      valptr = InBuf.values;
      lastbuf = InBuf.lastbuf;
      for(idx=4*InBuf.idx;InBuf.idx<InBuf.inbuf;InBuf.idx++)
      {
        int i,j,k,l;
        i = abs( (int) lblptr[idx++]);
        j = (int) lblptr[idx++];
        k = (int) lblptr[idx++];
        l = (int) lblptr[idx++];
        integral = (double) valptr[InBuf.idx];    
      
        int iloc=0,jloc=0,kloc=0,lloc=0;
        int I,J,K,L;
        for(int h=0;h<nirreps;h++)
        {
          if(sopi[h]!=0 && i>=symmo[h])
          {
            iloc=h;
          }
          if(sopi[h]!=0 && j>=symmo[h])
          {
            jloc=h;
          }
          if(sopi[h]!=0 && k>=symmo[h])
          {
            kloc=h;
          }
          if(sopi[h]!=0 && l>=symmo[h])
          {
            lloc=h;
          }
        }
        I=i-symmo[iloc];
        J=j-symmo[jloc];
        K=k-symmo[kloc];
        L=l-symmo[lloc];     
 
        if(iloc == jloc && kloc==lloc)
        {
          if(i==j && j==k && k==l)
          {
            matrix[iloc][I][J] += 2*density[kloc][K][L]*integral;
          }
          else if(i==j && k==l)
          {
            matrix[iloc][I][J] += 2*density[kloc][K][L]*integral;          
            matrix[kloc][K][L] += 2*density[iloc][I][J]*integral;
          }
          else if(i==k && j==l)
          {
            matrix[iloc][I][J] += 2*density[kloc][K][L]*integral;
            matrix[iloc][I][J] += 2*density[lloc][L][K]*integral;
            matrix[jloc][J][I] += 2*density[kloc][K][L]*integral;
            matrix[jloc][J][I] += 2*density[lloc][L][K]*integral;
          } 
          else if(i==j)
          {
            matrix[iloc][I][J] += 2*density[kloc][K][L]*integral;
            matrix[iloc][I][J] += 2*density[lloc][L][K]*integral;
            matrix[kloc][K][L] += 2*density[iloc][I][J]*integral;
            matrix[lloc][L][K] += 2*density[iloc][I][J]*integral;
          }
          else if(k==l)
          {
            matrix[iloc][I][J] += 2*density[kloc][K][L]*integral;
            matrix[jloc][J][I] += 2*density[kloc][K][L]*integral;
            matrix[kloc][K][L] += 2*density[iloc][I][J]*integral;
            matrix[kloc][K][L] += 2*density[jloc][J][I]*integral;
          }
          else
          {
            matrix[iloc][I][J] += 2*density[kloc][K][L]*integral;
            matrix[iloc][I][J] += 2*density[lloc][L][K]*integral;
            matrix[jloc][J][I] += 2*density[kloc][K][L]*integral;
            matrix[jloc][J][I] += 2*density[lloc][L][K]*integral;
            matrix[kloc][K][L] += 2*density[iloc][I][J]*integral;
            matrix[kloc][K][L] += 2*density[jloc][J][I]*integral;
            matrix[lloc][L][K] += 2*density[iloc][I][J]*integral;
            matrix[lloc][L][K] += 2*density[jloc][J][I]*integral;
          }
        }
        if(iloc == kloc && jloc==lloc)
        {
          if(i==j && j==k && k==l)
          {
            matrix[iloc][I][K] -= density[jloc][J][L]*integral;
          }
          else if(i==j && k==l)
          {
            matrix[iloc][I][K] -= density[jloc][J][L]*integral;
            matrix[kloc][K][I] -= density[lloc][L][J]*integral;
          }
          else if(i==k && j==l)
          {
            matrix[iloc][I][K] -= density[jloc][J][L]*integral;
            matrix[jloc][J][L] -= density[iloc][I][K]*integral; 
          } 
          else if(i==j)
          {
            matrix[iloc][I][K] -= density[jloc][J][L]*integral;
            matrix[kloc][K][I] -= density[lloc][L][J]*integral;
          }
          else if(k==l)
          {
            matrix[iloc][I][K] -= density[jloc][J][L]*integral;
            matrix[kloc][K][I] -= density[lloc][L][J]*integral;
          }
          else
          {
            matrix[iloc][I][K] -= density[jloc][J][L]*integral;
            matrix[kloc][K][I] -= density[lloc][L][J]*integral;
            matrix[lloc][L][J] -= density[kloc][K][I]*integral;
            matrix[jloc][J][L] -= density[iloc][I][K]*integral;
          }
        }
        if(iloc == lloc && jloc==kloc)
        {
          if(i==j && j==k && k==l)
          {
          }
	  else if(i==j && k==l)
          {
          }
          else if(i==k && j==l)
          {
            matrix[jloc][J][K] -= density[iloc][I][L]*integral;
            matrix[iloc][I][L] -= density[jloc][J][K]*integral;
          } 
          else if(i==j)
          {
            matrix[iloc][I][L] -= density[jloc][J][K]*integral;
            matrix[lloc][L][I] -= density[kloc][K][J]*integral;
          }
          else if(k==l)
          {
            matrix[jloc][J][K] -= density[iloc][I][L]*integral;
            matrix[kloc][K][J] -= density[lloc][L][I]*integral;
          }
          else
          {
            matrix[jloc][J][K] -= density[iloc][I][L]*integral;
            matrix[lloc][L][I] -= density[jloc][K][J]*integral;
            matrix[iloc][I][L] -= density[kloc][J][K]*integral;
            matrix[kloc][K][J] -= density[lloc][L][I]*integral;
          }
        }        
      }
    }while( !lastbuf );

    iwl_buf_close(&InBuf,1);
  }
}

//Calculate root mean square of 3D array
double rms3d(double ***oldmatrix,double ***newmatrix, int *sopi,int nirreps)
{
  double data = 0.0;
  double ***t;
  t = ( double *** ) malloc(nirreps*sizeof(double**));
  for(int h=0;h<nirreps;h++)
  {
    t[h] = block_matrix(sopi[h],sopi[h]);
  }
  for(int h=0;h<nirreps;h++)
  {
    for(int i=0;i<sopi[h];i++)
    {  
      for(int j=0;j<sopi[h];j++)
      {
        t[h][i][j] = newmatrix[h][i][j] - oldmatrix[h][i][j];
        data += t[h][i][j] * t[h][i][j];
        oldmatrix[h][i][j] = newmatrix[h][i][j];
      }
    }
  }
  data = sqrt(data);
  for(int i=0;i<nirreps;i++)
  {
    free_block(t[i]);
  }
  free(t);
  return data;
}


//Convert two electron integral to 4D array
void TEIarrayconverter(double *oned,double ****fourd,int nso)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project5

  int pr,qs,prqs,ps,qr,psqr;
  double value1,value2;
  for(int p=0;p<nso;p++)
  {
    for(int q=0;q<nso;q++)
    {
      for(int r=0;r<nso;r++)
      {
        for(int s=0;s<nso;s++)
        {
          pr = INDEX(p/2,r/2);
          qs = INDEX(q/2,s/2);
          prqs = INDEX(pr,qs);
          value1 = oned[prqs] * (p%2 == r%2) * (q%2 == s%2);
          ps = INDEX(p/2,s/2);
          qr = INDEX(q/2,r/2);
          psqr = INDEX(ps,qr);
          value2 = oned[psqr] * (p%2 == s%2) * (q%2 == r%2);
          fourd[p][q][r][s] = value1 - value2;
        }
      }
    }
  }
}

//Convert TEI to Molecular Orbital Basis
void buildteimo(double *arrAO,double *arrMO,double ***eigenvec,int *sopi,int nirreps,int *symmo,int nao, bool smart)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project4

  //Smart Algorithm -- Scales as N^5
  if(smart==1)
  {
    double **X,**Y,**TMP;
    TMP=block_matrix((nao*(nao+1)/2),(nao*(nao+1)/2));
    int I,J,K,L;
    int ij,kl;
    int ijkl,klij;
    for(int i=0,ij=0;i<nao;i++)
    {
      for(int j=0;j<=i;j++,ij++)
      {
        for(int hk=0;hk<nirreps;hk++)
        {
          for(int hl=0;hl<nirreps;hl++)
          {
            X=block_matrix(sopi[hk],sopi[hl]);
            Y=block_matrix(sopi[hk],sopi[hl]);
            for(int k=0;k<sopi[hk];k++)
            {
              K=k+symmo[hk];
              for(int l=0;l<sopi[hl];l++)
              {
                L=l+symmo[hl];
                kl=INDEX(K,L);
                ijkl=INDEX(ij,kl);
                X[k][l]=arrAO[ijkl];
              }
            }
            if(sopi[hk]!=0 && sopi[hl]!=0)
            {
              zero_mat(Y,sopi[hk],sopi[hl]);
              mmult(eigenvec[hk],1,X,0,Y,0,sopi[hk],sopi[hk],sopi[hl],0);
              zero_mat(X,sopi[hk],sopi[hl]);
              mmult(Y,0,eigenvec[hl],0,X,0,sopi[hk],sopi[hl],sopi[hl],0);
            }
            for(int k=0;k<sopi[hk];k++)
            {
              K=k+symmo[hk];
              for(int l=0;l<sopi[hl];l++)
              {
                L=l+symmo[hl];
                kl=INDEX(K,L);
                TMP[kl][ij]=X[k][l];
              }
            }
            free_block(X);
            free_block(Y);
          }  
        }
      }
    } 
    for(int k=0,kl=0;k<nao;k++)
    {
      for(int l=0;l<=k;l++,kl++)
      {
        for(int hi=0;hi<nirreps;hi++)
        {
          for(int hj=0;hj<nirreps;hj++)
          {
            X=block_matrix(sopi[hi],sopi[hj]);
            Y=block_matrix(sopi[hi],sopi[hj]);
            for(int i=0;i<sopi[hi];i++)
            {
              I=i+symmo[hi];
              for(int j=0;j<sopi[hj];j++)
              {
                J=j+symmo[hj];
                ij=INDEX(I,J);
                X[i][j]=TMP[kl][ij];
              }
            }
            if(sopi[hi]!=0 && sopi[hj]!=0)
            {
              zero_mat(Y,sopi[hi],sopi[hj]);
              mmult(eigenvec[hi],1,X,0,Y,0,sopi[hi],sopi[hi],sopi[hj],0);
              zero_mat(X,sopi[hi],sopi[hj]);
              mmult(Y,0,eigenvec[hj],0,X,0,sopi[hi],sopi[hj],sopi[hj],0);
            }
            for(int i=0;i<sopi[hi];i++)
            {
              I=i+symmo[hi];
              for(int j=0;j<sopi[hj];j++)
              {
                J=j+symmo[hj];
                ij=INDEX(I,J);
                klij=INDEX(kl,ij);
                arrMO[klij]=X[i][j];
              }
            }
            free_block(X);
            free_block(Y);
          }
        }
      }
    }  
    free_block(TMP);
  }
  //Noddy Algorithm -- Scales as N^8
  else
  {
    int I,J,K,L;
    int ijkl,pqrs;
    int ij,kl;
    int pq,rs;
    for(int a=0;a<nirreps;a++)
    {
      for(int i=0;i<sopi[a];i++)
      {
        for(int p=0;p<sopi[a];p++)
        {
          I = i+symmo[a];
          for(int b=0;b<nirreps;b++)
          {
            for(int j=0;j< ( (I-symmo[b]+1)<sopi[b] ? (I-symmo[b]+1):sopi[b] ) ;j++)
            {
              J=j+symmo[b];
              for(int q=0;q<sopi[b];q++)
              {
                for(int c=0;c<nirreps;c++)
                {
                  for(int k=0;k<( ((I-symmo[c]+1)<sopi[c] ? (I-symmo[c]+1):sopi[c]) );k++)
                  {
                    K=k+symmo[c];
                    for(int r=0;r<sopi[c];r++)
                    {
                      for(int d=0;d<nirreps;d++)
                      {
                        for(int l=0,limit=( (I==K)?(J-symmo[d]+1):(K-symmo[d]+1) );l< ( (limit<sopi[d])?(limit):(sopi[d]) );l++)
                        {
                          L=l+symmo[d];
                          for(int s=0;s<sopi[d];s++)
                          {
                            ij = INDEX(I,J); kl=INDEX(K,L);
                            pq = INDEX(p+symmo[a],q+symmo[b]);
                            rs = INDEX(r+symmo[c],s+symmo[d]);
                            ijkl = INDEX(ij,kl);
                            pqrs = INDEX(pq,rs);
                            arrMO[ijkl] += eigenvec[a][p][i]*eigenvec[b][q][j]*eigenvec[c][r][k]*eigenvec[d][s][l]*arrAO[pqrs];
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }  
  }
}

//Calculate the Moller Plesset 2 Energy
double calculateMP2(double *arrMO,double **eigenvalues,int *sopi,int nirreps,int *filled,int *symmo)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project4

  int ia,ja,jb,ib,iajb,ibja;
  int I,J,A,B;
  double energy=0.0;
  for(int hi=0;hi<nirreps;hi++)
  {
    for(int i=0;i<filled[hi];i++)
    {
      I=i+symmo[hi];
      for(int ha=0;ha<nirreps;ha++)
      {
        for(int a=filled[ha];a<sopi[ha];a++)
        {
          A=a+symmo[ha];
          ia = INDEX(I,A);
          for(int hj=0;hj<nirreps;hj++)
          {
            for(int j=0;j<filled[hj];j++)
            {
              J=j+symmo[hj];
              ja = INDEX(J,A);
              for(int hb=0;hb<nirreps;hb++)
              {
                for(int b=filled[hb];b<sopi[hb];b++)
                {
                  B=b+symmo[hb];
                  jb=INDEX(J,B);
                  ib=INDEX(I,B);
                  iajb=INDEX(ia,jb);
                  ibja=INDEX(ib,ja);
                  energy += arrMO[iajb] * ( (2.0*arrMO[iajb])-arrMO[ibja] )/(eigenvalues[hi][i] - eigenvalues[ha][a] + eigenvalues[hj][j] - eigenvalues[hb][b] );
                }
              }
            }
          }
        }
      }
    }
  }
  return energy;
}

//Build tau and tautilde matrices
void buildtaus(double ****tt,double ****t,double **t1,double ****t2,int *sopi,int nirreps,int *symmo)
{
  //equations 9 and 10 from J.F. Stanton, J. Gauss, J.D. Watts, and R.J. Bartlett J.Chem.Phys. volume 94, pg 4334-4345 (1991)
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project5

  //tt=tautilde, t=tau
  int I,J,K,L;
  for(int hi=0;hi<nirreps;hi++)
  {
    for(int i=0;i<2*sopi[hi];i++)
    {
      I= i+2*symmo[hi];
      for(int hj=0;hj<nirreps;hj++)
      {
        for(int j=0;j<2*sopi[hj];j++)
        {
          J=j+2*symmo[hj];
          for(int hk=0;hk<nirreps;hk++)
          {
            for(int k=0;k<2*sopi[hk];k++)
            {
              K=k+2*symmo[hk];
              for(int hl=0;hl<nirreps;hl++)
              {
                for(int l=0;l<2*sopi[hl];l++)
                {
                  L=l+2*symmo[hl];
                  tt[I][J][K][L] = t2[I][J][K][L]+(0.5)*(t1[I][K]*t1[J][L]-t1[J][K]*t1[I][L]);
                  t[I][J][K][L] = t2[I][J][K][L] + (t1[I][K]*t1[J][L])-(t1[J][K]*t1[I][L]);
                }
              }
            }
          }
        }
      }
    } 
  } 
}

//Build F or newF matrix
void buildF(double **matrix,double ***sfock,double **t1,double ****so,double ****tt,int *sopi,int *filled,int nirreps,int nso,int *symmo)
{
  //equations 3,4,5 from J.F. Stanton, J. Gauss, J.D. Watts, and R.J. Bartlett J.Chem.Phys. volume 94, pg 4334-4345 (1991)
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project5
  for(int i=0;i<nso;i++)
  {
    for(int j=0;j<nso;j++)
    {
      matrix[i][j] = 0.0;
    }
  }

  //equation 3
  int A,E,M,F,N,I;
  for(int ha=0;ha<nirreps;ha++)
  {
    for(int a=2*filled[ha];a<2*sopi[ha];a++)
    {
      A=a+2*symmo[ha];
      for(int he=0;he<nirreps;he++)
      {
        for(int e=2*filled[he];e<2*sopi[he];e++)
        {
          E=e+2*symmo[he];
          if(a!=e && ha==he)
          {      
            matrix[A][E] += sfock[ha][a][e]; 
          }
          for(int hm=0;hm<nirreps;hm++)
          {
            for(int m=0;m<2*filled[hm];m++)
            {
              M=m+2*symmo[hm];
              if(hm==he)
              {
                matrix[A][E]+=(-0.5*sfock[hm][m][e]*t1[A][M]);
              }
              for(int hf=0;hf<nirreps;hf++)
              {
                for(int f=2*filled[hf];f<2*sopi[hf];f++)
                {
                  F=f+2*symmo[hf];
                  matrix[A][E]+=t1[F][M]*so[M][A][F][E];
                  for(int hn=0;hn<nirreps;hn++)
                  {
                    for(int n=0;n<2*filled[hn];n++)
                    {
                      N=n+2*symmo[hn];
                      matrix[A][E] += (-0.5)*tt[A][F][M][N]*so[M][N][E][F];
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  //equation 4
  for(int hm=0;hm<nirreps;hm++)
  {
    for(int m=0;m<2*filled[hm];m++)
    {
      M=m+2*symmo[hm];
      for(int hi=0;hi<nirreps;hi++)
      {
        for(int i=0;i<2*filled[hi];i++)
        {
          I=i+2*symmo[hi];
          if(m!=i && hm==hi)
          {  
            matrix[M][I] += sfock[hm][m][i];
          }
          for(int he=0;he<nirreps;he++)
          {
            for(int e=2*filled[he];e<2*sopi[he];e++)
            {
              E=e+2*symmo[he];
              if(hm==he)
              {           
                matrix[M][I] += (0.5)*(t1[E][I]*sfock[hm][m][e]);
              }
              for(int hn=0;hn<nirreps;hn++)
              {
                for(int n=0;n<2*filled[hn];n++)
                {
                  N=n+2*symmo[hn];
                  matrix[M][I] += t1[E][N]*so[M][N][I][E];
                  for(int hf=0;hf<nirreps;hf++)
                  {
                    for(int f=2*filled[hf];f<2*sopi[hf];f++)
                    {
                      F=f+2*symmo[hf];
                      matrix[M][I] += (0.5)*tt[E][F][I][N]*so[M][N][E][F];
                    }
                  } 
                }
              }
            }
          }
        }
      }
    }
  }
  //equation 5
  for(int hm=0;hm<nirreps;hm++)
  {
    for(int m=0;m<2*filled[hm];m++)
    {
      M=m+2*symmo[hm];
      for(int he=0;he<nirreps;he++)
      {
        for(int e=2*filled[he];e<2*sopi[he];e++)
        {
          E=e+2*symmo[he];
          if(hm==he)
          {
            matrix[M][E]+=sfock[hm][m][e];
          }
          for(int hn=0;hn<nirreps;hn++)
          {
            for(int n=0;n<2*filled[hn];n++)
            {
              N=n+2*symmo[hn];
              for(int hf=0;hf<nirreps;hf++)
              {
                for(int f=2*filled[hf];f<2*sopi[hf];f++)
                {
                  F=f+2*symmo[hf];
                  matrix[M][E] += t1[F][N]*so[M][N][E][F];
                }
              }
            }
          }
        }
      }
    }
  }
}

//Build W or newW matrix
void buildW(double ****matrix,double ****so,double **t1,double ****t2,double ****tau,int nso,int *filled,int *sopi,int nirreps,int *symmo)
{
  //equations 6,7,8 from J.F. Stanton, J. Gauss, J.D. Watts, and R.J. Bartlett J.Chem.Phys. volume 94, pg 4334-4345 (1991)
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project5
  for(int i=0;i<nso;i++)
  {
    for(int j=0;j<nso;j++)
    {
      for(int k=0;k<nso;k++)
      {
        for(int l=0;l<nso;l++)
        {
          matrix[i][j][k][l]=0.0;
        }
      }
    }
  }

  //equation 6
  int M,N,I,J,E,F,A,B;
  for(int hm=0;hm<nirreps;hm++)
  {
    for(int m=0;m<2*filled[hm];m++)
    {
      M=m+2*symmo[hm];
      for(int hn=0;hn<nirreps;hn++)
      {
        for(int n=0;n<2*filled[hn];n++)
        {
          N=n+2*symmo[hn];
          for(int hi=0;hi<nirreps;hi++)
          {
            for(int i=0;i<2*filled[hi];i++)
            {
              I=i+2*symmo[hi];
              for(int hj=0;hj<nirreps;hj++)
              {
                for(int j=0;j<2*filled[hj];j++)
                {
                  J=j+2*symmo[hj];
                  matrix[M][N][I][J] += so[M][N][I][J];
                  for(int he=0;he<nirreps;he++)
                  {
                    for(int e=2*filled[he];e<2*sopi[he];e++)
                    {
                      E=e+2*symmo[he];
                      matrix[M][N][I][J] += t1[E][J]*so[M][N][I][E];
                      matrix[M][N][I][J] += (-1.0)*t1[E][I]*so[M][N][J][E];
                      for(int hf=0;hf<nirreps;hf++)
                      {
                        for(int f=2*filled[hf];f<2*sopi[hf];f++)
                        {
                          F=f+2*symmo[hf];
                          matrix[M][N][I][J] += (0.25)*tau[E][F][I][J]*so[M][N][E][F]; 
                        }
                      }
                    }
                  }
                }
              }  
            }
          }
        }
      }
    }
  }
  //equation 7
  for(int ha=0;ha<nirreps;ha++)
  {
    for(int a=2*filled[ha];a<2*sopi[ha];a++)
    {
      A=a+2*symmo[ha];
      for(int hb=0;hb<nirreps;hb++)
      {
        for(int b=2*filled[hb];b<2*sopi[hb];b++)
        {
          B=b+2*symmo[hb];
          for(int he=0;he<nirreps;he++)
          {
            for(int e=2*filled[he];e<2*sopi[he];e++)
            {
              E=e+2*symmo[he];
              for(int hf=0;hf<nirreps;hf++)
              {
                for(int f=2*filled[hf];f<2*sopi[hf];f++)
                {
                  F=f+2*symmo[hf];
                  matrix[A][B][E][F] += so[A][B][E][F];
                  for(int hm=0;hm<nirreps;hm++)
                  {
                    for(int m=0;m<2*filled[hm];m++)
                    {
                      M=m+2*symmo[hm];
                      matrix[A][B][E][F] -= ( t1[B][M]*so[A][M][E][F]-t1[A][M]*so[B][M][E][F] );
                      for(int hn=0;hn<nirreps;hn++)
                      {
                        for(int n=0;n<2*filled[hn];n++)
                        {
                          N=n+2*symmo[hn];
                          matrix[A][B][E][F] += (0.25)*tau[A][B][M][N]*so[M][N][E][F];
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  //equation 8
  for(int hm=0;hm<nirreps;hm++)
  {
    for(int m=0;m<2*filled[hm];m++)
    {
      M=m+2*symmo[hm];
      for(int hb=0;hb<nirreps;hb++)
      {
        for(int b=2*filled[hb];b<2*sopi[hb];b++)
        {
          B=b+2*symmo[hb];
          for(int he=0;he<nirreps;he++)
          {  
            for(int e=2*filled[he];e<2*sopi[he];e++)
            {
              E=e+2*symmo[he];
              for(int hj=0;hj<nirreps;hj++)
              {
                for(int j=0;j<2*filled[hj];j++)
                {
                  J=j+2*symmo[hj];
                  matrix[M][B][E][J] += so[M][B][E][J];
                  for(int hf=0;hf<nirreps;hf++)
                  {
                    for(int f=2*filled[hf];f<2*sopi[hf];f++)
                    {
                      F=f+2*symmo[hf];
                      matrix[M][B][E][J] += t1[F][J]*so[M][B][E][F];
                      for(int hn=0;hn<nirreps;hn++)
                      {
                        for(int n=0;n<2*filled[hn];n++)
                        {
                          N=n+2*symmo[hn];
                          matrix[M][B][E][J] += (-1.0)*( (0.5*t2[F][B][J][N])+(t1[F][J]*t1[B][N]) )*so[M][N][E][F];
                        }
                      }
                    }
                  }
                  for(int hn=0;hn<nirreps;hn++)
                  {
                    for(int n=0;n<2*filled[hn];n++)
                    {
                      N=n+2*symmo[hn];
                      matrix[M][B][E][J] += (-1.0)*(t1[B][N]*so[M][N][E][J]);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }  
}

//Build T1 Cluster amplitudes
void buildT1(double **matrix,double ***sfock,double **t1,double ****t2,double ****so,double **Fmat,int nso,int *filled,int *sopi,int nirreps,int *symmo)
{
  //equations 1 from J.F. Stanton, J. Gauss, J.D. Watts, and R.J. Bartlett J.Chem.Phys. volume 94, pg 4334-4345 (1991)
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project5
  for(int i=0;i<nso;i++)
  {
    for(int j=0;j<nso;j++)
    {
      matrix[i][j] = 0.0;
    }
  }

  int A,I,M,E,F,N;
  for(int z=0;z<nirreps;z++)
  {
    for(int a=2*filled[z];a<2*sopi[z];a++)
    {
      A=a+2*symmo[z];
      for(int y=0;y<nirreps;y++)
      {
        for(int i=0;i<2*filled[y];i++)
        {
          I=i+2*symmo[y];
          if(z==y)
          {
            matrix[A][I] += sfock[y][i][a];
          }
          for(int x=0;x<nirreps;x++)
          {
            for(int e=2*filled[x];e<2*sopi[x];e++)
            {
              E=e+2*symmo[x];
              matrix[A][I] += t1[E][I]*Fmat[A][E];
            }
          }
          for(int w=0;w<nirreps;w++)
          {
            for(int m=0;m<2*filled[w];m++)
            {
              M=m+2*symmo[w];
              matrix[A][I] -= t1[A][M]*Fmat[M][I];
              for(int u=0;u<nirreps;u++)
              {
                for(int e=2*filled[u];e<2*sopi[u];e++)
                {
                  E=e+2*symmo[u];
                  matrix[A][I] += t2[A][E][I][M]*Fmat[M][E];
                  for(int t=0;t<nirreps;t++)
                  {
                    for(int f=2*filled[t];f<2*sopi[t];f++)
                    {
                      F=f+2*symmo[t];
                      matrix[A][I] -= (0.5)*t2[E][F][I][M]*so[M][A][E][F];
                    }
                  }
                  for(int r=0;r<nirreps;r++)
                  {
                    for(int n=0;n<2*filled[r];n++)
                    {
                      N=n+2*symmo[r];
                      matrix[A][I] -= (0.5)*t2[A][E][M][N]*so[N][M][E][I];
                    }
                  }
                }
              } 
            }
          }
          for(int q=0;q<nirreps;q++)
          {
            for(int n=0;n<2*filled[q];n++)
            {
              N=n+2*symmo[q];
              for(int p=0;p<nirreps;p++)
              {
                for(int f=2*filled[p];f<2*sopi[p];f++)
                {
                  F=f+2*symmo[p];
                  matrix[A][I] -= t1[F][N]*so[N][A][I][F];
                }
              }
            }
          }      
          matrix[A][I] /= (sfock[y][i][i] - sfock[z][a][a]);
        }
      }
    }
  }
}
  
//Build T2 Cluster Amplitude
void buildT2(double ****matrix,double ***sfock,double ****so,double **t1,double ****t2,double ****t,double **Fmat,double ****Wmat,int nso,int *filled,int *sopi,int nirreps,int *symmo)
{
  //equation 2 from J.F. Stanton, J. Gauss, J.D. Watts, and R.J. Bartlett J.Chem.Phys. volume 94, pg 4334-4345 (1991)
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project5
  for(int i=0;i<nso;i++)
  {
    for(int j=0;j<nso;j++)
    {
      for(int k=0;k<nso;k++)
      {
        for(int l=0;l<nso;l++)
        {
          matrix[i][j][k][l] = 0.0;
        }
      }
    }
  }

  int A,B,I,J,E,M,F,N;
  for(int ha=0;ha<nirreps;ha++)
  {
    for(int a=2*filled[ha];a<2*sopi[ha];a++)
    {
      A=a+2*symmo[ha];
      for(int hb=0;hb<nirreps;hb++)
      {
        for(int b=2*filled[hb];b<2*sopi[hb];b++)
        {
          B=b+2*symmo[hb];
          for(int hi=0;hi<nirreps;hi++)
          {
            for(int i=0;i<2*filled[hi];i++)
            {
              I=i+2*symmo[hi];
              for(int hj=0;hj<nirreps;hj++)
              {
                for(int j=0;j<2*filled[hj];j++)
                {
                  J=j+2*symmo[hj];
                  matrix[A][B][I][J] += so[I][J][A][B];
                  for(int he=0;he<nirreps;he++)
                  {
                    for(int e=2*filled[he];e<2*sopi[he];e++)
                    {
                      E=e+2*symmo[he];
                      matrix[A][B][I][J] += t2[A][E][I][J]*Fmat[B][E];
                      matrix[A][B][I][J] -= t2[B][E][I][J]*Fmat[A][E];
                      for(int hm=0;hm<nirreps;hm++)
                      {
                        for(int m=0;m<2*filled[hm];m++)
                        {
                          M=m+2*symmo[hm];
                          matrix[A][B][I][J] -= (0.5)*(t2[A][E][I][J]*t1[B][M]*Fmat[M][E]);
                          matrix[A][B][I][J] += (0.5)*(t2[B][E][I][J]*t1[A][M]*Fmat[M][E]); 
                        }
                      }
                      for(int hf=0;hf<nirreps;hf++)
                      {
                        for(int f=2*filled[hf];f<2*sopi[hf];f++)  
                        {
                          F=f+2*symmo[hf];
                          matrix[A][B][I][J] += (0.5)*t[E][F][I][J]*Wmat[A][B][E][F];
                        }
                      }
                      matrix[A][B][I][J] += t1[E][I]*so[A][B][E][J];
                      matrix[A][B][I][J] -= t1[E][J]*so[A][B][E][I];
                    }
                  }
                  for(int hm=0;hm<nirreps;hm++)
                  {
                    for(int m=0;m<2*filled[hm];m++)
                    {
                      M=m+2*symmo[hm];
                      matrix[A][B][I][J] -= t2[A][B][I][M]*Fmat[M][J];
                      matrix[A][B][I][J] += t2[A][B][J][M]*Fmat[M][I];  
                      for(int he=0;he<nirreps;he++)
                      {
                        for(int e=2*filled[he];e<2*sopi[he];e++)
                        {
                          E=e+2*symmo[he];
                          matrix[A][B][I][J] -= (0.5)*(t2[A][B][I][M]*t1[E][J]*Fmat[M][E]);
                          matrix[A][B][I][J] += (0.5)*(t2[A][B][J][M]*t1[E][I]*Fmat[M][E]);
                          matrix[A][B][I][J] += (t2[A][E][I][M]*Wmat[M][B][E][J])-(t1[E][I]*t1[A][M]*so[M][B][E][J]);
                          matrix[A][B][I][J] -= (t2[A][E][J][M]*Wmat[M][B][E][I])-(t1[E][J]*t1[A][M]*so[M][B][E][I]);
                          matrix[A][B][I][J] -= (t2[B][E][I][M]*Wmat[M][A][E][J])-(t1[E][I]*t1[B][M]*so[M][A][E][J]);
                          matrix[A][B][I][J] += (t2[B][E][J][M]*Wmat[M][A][E][I])-(t1[E][J]*t1[B][M]*so[M][A][E][I]);
                        }
                      }
                      for(int hn=0;hn<nirreps;hn++)
                      {
                        for(int n=0;n<2*filled[hn];n++)
                        {
                          N=n+2*symmo[hn];
                          matrix[A][B][I][J] += (0.5)*t[A][B][M][N]*Wmat[M][N][I][J];
                        }
                      }
                      matrix[A][B][I][J] -= t1[A][M]*so[M][B][I][J];
                      matrix[A][B][I][J] += t1[B][M]*so[M][A][I][J];
                    }
                  }
                  matrix[A][B][I][J] /= (sfock[hi][i][i] + sfock[hj][j][j] - sfock[ha][a][a] - sfock[hb][b][b]);
                }
              }
            }
          }
        }
      }
    }
  }
}

//Calculate Coupled Clusters Energy
double calcEcc(double ***sfock,double ****so,double **t1,double ****t2,int *filled,int *sopi,int nirreps,int *symmo)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project5

  double data=0.0;
  int I,A,J,B;
  for(int hi=0;hi<nirreps;hi++)
  {
    for(int i=0;i<2*filled[hi];i++)
    {
      I=i+2*symmo[hi];
      for(int ha=0;ha<nirreps;ha++)
      {
        for(int a=2*filled[ha];a<2*sopi[ha];a++)
        {
          A=a+2*symmo[ha];
          if(hi==ha)
          {
            data+=sfock[hi][i][a]*t1[A][I];
          }
          for(int hj=0;hj<nirreps;hj++)
          {
            for(int j=0;j<2*filled[hj];j++)
            {
              J=j+2*symmo[hj];
              for(int hb=0;hb<nirreps;hb++)
              {
                for(int b=2*filled[hb];b<2*sopi[hb];b++)
                {
                  B=b+2*symmo[hb];
                  data+=(0.25)*(so[I][J][A][B]*t2[A][B][I][J]);
                  data+=(0.5)*(so[I][J][A][B]*t1[A][I]*t1[B][J]);
                }
              }
            }
          }
        }
      }
    }
  }  
  return data;
}

//Calculate Triples Correction to Coupled Clusters Energy
double calcEt(double ***sfock,double **t1,double ****t2,double ****so,int *filled,int *sopi,int nirreps,int *symmo)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project6
  //Utilizing full storage of triples

  double Et=0.0;
  int I,J,K,A,B,C,E,M;
  for(int hi=0;hi<nirreps;hi++)
  {
    for(int i=0;i<2*filled[hi];i++)
    {
      I=i+2*symmo[hi];
      for(int hj=0;hj<nirreps;hj++)
      {
        for(int j=0;j<2*filled[hj];j++)
        {
          J=j+2*symmo[hj];
          for(int hk=0;hk<nirreps;hk++)
          {
            for(int k=0;k<2*filled[hk];k++)
            {
              K=k+2*symmo[hk];
              for(int ha=0;ha<nirreps;ha++)
              {
                for(int a=2*filled[ha];a<2*sopi[ha];a++)
                {
                  A=a+2*symmo[ha];
                  for(int hb=0;hb<nirreps;hb++)
                  {
                    for(int b=2*filled[hb];b<2*sopi[hb];b++)
                    {
                      B=b+2*symmo[hb];
                      for(int hc=0;hc<nirreps;hc++)
                      {
                        for(int c=2*filled[hc];c<2*sopi[hc];c++)
                        {
                          C=c+2*symmo[hc];
                          double conT=0.0;
                          double disconT=0.0;
                          double D = sfock[hi][i][i]+sfock[hj][j][j]+sfock[hk][k][k]-sfock[ha][a][a]-sfock[hb][b][b]-sfock[hc][c][c];

                          //Calculate disconnected triples
                          disconT += ( t1[A][I]*so[J][K][B][C] - t1[B][I]*so[J][K][A][C] - t1[C][I]*so[J][K][B][A] );
                          disconT += (-t1[A][J]*so[I][K][B][C] + t1[B][J]*so[I][K][A][C] + t1[C][J]*so[I][K][B][A] );
                          disconT += (-t1[A][K]*so[J][I][B][C] + t1[B][K]*so[J][I][A][C] + t1[C][K]*so[J][I][B][A] );
                          disconT /= D;

                          //Calculate connected triples
                          for(int he=0;he<nirreps;he++)
                          {
                            for(int e=2*filled[he];e<2*sopi[he];e++)
                            {
                              E=e+2*symmo[he];
                              conT += (  t2[A][E][J][K]*so[E][I][B][C] - t2[B][E][J][K]*so[E][I][A][C] - t2[C][E][J][K]*so[E][I][B][A] );
                              conT += ( -t2[A][E][I][K]*so[E][J][B][C] + t2[B][E][I][K]*so[E][J][A][C] + t2[C][E][I][K]*so[E][J][B][A] );
                              conT += ( -t2[A][E][J][I]*so[E][K][B][C] + t2[B][E][J][I]*so[E][K][A][C] + t2[C][E][J][I]*so[E][K][B][A] );
                            }
                          }
                          for(int hm=0;hm<nirreps;hm++)
                          {
                            for(int m=0;m<2*filled[hm];m++)
                            {
                              M=m+2*symmo[hm];
                              conT -= (  t2[B][C][I][M]*so[M][A][J][K] - t2[A][C][I][M]*so[M][B][J][K] - t2[B][A][I][M]*so[M][C][J][K] );
                              conT -= ( -t2[B][C][J][M]*so[M][A][I][K] + t2[A][C][J][M]*so[M][B][I][K] + t2[B][A][J][M]*so[M][C][I][K] );
                              conT -= ( -t2[B][C][K][M]*so[M][A][J][I] + t2[A][C][K][M]*so[M][B][J][I] + t2[B][A][K][M]*so[M][C][J][I] );
                            }
                          }
                          conT /= D;
         
                          Et += conT*D*(conT+disconT); 
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  Et /= 36;
  return Et;
}
                                            
//Calculate Root mean square of 4D matrix
double rms4d(double ****oldmatrix, double ****newmatrix,int nso,bool flag)
{
  double data=0.0;
  for(int i=0;i<nso;i++)
  {
    for(int j=0;j<nso;j++)
    {
      for(int k=0;k<nso;k++)
      {
        for(int l=0;l<nso;l++)
        {
          data+=(newmatrix[i][j][k][l]-oldmatrix[i][j][k][l])*(newmatrix[i][j][k][l]-oldmatrix[i][j][k][l]);
          //Save new matrix?
          if(flag==1)
          {
            oldmatrix[i][j][k][l]=newmatrix[i][j][k][l];
          }
        }
      }
    }
  }
  data = sqrt(data);
  return data;
}

//Calculate Error Matrix for SCF DIIS Extrapolation
void calcErrorSCF(double ***err,double ***fock,double ***density,double ***S,int *sopi,int nirreps)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project8

  double***t1,***t2; //temporary variables
  t1=(double ***) malloc(nirreps*sizeof(double**));
  t2=(double ***) malloc(nirreps*sizeof(double**));
  for(int h=0;h<nirreps;h++)
  {
    t1[h] = block_matrix(sopi[h],sopi[h]);
    t2[h] = block_matrix(sopi[h],sopi[h]);
    mmult(fock[h],0,density[h],0,t1[h],0,sopi[h],sopi[h],sopi[h],0);
    mmult(t1[h],0,S[h],0,t2[h],0,sopi[h],sopi[h],sopi[h],0);
  }

  double***t3,***t4; //temporary variables
  t3=(double***) malloc(nirreps*sizeof(double**));
  t4=(double***) malloc(nirreps*sizeof(double**));
  for(int h=0;h<nirreps;h++)
  {
    t3[h] = block_matrix(sopi[h],sopi[h]);
    t4[h] = block_matrix(sopi[h],sopi[h]);
    mmult(S[h],0,density[h],0,t3[h],0,sopi[h],sopi[h],sopi[h],0);
    mmult(t3[h],0,fock[h],0,t4[h],0,sopi[h],sopi[h],sopi[h],0);
  }
  for(int h=0;h<nirreps;h++)
  {
    for(int i=0;i<sopi[h];i++)
    {
      for(int j=0;j<sopi[h];j++)
      {
        err[h][i][j] = t2[h][i][j]-t4[h][i][j];
      }
    }
  }
  for(int h=0;h<nirreps;h++)
  {
    free_block(t1[h]);
    free_block(t2[h]);
    free_block(t3[h]);
    free_block(t4[h]);
  }
  free(t1);
  free(t2);
  free(t3);
  free(t4);
}
//Load 2D matrix into 3D matrix for storage
void load2d(double ***matrix3d,double **matrix,int maxheight,int height,int nso)
{
  if(height<=maxheight)
  {
    for(int i=0;i<nso;i++)
    {
      for(int j=0;j<nso;j++)
      {
        matrix3d[height][i][j] = matrix[i][j];
      }
    }
  }
}
//Load 3D Matrix into 4D matrix for storage
void load3dsymm(double ****matrix4d,double ***matrix,int maxheight,int height,int *sopi,int nirreps)
{
  if(height<=maxheight)
  {
    for(int h=0;h<nirreps;h++)
    {
      for(int i=0;i<sopi[h];i++)
      {
        for(int j=0;j<sopi[h];j++)
        {
          matrix4d[h][height][i][j]=matrix[h][i][j];
        }
      }
    }
  }
}
//Load 4D matrix into 5D matrix for storage
void load4d(double *****matrix5d,double ****matrix,int maxheight,int height,int nso)
{
  if(height<=maxheight)
  {
    for(int i=0;i<nso;i++)
    {
      for(int j=0;j<nso;j++)
      {
        for(int k=0;k<nso;k++)
        {
          for(int l=0;l<nso;l++)
          {
            matrix5d[height][i][j][k][l] = matrix[i][j][k][l];
          }
        }
      }
    }
  }
}

//Build B matrix
void buildBSCF(double **matrix,double ****errormatrix,int *sopi,int nirreps,int bsize)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project8

  for(int i=0;i<bsize-1;i++)
  {
    for(int j=0;j<bsize-1;j++)
    {
      matrix[i][j]=0.0;
      for(int h=0;h<nirreps;h++)
      {
        for(int a=0;a<sopi[h];a++)
        { 
          for(int b=0;b<sopi[h];b++)
          {
             matrix[i][j] += errormatrix[h][i][a][b]*errormatrix[h][j][a][b];
          }
        }
      }
    }
  }

  for(int i=0;i<bsize-1;i++)
  {
    matrix[bsize-1][i]=-1.0;
    matrix[i][bsize-1]=-1.0;
  }
  matrix[bsize-1][bsize-1]=0.0;
}

//Calculate root mean square of 2D array
double rms2d(double **oldmatrix,double **newmatrix, int nso,bool flag)
{
  double data = 0.0;
  double **t;
  t = ( double ** ) init_matrix(nso,nso);

  for(int i=0;i<nso;i++)
  {
    for(int j=0;j<nso;j++)
    {
      t[i][j] = newmatrix[i][j] - oldmatrix[i][j];
      data += t[i][j] * t[i][j];
      //Save new matrix?
      if(flag==1)
      {
        oldmatrix[i][j] = newmatrix[i][j];
      }
    }
  }
  data = sqrt(data);
  free_matrix(t,nso);
  return data;
}

//Calculate error for coupled cluster amplitudes
void calcErrorCC(double **errorT1,double ****errorT2,double **t1,double **newt1, double****t2,double ****newt2,int nso)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project10

  for(int i=0;i<nso;i++)
  {
    for(int j=0;j<nso;j++)
    {
      errorT1[i][j] = newt1[i][j] - t1[i][j];
      for(int k=0;k<nso;k++)
      {
        for(int l=0;l<nso;l++)
        {
          errorT2[i][j][k][l] = newt2[i][j][k][l] - t2[i][j][k][l];
        }
      }
    }
  }
}
void buildBCC(double **matrix,double ***itT1,double*****itT2,int bsize,int nso)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project10

  for(int i=0;i<bsize-1;i++)
  {
    for(int j=0;j<bsize-1;j++)
    {
       matrix[i][j] = 0.0;
       for(int a=0;a<nso;a++)
       {
         for(int b=0;b<nso;b++)
         {
           matrix[i][j] += itT1[i][a][b]*itT1[j][a][b];
         }
       }
    }
  }
  for(int i=0;i<bsize-1;i++)
  {
    for(int j=0;j<bsize-1;j++)
    {
      for(int a=0;a<nso;a++)
      {
        for(int b=0;b<nso;b++)
        {
          for(int c=0;c<nso;c++)
          {
            for(int d=0;d<nso;d++)
            {
              matrix[i][j] += itT2[i][a][b][c][d] * itT2[j][a][b][c][d];
            }
          }
        }
      } 
    }
  }
  for(int i=0;i<bsize-1;i++)
  {
    matrix[bsize-1][i] = -1.0;
    matrix[i][bsize-1] = -1.0;
  }
  matrix[bsize-1][bsize-1]=0.0;
}

void modifyT(double **matrix1,double ****matrix2,double ***it1,double *****it2,int nso, int bsize, double *right)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project10

  for(int i=0;i<nso;i++)
  {
    for(int j=0;j<nso;j++)
    {
      matrix1[i][j] = 0.0;
      for(int k=0;k<nso;k++)
      {
        for(int l=0;l<nso;l++)
        {
          matrix2[i][j][k][l] = 0.0;
        }
      }
    }
  }
  for(int i=0;i<nso;i++)
  {
    for(int j=0;j<nso;j++)
    {  
      for(int b=0;b<bsize-1;b++)
      {
        matrix1[i][j] += it1[b][i][j] * right[b];
      }
    }
  }
  for(int i=0;i<nso;i++)
  {
    for(int j=0;j<nso;j++)
    {
      for(int k=0;k<nso;k++)
      {
        for(int l=0;l<nso;l++)
        {
          for(int b=0;b<bsize-1;b++)
          {
            matrix2[i][j][k][l] += it2[b][i][j][k][l] * right[b]; 
          }
        }
      }
    }
  }
}

void bubblesort(double *array,int n)
{
  int flag=1;
  double temp;
  for(int i=0;i<n && flag;i++)
  {
    flag=0;
    for(int j=0;j<(n-1);j++)
    {
      if(array[j+1]<array[j])
      {
        temp = array[j];
        array[j]=array[j+1];
        array[j+1]=temp;
        flag=1;
      }
    }
  }
}

void checkmp2(double ****t2,double ****so,int nirreps,int *sopi,int *filled,int *symmo,double Emp2)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project5

  double data=0.0;
  int I,J,A,B;
  for(int hi=0;hi<nirreps;hi++)
  {
    for(int i=0;i<2*filled[hi];i++)
    {
      I=i+2*symmo[hi];
      for(int hj=0;hj<nirreps;hj++)
      {
        for(int j=0;j<2*filled[hj];j++)
        {
          J=j+2*symmo[hj];
          for(int ha=0;ha<nirreps;ha++)
          {
            for(int a=2*filled[ha];a<2*sopi[ha];a++)
            {
              A=a+2*symmo[ha];
              for(int hb=0;hb<nirreps;hb++)
              {
                for(int b=2*filled[hb];b<2*sopi[hb];b++)
                {
                  B=b+2*symmo[hb];
                  data+=t2[A][B][I][J]*so[I][J][A][B];
                }
              }
            }
          }
        }
      }
    }
  }
  data *= 0.25;
  if(fabs(Emp2-data)>1e-10)
  {
    std::cout << "Energies of MP2 Methods do not match!" << std::endl;
  }
}

void buildinitialT(double **t1,double ****t2,double ****so,double **evals,int nso,int nirreps,int *sopi,int *filled,int *symmo)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project5

  int I,J,K,L;
  for(int i=0;i<nso;i++)
  {
    for(int j=0;j<nso;j++)
    {
      t1[i][j]=0.0;
      for(int k=0;k<nso;k++)
      {
        for(int l=0;l<nso;l++)
        {
          t2[i][j][k][l] = 0.0;
        }
      }
    }
  }
  for(int hi=0;hi<nirreps;hi++)
  {
    for(int i=0;i<2*filled[hi];i++)
    {
      I=i+2*symmo[hi];
      for(int hj=0;hj<nirreps;hj++)
      {
        for(int j=0;j<2*filled[hj];j++)
        {
          J=j+2*symmo[hj];
          for(int hk=0;hk<nirreps;hk++)
          {
            for(int k=2*filled[hk];k<2*sopi[hk];k++)
            {
              K=k+2*symmo[hk];
              for(int hl=0;hl<nirreps;hl++)
              {
                for(int l=2*filled[hl];l<2*sopi[hl];l++)
                {
                  L=l+2*symmo[hl];
                  t2[K][L][I][J] = so[I][J][K][L] / (evals[hi][(i/2)] + evals[hj][(j/2)] - evals[hk][(k/2)] - evals[hl][(l/2)]);
                }
              }
            }
          }
        }
      }
    }
  }    
}

void buildsofock(double ***sfock,double **evals,int nirreps,int *sopi)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project5

  for(int i=0;i<nirreps;i++)
  {
    for(int j=0;j<2*sopi[i];j++)
    {
      for(int k=0;k<2*sopi[i];k++)
      {
        sfock[i][j][k]=0.0;
        if(j==k)
        {
          sfock[i][j][k] = evals[i][j/2];
        }
      }
    }  
  }    
}

void buildHCIS(double **matrix,double ***sfock,double ****so,int nirreps,int *sopi,int *filled, int *symmo)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project12

  int I,J,A,B;
  int index1=0,index2=0;
  for(int hi=0;hi<nirreps;hi++)
  {
   for(int i=0;i<2*filled[hi];i++)
    {
      I=i+2*symmo[hi];
      for(int ha=0;ha<nirreps;ha++)
      {
        for(int a=2*filled[ha];a<2*sopi[ha];a++)
        {
          A=a+2*symmo[ha];
          index2=0;
          for(int hj=0;hj<nirreps;hj++)
          {
            for(int j=0;j<2*filled[hj];j++)
            {
              J=j+2*symmo[hj];
              for(int hb=0;hb<nirreps;hb++)
              {
                for(int b=2*filled[hb];b<2*sopi[hb];b++)
                {
                  B=b+2*symmo[hb];          
                  if(I==J)
                  {
                    if(ha==hj)
                    {
                      matrix[index1][index2]+=sfock[ha][a][b];
                    }
                  }
                  if(A==B)
                  {
                    if(hi==hb)
                    {
                      matrix[index1][index2]-=sfock[hi][i][j];
                    }
                  }
                  matrix[index1][index2]+=so[A][J][I][B];
                  index2++;
                }
              }  
            }
          }
          index1++;
        }
      } 
    }
  }
}

void calcsinglets(double **singlet,double ***sfock,double *tei,int nirreps,int *sopi,int *filled,int *symmo,int occupied,int unoccupied, double *eval,double **evec)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project12

  int jb,ai,ab,ji,ajib,ajbi;
  int I,J,A,B;
  int index1=0,index2=0;
  for(int hi=0;hi<nirreps;hi++)
  {
    for(int i=0;i<filled[hi];i++)
    {
      I=i+symmo[hi];
      for(int ha=0;ha<nirreps;ha++)
      {
        for(int a=filled[ha];a<sopi[ha];a++)
        {
          A=a+symmo[ha];
          index2=0;
          for(int hj=0;hj<nirreps;hj++)
          {
            for(int j=0;j<filled[hj];j++)
            {
              J=j+symmo[hj];
              for(int hb=0;hb<nirreps;hb++)
              {
                for(int b=filled[hb];b<sopi[hb];b++)
                {
                  B=b+symmo[hb];          
                  if(I==J)
                  {
                    if(hi==hj)
                    {
                      if(ha==hb)
                      {
                        singlet[index1][index2] += sfock[ha][2*a][2*b];
                      }
                    }
                  }
                  if(A==B)
                  {
                    if(ha==hb)
                    {
                      if(hi==hj)
                      {
                        singlet[index1][index2] -= sfock[hi][2*i][2*j];
                      }
                    }
                  }
                  ai=INDEX(A,I);
                  jb=INDEX(J,B);
                  ajib=INDEX(ai,jb);
                  ab=INDEX(A,B);
                  ji=INDEX(J,I);
                  ajbi=INDEX(ab,ji);  
                  singlet[index1][index2] += 2*tei[ajib]-tei[ajbi];
                  index2++;
                }
              }  
            }
          }
          index1++;
        }
      } 
    }
  }
  sq_rsp( (occupied/2)*(unoccupied/2),(occupied/2)*(unoccupied/2),singlet,eval,1,evec,1e-13);
}
void calctriplets(double **triplet,double ***sfock,double *tei,int nirreps,int *sopi,int *filled,int *symmo,int occupied,int unoccupied, double *eval,double **evec)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project12

  int index1=0;int index2=0;
  int ab,ji,ajbi; 
  int I,J,A,B;
  for(int hi=0;hi<nirreps;hi++)
  {
    for(int i=0;i<filled[hi];i++)
    {
      I=i+symmo[hi];
      for(int ha=0;ha<nirreps;ha++)
      {
        for(int a=filled[ha];a<sopi[ha];a++)
        {
          A=a+symmo[ha];
          index2=0;
          for(int hj=0;hj<nirreps;hj++)
          {
            for(int j=0;j<filled[hj];j++)
            {
              J=j+symmo[hj];
              for(int hb=0;hb<nirreps;hb++)
              {
                for(int b=filled[hb];b<sopi[hb];b++)
                {
                  B=b+symmo[hb];          
                  if(I==J)
                  {
                    if(ha==hb)
                    {
                      if(hi==hj)
                      { 
                        triplet[index1][index2]+=sfock[ha][2*a][2*b];
                      }                      
                    }
                  }
                  if(A==B)
                  {
                    if(ha==hb)
                    {
                      if(hi==hj)
                      {
                        triplet[index1][index2]-=sfock[hi][2*i][2*j];
                      }
                    }
                  }
                  ab = INDEX(A,B);
                  ji = INDEX(J,I);
                  ajbi = INDEX(ab,ji);
                  triplet[index1][index2] -= tei[ajbi];
                  index2++;
                }
              }  
            }
          }
          index1++;
        }
      }      
    }  
  }   
  sq_rsp((occupied/2)*(unoccupied/2),(occupied/2)*(unoccupied/2),triplet,eval,1,evec,1e-13);
}

void buildA(double **matrix,double **cis,int size)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project12
  for(int i=0;i<size;i++)
  {
    for(int j=0;j<size;j++)
    {
      matrix[i][j] = cis[i][j];
    }
  }
}
void buildB(double **matrix,double ****so,int nirreps,int *sopi,int *filled,int *symmo)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project12
  int index1=0, index2=0;
  int I,J,A,B;
  for(int hi=0;hi<nirreps;hi++)
  {
    for(int i=0;i<2*filled[hi];i++)
    {
      I=i+2*symmo[hi];
      for(int ha=0;ha<nirreps;ha++)
      {
        for(int a=2*filled[ha];a<2*sopi[ha];a++)
        {
          A=a+2*symmo[ha];         
          index2=0;
          for(int hj=0;hj<nirreps;hj++)
          {
            for(int j=0;j<2*filled[hj];j++)
            {
              J=j+2*symmo[hj];
              for(int hb=0;hb<nirreps;hb++)
              {
                for(int b=2*filled[hb];b<2*sopi[hb];b++)
                {
                   B=b+2*symmo[hb];
                   matrix[index1][index2] = so[A][B][I][J];
                   index2++; 
                }
              }
            }
          }
          index1++;
        }
      }
    }
  }
}

void buildTDHF(double **matrix,double **A,double **B,int size)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project12
  for(int i=0;i<size;i++)
  {
    for(int j=0;j<size;j++)
    {
      matrix[i][j] = A[i][j];
    }
  }
  for(int i=size;i<2*size;i++)
  {
    for(int j=size;j<2*size;j++)
    {
      matrix[i][j] = -A[i-size][j-size];
    }
  }
  for(int i=size;i<2*size;i++)
  {
    for(int j=0;j<size;j++)
    {
      matrix[i][j] = -B[i-size][j];
    }
  }
  for(int i=0;i<size;i++)
  {
    for(int j=size;j<2*size;j++)
    {
      matrix[i][j] = B[i][j-size];
    }
  }
}
void setupTDHF(double **A,double **B,double **AplusB,double **AminusB,double **ans,int occupied,int unoccupied)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project12
  for(int i=0;i<occupied*unoccupied;i++)
  {
    for(int j=0;j<occupied*unoccupied;j++)
    {
      AplusB[i][j] = A[i][j] + B[i][j];
      AminusB[i][j] = A[i][j] - B[i][j];
    }
  }
  mmult(AplusB,0,AminusB,0,ans,0,occupied*unoccupied,occupied*unoccupied,occupied*unoccupied,0);
}

void comparison(double **ab,double **comp, int initialL,int N,int L,int *location,double *locVal)
{
  //Identifies and flags position of desired eigenvectors for Davidson Lui Algorithm
  bool checker=1;
  double diffsum,plussum;
  for(int i=0;i<initialL;i++)
  {
    locVal[i]=100000000.0;
    location[i]=-1;
    for(int j=0;j<L;j++)
    {
      diffsum=0.0;
      plussum=0.0;
      for(int k=0;k<N;k++)
      {
        diffsum += pow( (ab[j][k] - comp[k][i]),2 );
        plussum += pow( (ab[j][k] + comp[k][i]),2 );
      }
      diffsum = sqrt(diffsum);
      plussum = sqrt(plussum);
      for(int k=0;k<i;k++)
      {
        if(location[k]==j)
        {
          checker=0;
        }
      }
      if( ((diffsum<locVal[i]) || (plussum < locVal[i])) && checker )
      {
        if(diffsum<plussum)
        {
          locVal[i]=diffsum;
          location[i]=j;
        }
        else
        {
          locVal[i]=plussum;
          location[i]=j;
        }
      }
      checker=1;
    }
  }
}
void gschmidt(double **matrix,int rows,int columns)
{
  //http://en.wikipedia.org/wiki/Gram_schmidt

  double **u;
  double **v;
  v=block_matrix(columns,rows);
  u=block_matrix(columns,rows);
  zero_mat(v,columns,rows);  
  zero_mat(u,columns,rows);
  double *rtrn;
  rtrn=init_array(rows);
  for(int i=0;i<rows;i++)
  {
    for(int j=0;j<columns;j++)
    {
      v[j][i]=matrix[i][j];
    }
  }
  double col;
  double rw;
  col=rows;
  rw=columns;
  for(int i=0;i<rw;i++)
  {
    for(int j=0;j<col;j++)
    {
      u[i][j]=v[i][j];
    }
  }
  for(int i=0;i<rw;i++)
  {
    for(int j=0;j<i;j++)
    {
      if(i>0)
      {
        proj(u[j],v[i],col,rtrn);
        for(int k=0;k<col;k++)
        {
          u[i][k] -= rtrn[k];
        }
      }
    }
  }
  for(int i=0;i<rw;i++)
  {
    for(int j=0;j<col;j++)
    {
      v[i][j]=u[i][j];
    }
  }
  for(int i=0;i<rows;i++)
  {
    for(int j=0;j<columns;j++)
    {
      matrix[i][j]=v[j][i];
    }
  }
  free(rtrn);
  free_block(u);
  free_block(v);
}
void proj(double *u,double *v,int n,double *rtrn)
{
  //http://en.wikipedia.org/wiki/Gram_schmidt

  double vu=0.0;
  double uu=0.0;
  double scale;
  for(int i=0;i<n;i++)
  {
    rtrn[i]=u[i];
    vu+=v[i]*u[i];
    uu+=u[i]*u[i];
  }
  if(uu>0.0001)
  {
    scale=vu/uu;
  }
  else
  {
    scale=vu;
  }
  for(int i=0;i<n;i++)
  {
    rtrn[i]=scale*u[i];
  }  
}
void davidsonlui(int occupied,int unoccupied,double **singlets,int initialL,double *dlvals,int high)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project13
    
  int *pos;
  double **guess;
  double **G;
  double **temporary;
  double *real1,*imag1;
  double **left1,**right1;
  double **workspace1;
  double **correction;
  double **residuals;
  double **sigma;
  double **alphab;
  double *norm;
  bool *include;
  double **modguess;
  double **saveab;
  int *location;
  double *locVal;
  double *prevval;


  int L=initialL;
  int prevL;
  int N=(occupied/2)*(unoccupied/2);
  char decider1;
  char decider2;
  int inc;
  int sz;
  int worksize;
  int info;
  saveab = block_matrix(N,initialL);
  prevval = init_array(initialL);

  int count=0;
  int counter=0,counter1=0,counter2=0;

  while( (counter1!=initialL  || counter2!=initialL) )
  {
    if(count>0)
    {
      for(int i=0;i<initialL;i++)
      {
        prevval[i]=real1[location[i]];
      }
      free(real1);
    }

    pos = (int *) malloc(L*sizeof(int));
    guess = block_matrix(N,L);
    zero_mat(guess,N,L);
    G = block_matrix(L,L);
    zero_mat(G,L,L);
    temporary = block_matrix(N,L);
    zero_mat(temporary,N,L);
    real1=init_array(L);
    imag1=init_array(L);
    left1=block_matrix(L,L);
    zero_mat(left1,L,L);
    right1=block_matrix(L,L);
    zero_mat(right1,L,L);
    workspace1=block_matrix(10*L,10*L);
    zero_mat(workspace1,10*L,10*L);
    correction = block_matrix(L,N);
    zero_mat(correction,L,N);
    residuals = block_matrix(L,N);
    zero_mat(residuals,L,N);
    sigma = block_matrix(L,N);
    zero_mat(sigma,L,N);
    alphab = block_matrix(L,N);
    zero_mat(alphab,L,N);
    norm = init_array(initialL);
    include = (bool *) malloc(initialL*sizeof(bool));      
    location = init_int_array(initialL);
    locVal = init_array(initialL);
      
    if(count==0)
    {
      double *holder;
      int *hloc;
      holder=init_array((occupied/2)*(unoccupied/2));
      hloc = init_int_array((occupied/2)*(unoccupied/2));
      for(int i=0;i<(occupied/2)*(unoccupied/2);i++)
      {
        holder[i] = singlets[i][i];
        hloc[i]=i;
      }
      int tempint;
      double tempdouble;
      bool hflag=1;  
      for(int i=1;i<=(occupied/2)*(unoccupied/2) && hflag;i++)
      {
        hflag=0;
        for(int j=0;j<( (occupied/2)*(unoccupied/2)-1 );j++)
        {
          if(holder[j+1]>holder[j])
          {
            tempdouble=holder[j];
            tempint=hloc[j];
            holder[j]=holder[j+1];
            hloc[j]=hloc[j+1];
            holder[j+1]=tempdouble;
            hloc[j+1]=tempint;
            hflag=1;
          }
        }
      }
      if(high==1)
      {
        for(int i=0;i<initialL;i++)
        {
          guess[hloc[i]][i]=1.0;
        }
      }
      else
      {
        int *reverseloc;
        reverseloc=init_int_array((occupied/2)*(unoccupied/2));
        for(int i=0;i<(occupied/2)*(unoccupied/2);i++)
        {
          reverseloc[i]=hloc[((occupied/2)*(unoccupied/2))-i-1];
        }
        for(int i=0;i<initialL;i++)
        {
          guess[hloc[(occupied/2)*(unoccupied/2)-i-1]][i]=1.0;
        }
        free(reverseloc);
      }
      free(hloc);
      free(holder);
      for(int i=0;i<N;i++)
      {
        for(int j=0;j<L;j++)
        {
          saveab[i][j] = guess[i][j];
        }
      }  
    }  
    else
    {        
      for(int i=0;i<N;i++)
      {
        for(int j=0;j<L;j++)
        {
          guess[i][j]=modguess[i][j];
        }
      }
      free_block(modguess);
    }
    
    mmult(singlets,0,guess,0,temporary,0,N,N,L,0);
    mmult(guess,1,temporary,0,G,0,L,N,L,0);
    
    decider1='n';
    decider2='v';
    sz=L;
    worksize=10*L;
    dgeev_(&decider1,&decider2,&sz,&G[0][0],&sz,&real1[0],&imag1[0],&left1[0][0],&sz,&right1[0][0],&sz,&workspace1[0][0],&worksize,&info);

    mmult(singlets,0,guess,0,sigma,1,N,N,L,0);
    mmult(right1,0,guess,1,alphab,0,L,L,N,0);
    mmult(right1,0,sigma,0,residuals,0,L,L,N,0);
    comparison(alphab,saveab,initialL,N,L,location,locVal);

    for(int i=0;i<initialL;i++)
    {
      for(int j=0;j<N;j++)
      {
        saveab[j][i]=alphab[location[i]][j];        
      }
    }

    //subspace collapse
    if(L>N)
    {
      free_block(guess);
      L=initialL;
      guess=block_matrix(N,L);
      for(int i=0;i<N;i++)
      {
        for(int j=0;j<initialL;j++)
        {
          guess[i][j] = alphab[location[j]][i];
        }
      }
      free_block(sigma);
      sigma = block_matrix(L,N);
      mmult(singlets,0,guess,0,sigma,1,N,N,L,0);
      free_block(alphab);
      alphab = block_matrix(L,N);
      mmult(right1,0,guess,1,alphab,0,L,L,N,0);
      free_block(residuals);
      residuals = block_matrix(L,N);
      mmult(right1,0,sigma,0,residuals,0,L,L,N,0);
      comparison(alphab,saveab,initialL,N,L,location,locVal);
      free_block(correction);
      correction=block_matrix(L,N);
    }
 
    for(int i=0;i<L;i++)
    {
      C_DAXPY(N,-real1[i],alphab[i],1,residuals[i],1);
    }
    for(int k=0;k<L;k++)
    { 
      for(int I=0;I<N;I++)
      {
      if(fabs(real1[k]-singlets[I][I]) > 0.00001 )
        {
          correction[k][I] =pow( (real1[k]-singlets[I][I]),-1 )*residuals[k][I];
        }
        else
        {
          correction[k][I] = 0.0;
        }
      }  
    }
    double sum;
    for(int i=0;i<L;i++)
    {
      sum=0.0;
      for(int j=0;j<N;j++)
      { 
        sum+=correction[i][j]*correction[i][j];
      }
      sum=sqrt(sum);
      for(int j=0;j<N;j++)
      {
        if(sum > 0.00001 )
        {
          correction[i][j] /= sum;
        }
      }
    }
    for(int i=0;i<initialL;i++)
    {
      int I;
      I=location[i];
      for(int j=0;j<N;j++)
      {
        norm[i] +=( correction[I][j] * correction[I][j] );
      }
      norm[i]=sqrt(norm[i]);
    }

    inc=0;
    for(int i=0;i<initialL;i++)
    {
      include[i]=0;
      if(norm[i]>1e-3)
      {
        include[i]=1;
        inc++;
      }
    }
    prevL=L;
    L+=inc;
    modguess = block_matrix(N,L);
    zero_mat(modguess,N,L);
    for(int i=0;i<prevL;i++)
    {
      for(int j=0;j<N;j++)
      {
        modguess[j][i]=guess[j][i];
      }
    }
    counter=prevL;
    for(int i=prevL;i<prevL+initialL;i++)
    {
      if(include[i-prevL]==1)
      {
        for(int j=0;j<N;j++)
        {
          modguess[j][counter] = correction[location[i-prevL]][j];
        }
        counter++;
      }
    }
    gschmidt(modguess,N,L);
    for(int i=0;i<L;i++)
    {
      sum=0.0;
      for(int j=0;j<N;j++)
      {
        sum+=modguess[j][i]*modguess[j][i];
      }
      sum=sqrt(sum);
      for(int j=0;j<N;j++)
      {
        if(fabs(sum)> 1e-3)
        { 
          modguess[j][i] /= sum;
        }
      }
    }

    counter1=0;
    counter2=0;
    if(count>0)
    {
      for(int i=0;i<initialL;i++)
      {
        if(fabs(locVal[i])<1e-15)
        {
          counter1++;
        }
        if(fabs( prevval[i]-real1[location[i]] ) < 1e-15)
        {
          counter2++;
        }
      }
    }
    for(int i=0;i<initialL;i++)
    {
      dlvals[i]=real1[location[i]];
    }

    free(pos);
    free_block(guess);
    free_block(G);
    free_block(temporary);
    free(imag1);
    free_block(left1);
    free_block(right1);
    free_block(workspace1);
    free_block(correction);
    free_block(residuals);
    free_block(sigma);
    free_block(alphab);
    free(norm);
    free(include);

    count++;
  }
  free(real1);
  free(location);
  free(locVal);
  free(prevval);
  free_block(modguess);
  free_block(saveab);
}

void clockPrint(std::string title)
{
  std::string mix=title+"(Use) %ld %ld :(Sys) %ld %ld\n";
  struct rusage now;
  getrusage(RUSAGE_SELF,&now);
  timeval timeUCPU=now.ru_utime;
  timeval timeSCPU=now.ru_stime;
  fprintf(stderr,mix.c_str(),timeUCPU.tv_sec,timeUCPU.tv_usec,timeSCPU.tv_sec,timeSCPU.tv_usec);
}
void findfilled(double *zval,double **e_vals,int natom,int nirreps,int *sopi,int &occupied,int &unoccupied,int *filled)
{
  int sum=0;
  for(int i=0;i<nirreps;i++)
  {
    sum+=sopi[i];
  }
  double *fockval;
  fockval=init_array(sum);
  int *loc;
  loc = init_int_array(sum);

  int count=0;
  for(int h=0;h<nirreps;h++)
  {
    for(int i=0;i<sopi[h];i++)
    {
      fockval[count]=e_vals[h][i];
      loc[count]=h;
      count++;
    }
  }
  bool flag=1;
  int inttemp;
  double doubletemp;
  for(int i=0;i<=sum && flag;i++)
  {
    flag=0;
    for(int j=0;j<sum-1;j++)
    {
      if(fockval[j+1]<fockval[j])
      {
        doubletemp=fockval[j];
        inttemp=loc[j];
        fockval[j]=fockval[j+1];
        loc[j]=loc[j+1];
        fockval[j+1]=doubletemp;
        loc[j+1]=inttemp;
        flag=1;
      }    
    }
  }
  occupied=0;
  for(int i=0;i<natom;i++)
  {
    occupied += (int) zval[i];
  }
  occupied /= 2;
  for(int i=0;i<nirreps;i++)
  {
    filled[i]=0;
  }
  for(int i=0;i<occupied;i++)
  {
    filled[loc[i]]+=1;
  }
  for(int i=0;i<nirreps;i++)
  {
    unoccupied += sopi[i];
  }
  occupied *= 2;
  unoccupied *=2;
  unoccupied -= occupied;

  free(fockval);
  free(loc);
}

double dist(double x1, double x2, double y1, double y2, double z1, double z2)
{
  double distance;
  distance = sqrt( ( pow( (x1-x2),2 ) + pow( (y1-y2),2 ) + pow( (z1-z2),2 )  ) );
  return distance;
}

void moleculedata(int atomno,int *zno,double **geom)
{
  //http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project1

  double *xcoord,*ycoord,*zcoord;
  double **dmatrix;
  double **ex,**ey,**ez;
  double ***phi;
  double theta;
  double tau;
  double ejkl_x, ejkl_y, ejkl_z, exx, eyy, ezz;
  double eijk_x,eijk_y,eijk_z;
  double cross_x, cross_y, cross_z, norm, dot, sign;
  double xcm=0.0, ycm=0.0, zcm=0.0, mi=0.0 ,M=0.0;
  double conv;
  
  xcoord=init_array(atomno);
  ycoord=init_array(atomno);
  zcoord=init_array(atomno);
  for(int i=0;i<atomno;i++)
  {
    zcoord[i]=geom[i][0];
  }
  for(int i=0;i<atomno;i++)
  {
    xcoord[i]=geom[i][1];
  }
  for(int i=0;i<atomno;i++)
  {
    ycoord[i]=geom[i][2];
  }

  //Allocate memory for dmatrix
  dmatrix=(double **) malloc(atomno*sizeof(double*));
  for(int i=0;i<atomno;i++)
  {
    dmatrix[i]=(double *) malloc(atomno*sizeof(double));
  }

  //Allocate memory for unit vectors and bond angle
  //Makes  1-D array and allocates memory--sizeof(double) allocates enough for doubles, atomno tells number of doubles to save for
  ex = (double **) malloc(atomno * sizeof(double));  
  ey = (double **) malloc(atomno * sizeof(double));
  ez = (double **) malloc(atomno * sizeof(double));
  //Allocates memory for second dimension in array
  for(int i=0;i<atomno;i++)
  {
    ex[i]=(double*) malloc(atomno * sizeof(double));
    ey[i]=(double*) malloc(atomno * sizeof(double));
    ez[i]=(double*) malloc(atomno * sizeof(double));
  }
  //Allocate memory for 3-D array for bond angle (phi)
  phi = (double ***) malloc(atomno * sizeof(double **));
  for(int i=0; i<atomno; i++)
  {
    phi[i]=(double **) malloc(atomno * sizeof(double *));
    for(int j=0; j<atomno; j++)
    {
      phi[i][j] = (double *) malloc(atomno * sizeof(double));
    }
  }
	  
  //Fill Distance Matrix
  for(int i=0;i<atomno;i++)
  {
    for(int j=0;j<atomno;j++)
    {
      dmatrix[i][j]=dist(xcoord[i],xcoord[j],ycoord[i],ycoord[j],zcoord[i],zcoord[j]);
    }
  }

  //Calculating interatomic unit vectors
  for(int i=0;i<atomno;i++)
  {
    for(int j=0;j<atomno;j++)
    {
      if(i != j)
      {
        ex[i][j] = -( xcoord[i]-xcoord[j] )/dmatrix[i][j];
	ey[i][j] = -( ycoord[i]-ycoord[j] )/dmatrix[i][j];
	ez[i][j] = -( zcoord[i]-zcoord[j] )/dmatrix[i][j];
      }
    }		
  }
  
  //Calculating bond angles
  for(int i=0; i<atomno; i++)
  {
    for(int j=0; j<atomno; j++)
    {
      for(int k=0; k<atomno; k++)
      {
        if(i!=j && i!=k && j!=k)
        {
          phi[i][j][k]=acos( (ex[j][i]*ex[j][k]) + (ey[j][i]*ey[j][k]) + (ez[j][i]*ez[j][k]) );
        }
	else
	{
	  phi[i][j][k]=0.0;
	}
      }
    }
  }

  //Display
  int counter=0;
  std::cout << "Number of atoms: " << atomno << std::endl;
  std::cout << "Coordinates:" << std::endl;
  std::cout << "zno" << std::setw(10) << "x coord" << std::setw(12) << "y coord" << std::setw(12) << "z coord" << std::endl;
  for(int i=0;i<atomno;i++)
  {
    std::cout << zno[i] << std::setw(12) << std::setprecision(7) << xcoord[i] << std::setw(12) << std::setprecision(7) << ycoord[i] << std::setw(12) <<std::setprecision(7) << zcoord[i] << std::endl;
  }
  std::cout << std::endl << "Interatomic Distances (bohr)" << std::endl;
  for(int i=0;i<atomno;i++)
  {
    for(int j=0;j<atomno;j++)
    {
      if(i != j)
      {
        std::cout << i << std::setw(5) << j << std::setw(10) << std::setprecision(7) << dmatrix[i][j] << std::endl;
      }
    }
  }
  std::cout << std::endl << "Bond Angles:" << std::endl;
  for(int i=0;i<atomno;i++)
  {
    for(int j=0; j<atomno;j++)
    {
      for(int k=0;k<atomno;k++)
      {
        if(i<j && j<k && dmatrix[i][j] < 4.0 && dmatrix[j][k] < 4.0)
	{
	  std::cout << i << "-" << j << "-" << k << std::setw(10) << phi[i][j][k]*(180/acos(-1.0)) << std::endl;
          counter++;
	}
      }
    }
  }
  if(counter==0)
  {
    std::cout << "None" << std::endl;
  }
  //Calculate and Display Out-of-Plane Angles
  std::cout << std::endl << "Out-of-Plane Angles: " << std::endl;
  counter=0;
  for(int i=0;i<atomno;i++)
  {
    for(int k=0;k<atomno;k++)
    {
      for(int j=0;j<atomno;j++)
      {
        for(int l=0;l<j;l++)
        {
          //Computing cross product components
	  ejkl_x = (ey[k][j] * ez[k][l] - ez[k][j] * ey[k][l]);
	  ejkl_y = (ez[k][j] * ex[k][l] - ex[k][j] * ez[k][l]);
	  ejkl_z = (ex[k][j] * ey[k][l] - ey[k][j] * ex[k][l]);
          //Dot product with e[k][i]
	  exx = ejkl_x * ex[k][i];
	  eyy = ejkl_y * ey[k][i];
	  ezz = ejkl_z * ez[k][i];

	  theta = (exx + eyy + ezz) / sin( phi[j][k][l] );

	  if(theta < -1.0)
	  {
	    theta = asin(-1.00);
	  }
	  else if(theta > 1.0)
	  {
	    theta = asin(1.00);
	  }
	  else
	  {
	    theta = asin(theta);
	  }
          				  
	  if(i != j && i != k && i != l && j != k && j != l && k != l && dmatrix[i][k] < 4.0 && dmatrix[j][k] < 4.0 && dmatrix[k][l] < 4.0)
	  {
	    std::cout << i << "- " << j << "- " << k << "- " << l << std::setw(12) << std::setprecision(6) << theta*(180/acos(-1.0)) << std::endl;
            counter++;
	  }
        }
      }
    }
  }
  if(counter==0)
  {
    std::cout << "None" << std::endl;
  }
  //Calculate and Display Torsional/Dihedral Angles
  std::cout << std::endl << "Torsional / Dihedral Angles:" << std::endl;
  counter=0;
  for(int i=0;i<atomno;i++)
  {
    for(int j=0;j<i;j++)
    {
      for(int k=0;k<j;k++)
      {
        for(int l=0;l<k;l++)
	{
  	  eijk_x = (ey[j][i] * ez[j][k] - ez[j][i] * ey[j][k]);
	  eijk_y = (ez[j][i] * ex[j][k] - ex[j][i] * ez[j][k]);
	  eijk_z = (ex[j][i] * ey[j][k] - ey[j][i] * ex[j][k]);

	  ejkl_x = (ey[k][j] * ez[k][l] - ez[k][j] * ey[k][l]);
	  ejkl_y = (ez[k][j] * ex[k][l] - ex[k][j] * ez[k][l]);
	  ejkl_z = (ex[k][j] * ey[k][l] - ey[k][j] * ex[k][l]);

          exx = eijk_x * ejkl_x;
	  eyy = eijk_y * ejkl_y;
	  ezz = eijk_z * ejkl_z;

	  tau = (exx + eyy + ezz) / ( sin(phi[i][j][k])*sin(phi[j][k][l]) );
	  if(tau < -1.0)
	  {
	    tau = acos(-1.0);
	  }
	  else if(tau > 1.0)
	  {
	    tau = acos(1.0);
	  }
          else
	  {
	    tau = acos(tau);	
	  }
	  //Check the sign
	  cross_x = eijk_y * ejkl_z - eijk_z * ejkl_y;
	  cross_y = eijk_z * ejkl_x - eijk_x * ejkl_z;
	  cross_z = eijk_x * ejkl_y - eijk_y * ejkl_x;
	  norm = (cross_x * cross_x)+(cross_y * cross_y)+(cross_z * cross_z); //find normal vector for normalization
	  cross_x /= norm;  //normalize
	  cross_y /= norm;  //normalize
	  cross_z /= norm;  //normalize
	  sign = 1.0;
	  dot = cross_x*ex[j][k] + cross_y*ey[j][k]+cross_z*ez[j][k];
	  if(dot < 0.0)
	  {
	    sign = -1.0;
	  }
	  if(j<i && k<j && l<k && dmatrix[i][j] < 4.0 && dmatrix[j][k] < 4.0 && dmatrix[k][l] < 4.0)
	  {
	    std::cout << i << "- " << j << "- " << k << "- " << l << std::setw(12) << std::setprecision(6) << sign*tau*(180/acos(-1.0)) << std::endl; 
            counter++;
	  }
	}
      }
    }
  }
  if(counter==0)
  {
    std::cout << "None" << std::endl;
  }
  //Calculate Center-of-Mass
  for(int i=0; i<atomno; i++)
  {
    M = M + masses[int(zno[i])];
  }

  for(int i=0; i<atomno;i++)
  {
    mi=masses[int(zno[i])];
    xcm += mi*xcoord[i];
    ycm += mi*ycoord[i];
    zcm += mi*zcoord[i];
  }
  xcm /= M;
  ycm /= M;
  zcm /= M;
  std::cout << std::endl << "Center of Mass:" << std::endl;
  std::cout << xcm << std::setw(15) << ycm << std::setw(5) << zcm << std::endl; 

  // shift the molecule to the COM
  for(int i=0; i < atomno; i++) 
  {
    xcoord[i] -= xcm; 
    ycoord[i] -= ycm; 
    zcoord[i] -= zcm;
  }

  double **inertia;
  //allocate memory for inertia matrix
  inertia = (double **) malloc(3*sizeof(double *));
  for(int i=0;i<3;i++)
  {
    inertia[i]= (double *) malloc(3*sizeof(double));
  }
  //initialize inertia matrix to zero
  for(int i=0;i<3;i++)
  {
    for(int j=0;j<3;j++)
    {
      inertia[i][j]=0.0;
    }
  }

  //allocate memory for e_vals and e_vecs
  double *e_vals, **e_vecs;
  e_vals=(double *) malloc(3*sizeof(double));
  e_vecs=(double**) malloc(3*sizeof(double *));
  for(int i=0;i<3;i++)
  {
    e_vecs[i]=(double*) malloc(3*sizeof(double));
  }

  for(int i=0;i<atomno;i++) 
  {
    mi = masses[int(zno[i])];
    inertia[0][0] += mi * (ycoord[i]*ycoord[i] + zcoord[i]*zcoord[i]);
    inertia[1][1] += mi * (xcoord[i]*xcoord[i] + zcoord[i]*zcoord[i]);
    inertia[2][2] += mi * (xcoord[i]*xcoord[i] + ycoord[i]*ycoord[i]);
    inertia[0][1] += mi * xcoord[i] * ycoord[i];
    inertia[0][2] += mi * xcoord[i] * zcoord[i];
    inertia[1][2] += mi * ycoord[i] * zcoord[i];
  }
  inertia[1][0] = inertia[0][1];
  inertia[2][0] = inertia[0][2];
  inertia[2][1] = inertia[1][2];

  sq_rsp(3,3,inertia,e_vals,0,e_vecs,1e-13);

  std::cout << std::endl << "Moment of Inertia Tensor (amu bohr^2):" << std::endl;
  for(int i=0;i<3;i++)
  {
    for(int j=0;j<3;j++)
    {
      std::cout << std::setprecision(5)<< std::setw(10) << inertia[i][j];
    }
    std::cout << std::endl;
  }
  conv = 0.529177249 * 0.529177249;

  std::cout << std::endl << "Principle Moments of Inertia (amu bohr^2):" << std::endl;
  for(int i=0; i<3;i++)
  {
    std::cout << std::setprecision(5) << e_vals[i] << std::setw(8);
  }
  std::cout << std::endl;


  std::cout << std::endl << "Principle Moments of Inertia (amu AA^2):" << std::endl;
  for(int i=0; i<3;i++)
  {
    std::cout << std::setprecision(5) << e_vals[i]*conv << std::setw(8);
  }
  std::cout << std::endl;

  conv = 1.6605402E-24 * 0.529177249E-8 * 0.529177249E-8;
  std::cout << std::endl << "Principle Moments of Inertia (g*cm^2):" << std::endl;
  for(int i=0; i<3;i++)
  {
    std::cout << std::setprecision(5) << std::setw(12) << e_vals[i]*conv;
  }
  std::cout << std::endl <<std::endl;

  //classify
  if(atomno==2)
  {
    std::cout << "The molecule is diatomic." << std::endl;
  }
  else if( (abs(e_vals[0]-e_vals[1]) < 1e-4) && (abs(e_vals[1]-e_vals[2])<1e-4))
  {
    std::cout << "The molecule is a spherical top." << std::endl;
  }
  else if((abs(e_vals[0] - e_vals[1]) > 1e-4) && (abs(e_vals[1] - e_vals[2]) < 1e-4))
  {
    std::cout << "The molecule is a prolate spherical top"<< std::endl;
  }
  else
  {
    std::cout << "The molecule is an asymmetric top." << std::endl;
  }
  std::cout << std::endl;

  conv = 6.6260755E-34/(8.0 * acos(-1.0) * acos(-1.0));
  conv /= 1.6605402E-27 * 0.529177249E-10 * 0.529177249E-10;
  conv *= 1e-6;
  std::cout << "Rotational constants (MHz):" << std::endl;
  for(int i=0; i<3;i++)
  {
    std::cout << std::setprecision(5) << std::setw(12) <<conv/e_vals[i];
  }
  std::cout << std::endl << std::endl;

  conv = 6.6260755E-34/(8.0 * acos(-1.0) *acos(-1.0));
  conv /= 1.6605402E-27 * 0.529177249E-10 * 0.529177249E-10;
  conv /= 2.99792458E10;

  std::cout << "Rotational constants (cm-1):" << std::endl;
  for(int i=0; i<3;i++)
  {
    std::cout << std::setprecision(5) << std::setw(12)<< conv/e_vals[i];
  }
  std::cout << std::endl<<std::endl;;

  for(int i=0;i<3;i++)
  {
    free(e_vecs[i]);
    free(inertia[i]);
  }
  free(e_vecs);
  free(e_vals);
  free(inertia);
  for(int i=0;i<atomno;i++)
  {
    for(int j=0;j<atomno;j++)
    {
      free(phi[i][j]);
    }
    free(phi[i]);
  }
  free(phi);

  for(int i=0;i<atomno;i++)
  {
    free(dmatrix[i]);
  }
  free(dmatrix);
  free(xcoord);
  free(ycoord);
  free(zcoord);
}
