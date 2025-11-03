#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "tools.h"

/*-----------------------------------------------------------
/
/   lchk - license check
/
/------------------------------------------------------------*/
long lchk(long lic)
{
  long       ok;
  long      i, test1, test2;
  time_t    now;

  test1=130;
  test2=190;
  i=10000000;
  test1=test1*i;
  test2=test2*i;
  now=time(NULL);
  ok=0;
  if(now>test1 && now<test2) {
    ok=1;
  }
if(lic==1) {
  printf("\nControl value = %ld \n ",now);
  printf("Beg value = %s ",ctime(&test1));
  printf("Now value = %s ",ctime(&now));
  printf("End value = %s \n",ctime(&test2));
}
  return ok;
}

/*-----------------------------------------------------------
/
/   rmsd - RMSD value calculated on selected atoms
/
/------------------------------------------------------------*/
float rmsd(long n_atoms, long atoms[2][MAXRES+1],
           pdb_struct pdb[2], float r[3][3], float v[3])
{
  long    i,j;
  double XASQ, XBSQ, XNI, T, RTSUM;
  double CMA[4], CMB[4], UMAT[4][4], MAT[4][4];
  float  rms, x0, x1, y0, y1, z0, z1;

/*
  SUBROUTINE TO FIT THE COORD SET XA(N,3) TO THE SET XB(N,3)
  IN THE SENSE OF:
     R*XA + V => XB
  R IS A UNITARY 3.3 RIGHT HANDED ROTATION MATRIX
  AND V IS THE OFFSET VECTOR.

     r(0,0)*xa(0) + r(1,0)*xa(1) + r(2,0)*xa(2) + v(0) => xb(0)
     r(0,1)*xa(0) + r(1,1)*xa(1) + r(2,1)*xa(2) + v(1) => xb(1)
     r(0,2)*xa(0) + r(1,2)*xa(1) + r(2,2)*xa(2) + v(2) => xb(2)
*/

  rms=999.99;
  if(n_atoms<3) return rms;

  RTSUM=0.0;
  XASQ=0.0;
  XBSQ=0.0;
  XNI=1.0/n_atoms;

  for(i=1;i<=3;i++) {
    CMA[i]=0.0;
    CMB[i]=0.0;
    for(j=1;j<=3;j++) {
      MAT[i][j]=0.0;
    }
  }
  for(j=1;j<=n_atoms;j++) {
    x0=pdb[0].atom[atoms[0][j]].R.x;
    x1=pdb[1].atom[atoms[1][j]].R.x;
    y0=pdb[0].atom[atoms[0][j]].R.y;
    y1=pdb[1].atom[atoms[1][j]].R.y;
    z0=pdb[0].atom[atoms[0][j]].R.z;
    z1=pdb[1].atom[atoms[1][j]].R.z;

    MAT[1][1]=MAT[1][1]+x0*x1;
    MAT[1][2]=MAT[1][2]+x0*y1;
    MAT[1][3]=MAT[1][3]+x0*z1;
    MAT[2][1]=MAT[2][1]+y0*x1;
    MAT[2][2]=MAT[2][2]+y0*y1;
    MAT[2][3]=MAT[2][3]+y0*z1;
    MAT[3][1]=MAT[3][1]+z0*x1;
    MAT[3][2]=MAT[3][2]+z0*y1;
    MAT[3][3]=MAT[3][3]+z0*z1;

    XASQ=XASQ+x0*x0+y0*y0+z0*z0;
    XBSQ=XBSQ+x1*x1+y1*y1+z1*z1;
    CMA[1]=CMA[1]+x0;
    CMB[1]=CMB[1]+x1;
    CMA[2]=CMA[2]+y0;
    CMB[2]=CMB[2]+y1;
    CMA[3]=CMA[3]+z0;
    CMB[3]=CMB[3]+z1;
  }
  for(i=1;i<=3;i++) {
    XASQ=XASQ-CMA[i]*CMA[i]*XNI;
    XBSQ=XBSQ-CMB[i]*CMB[i]*XNI;
    for(j=1;j<=3;j++) {
      UMAT[i][j]=(MAT[i][j]-CMA[i]*CMB[j]*XNI)*XNI;
    }
  }

  RTSUM=qkfit(UMAT, MAT);
/*
     CALCULATE rms DEVIATION, ROTATION MATRIX r AND THE OFFSET VECTOR v
*/
  rms=(XASQ+XBSQ)*XNI-2.0*RTSUM;
  if(rms<=0.0) rms=0.0;
  rms=sqrt(rms);

  for(i=0;i<3;i++) {
    r[0][i]=MAT[0][i];
    r[1][i]=MAT[1][i];
    r[2][i]=MAT[2][i];
    T=MAT[0][i]*CMA[1]+MAT[1][i]*CMA[2]+MAT[2][i]*CMA[3];
    v[i]=(CMB[i+1]-T)*XNI;
  }

  return rms;
}

/*
*************************************************************

     THIS SUBROUTINE IS A COMBINATION OF MCLACHLAN'S AND KABSCH'S
     TECHNIQUES. SEE
     KABSCH, W. ACTA CRYST A34, 827,1978
     MCLACHAN, A.D., J. MOL. BIO. NNN, NNNN 1978

*************************************************************
*/
double qkfit(double UMAT[4][4], double R[4][4])
{
  long    i, j, ISIG;
  long    JJ[4], LL[4];
  double B1, B2, B13, B23, B33, S, RTSUM;
  double UTR[7], A[4][4], B[4][4];

  JJ[1]=0; JJ[2]=1; JJ[3]=3;
  LL[1]=1; LL[2]=3; LL[3]=6;

  ISIG=1;
  RTSUM=0.0;
/*
     FORM USQ = (UT)*U   (IN UPPER TRIANGULAR SYMMETRIC STORAGE MODE)
*/
  UTR[1]=UMAT[1][1]*UMAT[1][1]+UMAT[2][1]*UMAT[2][1]+UMAT[3][1]*UMAT[3][1];
  UTR[2]=UMAT[1][1]*UMAT[1][2]+UMAT[2][1]*UMAT[2][2]+UMAT[3][1]*UMAT[3][2];
  UTR[3]=UMAT[1][2]*UMAT[1][2]+UMAT[2][2]*UMAT[2][2]+UMAT[3][2]*UMAT[3][2];
  UTR[4]=UMAT[1][1]*UMAT[1][3]+UMAT[2][1]*UMAT[2][3]+UMAT[3][1]*UMAT[3][3];
  UTR[5]=UMAT[1][2]*UMAT[1][3]+UMAT[2][2]*UMAT[2][3]+UMAT[3][2]*UMAT[3][3];
  UTR[6]=UMAT[1][3]*UMAT[1][3]+UMAT[2][3]*UMAT[2][3]+UMAT[3][3]*UMAT[3][3];
/*
     CALCULATE EIGENVALUES AND VECTORS
*/
  eigen(UTR, A, JJ);
  esort(UTR, A, LL);
/*
     SET A3 = A1 CROSS A2
     ROOTS ARE IN ORDER R(1) >= R(2) >= R(3) >= 0
*/
  A[1][3]=A[2][1]*A[3][2]-A[3][1]*A[2][2];
  A[2][3]=A[3][1]*A[1][2]-A[1][1]*A[3][2];
  A[3][3]=A[1][1]*A[2][2]-A[2][1]*A[1][2];
/*
     VECTOR SET B=U*A
*/
  for(i=1;i<=3;i++) {
    for(j=1;j<=3;j++) {
      B[j][i]=UMAT[j][1]*A[1][i]+UMAT[j][2]*A[2][i]+UMAT[j][3]*A[3][i];
    }
  }
/*
     NORMALIZE B1 AND B2 AND CALCULATE B3 = B1 CROSS B2
*/
  B1=B[1][1]*B[1][1]+B[2][1]*B[2][1]+B[3][1]*B[3][1];
  B2=B[1][2]*B[1][2]+B[2][2]*B[2][2]+B[3][2]*B[3][2];
  if(B1>0.0) B1=1.0/sqrt(B1);
  if(B2>0.0) B2=1.0/sqrt(B2);

  B[1][1]=B[1][1]*B1;
  B[1][2]=B[1][2]*B2;
  B[2][1]=B[2][1]*B1;
  B[2][2]=B[2][2]*B2;
  B[3][1]=B[3][1]*B1;
  B[3][2]=B[3][2]*B2;
/*
     CHECK FOR LEFT HANDED ROTATION
*/
  B13=B[2][1]*B[3][2]-B[3][1]*B[2][2];
  B23=B[3][1]*B[1][2]-B[1][1]*B[3][2];
  B33=B[1][1]*B[2][2]-B[2][1]*B[1][2];

  S=B13*B[1][3]+B23*B[2][3]+B33*B[3][3];
  if(S<0.0) ISIG=-1;
  B[1][3]=B13;
  B[2][3]=B23;
  B[3][3]=B33;
/*
     CALCULATE ROTATION MATRIX r
     CHANGE SIGN OF EVAL #3 IF LEFT HANDED
*/
  for(i=1;i<=3;i++) {
    for(j=1;j<=3;j++) {
      R[i-1][j-1]=B[i][1]*A[j][1]+B[i][2]*A[j][2]+B[i][3]*A[j][3];
    }
  }

  if(UTR[6]>0.0) RTSUM=sqrt(UTR[6]);
  if(ISIG<0) RTSUM=-RTSUM;
  if(UTR[3]>0.0) RTSUM=RTSUM+sqrt(UTR[3]);
  if(UTR[1]>0.0) RTSUM=RTSUM+sqrt(UTR[1]);

  return RTSUM;
}

/*
*************************************************************

     SUBROUTINE TO COMPUTE EIGENVALUES & EIGENVECTORS OF A REAL
     SYMMETRIC MATRIX, STOLEN FROM IBM SSP MANUAL (SEE P165)
     DESCRIPTION OF PARAMETERS:
     A - ORIGINAL MATRIX STORED COLUMNWISE AS UPPER TRIANGLE ONLY,
     I.E. "STORAGE MODE" = 1.  EIGENVALUES ARE WRITTEN INTO DIAGONAL
     ELEMENTS OF A  I.E.  A(1)  A(3)  A(6)  FOR A 3*3 MATRIX.
     R - RESULTANT MATRIX OF EIGENVECTORS STORED COLUMNWISE IN SAME
     ORDER AS EIGENVALUES.
     3 - ORDER OF MATRICES A & R.

*************************************************************
*/
void eigen(double A[7], double R[4][4], long JJ[4])
{
  long    I, J, L, M, IND, LL, LM, LQ, MM, MQ, IQ, IM, IL;
  double RANGE, SINX, SINX2, COSX, COSX2, SINCS, ANORM, THR, X, Y;

  RANGE=0.000000001;
  R[1][1]=1.0; R[1][2]=0.0; R[1][3]=0.0;
  R[2][1]=0.0; R[2][2]=1.0; R[2][3]=0.0;
  R[3][1]=0.0; R[3][2]=0.0; R[3][3]=1.0;
/*
     INITIAL AND FINAL NORM (ANORM)
*/
  ANORM=2.0*(A[2]*A[2]+A[4]*A[4])+A[3]*A[3]+A[5]*A[5];
  if(ANORM<=0.0) return;
/*
     INITIALIZE INDICATORS AND COMPUTE THRESHOLD
*/
  IND=0;
  THR=sqrt(2.0*ANORM);
  ANORM=THR*RANGE/3.0;

  while (THR>ANORM) {
    THR=THR/3.0;

  do {
    IND=0;
  for(L=1;L<=2;L++) {
/*
     COMPUTE SIN & COS
*/
  J=L+1;
  for(M=J;M<=3;M++) {
    MQ=JJ[M];
    LQ=JJ[L];
    LM=L+MQ;

  if(fabs(A[LM])>=THR) {
    IND=1;
    LL=L+LQ;
    MM=M+MQ;
    X=0.5*(A[LL]-A[MM]);
    Y=A[LM]*A[LM]/(A[LM]*A[LM]+X*X);
    SINX2=Y/(2.0*(1.0+sqrt(1.0-Y)));
    SINX=sqrt(SINX2);
    COSX2=1.0-SINX2;
    COSX=sqrt(COSX2);

    if(A[LM]>0.0) SINX=-SINX;
    if(X<0.0) SINX=-SINX;
    SINCS=SINX*COSX;
/*
     ROTATE L & M COLUMNS
*/
    for(I=1;I<=3;I++) {
      if(I!=L && I!=M) {
        IQ=JJ[I];
        if(I<M) IM=I+MQ;
        else IM=M+IQ;
        if(I<L) IL=I+LQ;
        else IL=L+IQ;
        X=A[IL]*COSX-A[IM]*SINX;
        A[IM]=A[IL]*SINX+A[IM]*COSX;
        A[IL]=X;
      }
      X=R[I][L]*COSX-R[I][M]*SINX;
      R[I][M]=R[I][L]*SINX+R[I][M]*COSX;
      R[I][L]=X;
    }
        
    X=2.0*A[LM]*SINCS;
    Y=A[LL];
    A[LM]=(Y-A[MM])*SINCS+A[LM]*(COSX2-SINX2);
    A[LL]=Y*COSX2+A[MM]*SINX2-X;
    A[MM]=Y*SINX2+A[MM]*COSX2+X;
  }
/*
     TESTS FOR COMPLETION
     TEST FOR M = LAST COLUMN
*/
  }
/*
     TEST FOR L = PENULTIMATE COLUMN
*/
  }
/*
     END of the LOOP:
     INDICATOR
*/
  }
  while (IND==1);
/*
     END of the LOOP:
     COMPARE THRESHOLD WITH FINAL NORM
*/
  }

  return;
}

/*
*************************************************************

     SORT EIGENVALUES AND EIGENVECTORS IN DESCENDING ORDER OF
     EIGENVALUES

*************************************************************
*/
void esort(double A[7], double R[4][4], long LL[4])
{
  long    I, J;
  double X;
      
  for(I=1;I<=3;I++) {
    for(J=I;J<=3;J++) {
      if(A[LL[I]]<A[LL[J]]) {
        X=A[LL[I]]; A[LL[I]]=A[LL[J]]; A[LL[J]]=X;
        X=R[1][I]; R[1][I]=R[1][J]; R[1][J]=X;
        X=R[2][I]; R[2][I]=R[2][J]; R[2][J]=X;
        X=R[3][I]; R[3][I]=R[3][J]; R[3][J]=X;
      }
    }
  }

  return;
}

/*-----------------------------------------------------------
/
/   rmsd_any - RMSD value calculated on ANY selected atoms
/
/------------------------------------------------------------*/
float rmsd_any(long n_atoms, long atoms1[MAXATOMS+1], long atoms2[MAXATOMS+1],
               atom_coords coords1[MAXATOMS+1], atom_coords coords2[MAXATOMS+1],
               float r[3][3], float v[3])
{
  long    i,j;
  double XASQ, XBSQ, XNI, T, RTSUM;
  double CMA[4], CMB[4], UMAT[4][4], MAT[4][4];
  float  rms, x0, x1, y0, y1, z0, z1;

/*
  SUBROUTINE TO FIT THE COORD SET XA(N,3) TO THE SET XB(N,3)
  IN THE SENSE OF:
     R*XA + V => XB
  R IS A UNITARY 3.3 RIGHT HANDED ROTATION MATRIX
  AND V IS THE OFFSET VECTOR.

     r(0,0)*xa(0) + r(1,0)*xa(1) + r(2,0)*xa(2) + v(0) => xb(0)
     r(0,1)*xa(0) + r(1,1)*xa(1) + r(2,1)*xa(2) + v(1) => xb(1)
     r(0,2)*xa(0) + r(1,2)*xa(1) + r(2,2)*xa(2) + v(2) => xb(2)
*/

  rms=999.99;
  if(n_atoms<3) return rms;

  RTSUM=0.0;
  XASQ=0.0;
  XBSQ=0.0;
  XNI=1.0/n_atoms;

  for(i=1;i<=3;i++) {
    CMA[i]=0.0;
    CMB[i]=0.0;
    for(j=1;j<=3;j++) {
      MAT[i][j]=0.0;
    }
  }
  for(j=1;j<=n_atoms;j++) {
    x0=coords1[atoms1[j]].R.x;
    x1=coords2[atoms2[j]].R.x;
    y0=coords1[atoms1[j]].R.y;
    y1=coords2[atoms2[j]].R.y;
    z0=coords1[atoms1[j]].R.z;
    z1=coords2[atoms2[j]].R.z;

    MAT[1][1]=MAT[1][1]+x0*x1;
    MAT[1][2]=MAT[1][2]+x0*y1;
    MAT[1][3]=MAT[1][3]+x0*z1;
    MAT[2][1]=MAT[2][1]+y0*x1;
    MAT[2][2]=MAT[2][2]+y0*y1;
    MAT[2][3]=MAT[2][3]+y0*z1;
    MAT[3][1]=MAT[3][1]+z0*x1;
    MAT[3][2]=MAT[3][2]+z0*y1;
    MAT[3][3]=MAT[3][3]+z0*z1;

    XASQ=XASQ+x0*x0+y0*y0+z0*z0;
    XBSQ=XBSQ+x1*x1+y1*y1+z1*z1;
    CMA[1]=CMA[1]+x0;
    CMB[1]=CMB[1]+x1;
    CMA[2]=CMA[2]+y0;
    CMB[2]=CMB[2]+y1;
    CMA[3]=CMA[3]+z0;
    CMB[3]=CMB[3]+z1;
  }
  for(i=1;i<=3;i++) {
    XASQ=XASQ-CMA[i]*CMA[i]*XNI;
    XBSQ=XBSQ-CMB[i]*CMB[i]*XNI;
    for(j=1;j<=3;j++) {
      UMAT[i][j]=(MAT[i][j]-CMA[i]*CMB[j]*XNI)*XNI;
    }
  }

  RTSUM=qkfit(UMAT, MAT);
/*
     CALCULATE rms DEVIATION, ROTATION MATRIX r AND THE OFFSET VECTOR v
*/
  rms=(XASQ+XBSQ)*XNI-2.0*RTSUM;
  if(rms<=0.0) rms=0.0;
  rms=sqrt(rms);

  for(i=0;i<3;i++) {
    r[0][i]=MAT[0][i];
    r[1][i]=MAT[1][i];
    r[2][i]=MAT[2][i];
    T=MAT[0][i]*CMA[1]+MAT[1][i]*CMA[2]+MAT[2][i]*CMA[3];
    v[i]=(CMB[i+1]-T)*XNI;
  }

  return rms;
}

/*-----------------------------------------------------------
/
/   calc_cb - calculate CB's positions 
/
/------------------------------------------------------------*/
void calc_cb(float cb_v, float m_v, float o_v, pdb_struct pdb[2], long in_CB) 
{
  long atoms[2][MAXRES+1];
  long n;
  float r[3][3], v[3];
  float rms=0.0, cb_x, cb_y, cb_z, m_x, m_y, m_z;
  float ca_cb_x, ca_cb_y, ca_cb_z;
  float ca_m_x, ca_m_y, ca_m_z;
  float ca_o_x, ca_o_y, ca_o_z;

  ca_cb_x=rms;
  
/*
#!/usr/bin/perl
#
#  calc_N_C_CA_CB.pl
#
#  Calculates N, C, CA, CB positions from given distances: 
#    CA-CB, N-CB, N-CA, C-CB, C-CA, C-N
#

$|=1;

####################
# Information derived from the experiments with high 
# resolution structures and
# R.A.Engh and R.Huber, Acta. Cryst., A47 292-300 (1991)
#
# Fixed distances:

$ca_cb=1.532;    $ca_cb=1.524;    1.524
$n_cb=2.443;     $n_cb=2.450;     CB-CA-N ~110.400
$n_ca=1.454;     $n_ca=1.458;     1.458
$c_cb=2.527;     $c_cb=2.490;     C-CA-CB ~110.200
$c_ca=1.538;     $c_ca=1.525;     1.525
$c_n=2.459;      $c_n=2.460;      C-CA-N  ~111.200

####################
# Fixed coordinates:

$ca_x=0.0;
$ca_y=0.0;   
$ca_z=0.0;     

$cb_x=$ca_cb;
$cb_y=0.0;   
$cb_z=0.0;
     
$n_z=0.0;
     
####################

$n_x=($ca_cb*$ca_cb - $n_cb*$n_cb + $n_ca*$n_ca)/(2.0*$ca_cb);
$n_y=-sqrt($n_ca*$n_ca-$n_x*$n_x);

$c_x=($ca_cb*$ca_cb - $c_cb*$c_cb + $c_ca*$c_ca)/(2.0*$ca_cb);
$c_y=($c_ca*$c_ca - $c_n*$c_n + $n_y*$n_y - 2.0*$c_x*$n_x + $n_x*$n_x)*0.5/$n_y;
$c_z=-sqrt($c_ca*$c_ca-$c_x*$c_x-$c_y*$c_y);

printf "N:  %10.7f   %10.7f   %10.7f\n",$n_x,$n_y,$n_z;
printf "C:  %10.7f   %10.7f   %10.7f\n",$c_x,$c_y,$c_z;
printf "CA: %10.7f   %10.7f   %10.7f\n",$ca_x,$ca_y,$ca_z;
printf "CB: %10.7f   %10.7f   %10.7f\n",$cb_x,$cb_y,$cb_z;

exit;
*/

  n=3;
  atoms[0][1]=MAXRES;
  atoms[1][1]=MAXRES;
  atoms[0][2]=MAXRES-2;
  atoms[1][2]=MAXRES-2;
  atoms[0][3]=MAXRES-1;
  atoms[1][3]=MAXRES-1;
  atoms[0][4]=MAXRES-3;
  atoms[1][4]=MAXRES-3;

// CB calculations  
if(in_CB==0) {
  pdb[0].atom[atoms[0][1]].R.x=-0.4918763;                 //  N  
  pdb[0].atom[atoms[0][1]].R.y=-1.3682740;
  pdb[0].atom[atoms[0][1]].R.z=0.0000000;
  pdb[0].atom[atoms[0][2]].R.x=-0.5461035;                 //  C  
  pdb[0].atom[atoms[0][2]].R.y=0.7689804;
  pdb[0].atom[atoms[0][2]].R.z=-1.2148597;
  pdb[0].atom[atoms[0][3]].R.x=0.0000000;                  // CA  
  pdb[0].atom[atoms[0][3]].R.y=0.0000000;
  pdb[0].atom[atoms[0][3]].R.z=0.0000000;

  pdb[0].atom[atoms[0][4]].R.x=1.5320000;                  // CB  
  pdb[0].atom[atoms[0][4]].R.y=0.0000000;
  pdb[0].atom[atoms[0][4]].R.z=0.0000000;

  rms=rmsd(n, atoms, pdb, r, v);

  cb_x=r[0][0]*pdb[0].atom[atoms[0][4]].R.x +
       r[1][0]*pdb[0].atom[atoms[0][4]].R.y +
       r[2][0]*pdb[0].atom[atoms[0][4]].R.z + v[0];
  cb_y=r[0][1]*pdb[0].atom[atoms[0][4]].R.x +
       r[1][1]*pdb[0].atom[atoms[0][4]].R.y +
       r[2][1]*pdb[0].atom[atoms[0][4]].R.z + v[1];
  cb_z=r[0][2]*pdb[0].atom[atoms[0][4]].R.x +
       r[1][2]*pdb[0].atom[atoms[0][4]].R.y +
       r[2][2]*pdb[0].atom[atoms[0][4]].R.z + v[2];
}
else {
  cb_x=pdb[1].atom[MAXRES-4].R.x;
  cb_y=pdb[1].atom[MAXRES-4].R.y;
  cb_z=pdb[1].atom[MAXRES-4].R.z;
}

// Mid polong calculations  
  m_x=(pdb[1].atom[atoms[1][1]].R.x+
       pdb[1].atom[atoms[1][2]].R.x+
       pdb[1].atom[atoms[1][3]].R.x+
       cb_x)/4.0;
  m_y=(pdb[1].atom[atoms[1][1]].R.y+
       pdb[1].atom[atoms[1][2]].R.y+
       pdb[1].atom[atoms[1][3]].R.y+
       cb_y)/4.0;
  m_z=(pdb[1].atom[atoms[1][1]].R.z+
       pdb[1].atom[atoms[1][2]].R.z+
       pdb[1].atom[atoms[1][3]].R.z+
       cb_z)/4.0;

// CA-CB, CA-M, and CA-O vectors calculations
  ca_cb_x=cb_x-pdb[1].atom[atoms[1][3]].R.x;
  ca_cb_y=cb_y-pdb[1].atom[atoms[1][3]].R.y;
  ca_cb_z=cb_z-pdb[1].atom[atoms[1][3]].R.z;

  ca_m_x=m_x-pdb[1].atom[atoms[1][3]].R.x;
  ca_m_y=m_y-pdb[1].atom[atoms[1][3]].R.y;
  ca_m_z=m_z-pdb[1].atom[atoms[1][3]].R.z;

  ca_o_x=pdb[1].atom[atoms[1][4]].R.x-pdb[1].atom[atoms[1][3]].R.x;
  ca_o_y=pdb[1].atom[atoms[1][4]].R.y-pdb[1].atom[atoms[1][3]].R.y;
  ca_o_z=pdb[1].atom[atoms[1][4]].R.z-pdb[1].atom[atoms[1][3]].R.z;

// BMO hypothetical polong calculations  
  pdb[0].atom[atoms[1][3]].R.x=cb_v*ca_cb_x + m_v*ca_m_x + o_v*ca_o_x + 
                               pdb[1].atom[atoms[1][3]].R.x;
  pdb[0].atom[atoms[1][3]].R.y=cb_v*ca_cb_y + m_v*ca_m_y + o_v*ca_o_y + 
                               pdb[1].atom[atoms[1][3]].R.y;
  pdb[0].atom[atoms[1][3]].R.z=cb_v*ca_cb_z + m_v*ca_m_z + o_v*ca_o_z + 
                               pdb[1].atom[atoms[1][3]].R.z;

  return;
}
