#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tools.h"

/*-------------------------------------------------------------
/
/   rmsd_isp - ISP procedure for RMSD calculation with a given
/              DISTANCE cutoff
/
/------------------------------------------------------------*/
float rmsd_isp(long n_atoms, long atoms[2][MAXRES+1], pdb_struct pdb[2],
              float r[3][3], float v[3], part_struct *part)
{
  typedef struct {
    long   atoms[2][MAXRES+1];
  } iter_struct;
  /* Large array moved off stack to avoid stack overflow on systems
     with small default thread stack sizes. */
  static iter_struct  calphas[MAXITR+1];
  long          i, k, k1, k2, n, nr, ok, out, restart;
  float        rms;

  for(i=0;i<=n_atoms;i++) {
    calphas[0].atoms[0][i]=atoms[0][i];
    calphas[0].atoms[1][i]=atoms[1][i];
  }
  n=n_atoms;

  k=part->isp_iter_start;
  nr=0;
  out=0;
  restart=0;
/*
  printf("\n");
*/
  while(k<MAXITR && out!=1) {
    part->isp_iter=k;
    n=calphas[k].atoms[0][0];
    if(n<3) {
      out=1;
      if(restart<=10) { 
        nr=21;
        for(i=0;i<=20;i++) {
          if(part->opt_sup1[i][0]>=3 && out==1) {
            nr=i+restart;
            out=0;
          }
        }
        if(out==0) {
          if(nr>20) nr=20;
          n=part->opt_sup1[nr][0];
          for(i=0;i<=n;i++) {
            calphas[k].atoms[0][i]=part->opt_sup1[nr][i];
            calphas[k].atoms[1][i]=part->opt_sup2[nr][i];
          }
        }
      }
      restart++;
    }
    if(out==0) {
      n=calphas[k].atoms[0][0];
      rms=rmsd(n, calphas[k].atoms, pdb, r, v);
      if(part->isp_iter==0) part->rmsd=rms;
/*
      printf("TEST (loc) RMSD = %7.3f A  Number of atoms = %4ld  Cutoff  = %7.3f A  Nr = %ld \n",
                rms,n,part->isp_dist,nr);
*/
      n=calphas[0].atoms[0][0];
      rms=rmsd_dist(n, calphas[0].atoms, pdb, r, v, part);
      n=0;
      for(i=1;i<=calphas[0].atoms[0][0];i++) {
        if(calphas[0].atoms[0][i]>0) {
          n++;
          calphas[k+1].atoms[0][n]=atoms[0][i];
          calphas[k+1].atoms[1][n]=atoms[1][i];
        }
        else {
          calphas[0].atoms[0][i]=atoms[0][i];
          calphas[0].atoms[1][i]=atoms[1][i];
        }
      }
      calphas[k+1].atoms[0][0]=n;
      calphas[k+1].atoms[1][0]=n;
      if(n<3) {
        ok=0;
      }
      else {
        ok=1;
      }
      k1=1;
      while(k1<=k && ok==1) {
        ok=0;
        if(n==calphas[k1].atoms[0][0]) {
          for(k2=0;k2<=n;k2++) {
            if(calphas[k1].atoms[0][k2]!=calphas[k+1].atoms[0][k2] ||
               calphas[k1].atoms[1][k2]!=calphas[k+1].atoms[1][k2]) {
              ok=1;
              k2=n;
            }
          }
        }
        else {
          ok=1;
        }
        if(ok==1) {
          k1++;
        }
      }
      if(ok==0 && part->isp_iter+restart>20) { 
        out=1;
      }
      else {
        k++;
      }
    }
  }
  n=part->opt_sup1[0][0];

  if(n<3) {
    if(part->error[0]==0) {
      part->error[0]=1;
      printf("\n# WARNING! The change of the distance cutoff DIST may give you better result.\n");
    }
    return 999.99;
  }
  for(i=0;i<=n;i++) {
    calphas[0].atoms[0][i]=part->opt_sup1[0][i];
    calphas[0].atoms[1][i]=part->opt_sup2[0][i];
  }
  part->opt_sup_rms[0]=rmsd(n, calphas[0].atoms, pdb, r, v);

  for(i=0;i<3;i++) {
    r[0][i]=part->opt_r[0][0][i];
    r[1][i]=part->opt_r[0][1][i];
    r[2][i]=part->opt_r[0][2][i];
    v[i]=part->opt_v[0][i];
  }

  if(part->opt_all_rms[0]*2.0<part->rmsd && part->opt_all_rms[0]>0.02) {
    printf("\n# ERROR! Check ATOM records and the distances between the atoms.");
    printf("\n#        The accuracy of the calculations is insufficient!\n");
    printf("\nSUMMARY(LGA)    0    0    0.0      0  999.99     0.00      0.000     0.000\n");

    exit(0);
  }
/*
      printf("TEST (out) RMSD = %7.3f A  Number of atoms = %4ld  Cutoff  = %7.3f A  Nr = %ld \n",
                part->opt_sup_rms[0],n,part->isp_dist,nr);
*/
  return part->opt_all_rms[0];
}

/*--------------------------------------------------------------------
/
/   rmsd_dist - calculates RMSD and number of residues that can fit
/               under specified DISTANCE cutoff
/
/--------------------------------------------------------------------*/
float rmsd_dist(long n_atoms, long atoms[2][MAXRES+1], pdb_struct pdb[2],
                float r[3][3], float v[3], part_struct *part)
{
  long    i,j,k,n,opt_k[21],i1;
  /* Large arrays moved off stack to avoid stack overflow. */
  static long  opt1[21][MAXRES], opt2[21][MAXRES];
  double frms, rms[21], trms, s, t, xv, yv, zv, x, y, z;

  k=0;
  trms=0.0;
  for(i=0;i<=20;i++) {
    opt_k[i]=0;
    rms[i]=999.999;
  }

  if(part->isp_iter<11) {
    part->isp_dist=part->dist_cutoff*part->isp_cutoff[part->isp_iter];
    while(part->isp_dist>100.0) part->isp_dist=part->isp_dist/2.0;
  }                                            
  else {
    part->isp_dist=part->dist_cutoff;
  } 
  for(i=0;i<n_atoms;i++) {
    i1=i+1;
    xv=pdb[1].atom[atoms[1][i1]].R.x - v[0];
    yv=pdb[1].atom[atoms[1][i1]].R.y - v[1];
    zv=pdb[1].atom[atoms[1][i1]].R.z - v[2];
    x=pdb[0].atom[atoms[0][i1]].R.x -
      (r[0][0]*xv+r[0][1]*yv+r[0][2]*zv);
    y=pdb[0].atom[atoms[0][i1]].R.y -
      (r[1][0]*xv+r[1][1]*yv+r[1][2]*zv);
    z=pdb[0].atom[atoms[0][i1]].R.z -
      (r[2][0]*xv+r[2][1]*yv+r[2][2]*zv);
    t=x*x + y*y + z*z;
    trms=trms + t;
    s=sqrt(t);
    for(k=0;k<=20;k++) {
      if(s<=part->gdt_cutoff[k]) {
        opt_k[k]++;
        j=opt_k[k];
        rms[k]=rms[k] + t;
        opt1[k][j]=atoms[0][i1];
        opt2[k][j]=atoms[1][i1];
      }
    }
    if(s>part->isp_dist && part->isp==1) {
      atoms[0][i1]=0;
      atoms[1][i1]=0;
    }
  }

  if(n_atoms>0)
    frms=sqrt(trms/n_atoms);
  else
    frms=999.999;

  for(k=0;k<=20;k++) {
    j=opt_k[k];
    if(j>0) {
      rms[k]=sqrt(rms[k]/j);
      for(i=1;i<=j;i++) {
        n=pdb[1].atom[opt2[k][i]].res_n_local;
        if(part->opt_set_nb[k][n]<j) {
          part->opt_set_nb[k][n]=j;
        }
      }
      if(j>part->opt_sup1[k][0]) {
        for(i=1;i<=j;i++) {
          part->opt_sup1[k][i]=opt1[k][i];
          part->opt_sup2[k][i]=opt2[k][i];
        }
        part->opt_sup1[k][0]=j;
        part->opt_sup2[k][0]=j;
        part->opt_sup_rms[k]=rms[k];
        part->opt_all_rms[k]=frms;
        for(n=0;n<3;n++) {
          part->opt_r[k][0][n]=r[0][n];
          part->opt_r[k][1][n]=r[1][n];
          part->opt_r[k][2][n]=r[2][n];
          part->opt_v[k][n]=v[n];
        }
      }
      else if(j==part->opt_sup1[k][0] && rms[k]<=part->opt_sup_rms[k]) {
        for(i=1;i<=j;i++) {
          part->opt_sup1[k][i]=opt1[k][i];
          part->opt_sup2[k][i]=opt2[k][i];
        }
        part->opt_sup_rms[k]=rms[k];
        part->opt_all_rms[k]=frms;
        for(n=0;n<3;n++) {
          part->opt_r[k][0][n]=r[0][n];
          part->opt_r[k][1][n]=r[1][n];
          part->opt_r[k][2][n]=r[2][n];
          part->opt_v[k][n]=v[n];
        }
      }
    }
  }

  return (float)frms;
}

/*---------------------------------------------------------------------
/
/   lcs_run - for given rms cutoff calculates the max lenght 
/             of the sequence that can fit under RMS cutoff
/
/--------------------------------------------------------------------*/
float lcs_run(long n_atoms, long atoms[2][MAXRES+1], pdb_struct pdb[2],
            float r[3][3], float v[3], part_struct *part)
{
  long    i, j, k1, k2, k_opt, n, n1, n2, out;
  long    j1, j2, mn, max_n, next;
  long    m_k1_opt[500], m_k2_opt[500];
  float  m_rms_local[500], m_rms_all[500];
  float  rms_local, rms_all, lcs;
  /* Large array moved off stack to avoid stack overflow. */
  static long  calphas[2][MAXRES+1];
  char   si1[10], si2[10], sj1[10], sj2[10];

  if(part->lcs_gdt_print==1) 
    printf("\nLCS - RMSD CUTOFF %6.2f      length       segment         l_RMS    g_RMS",
            part->rms_lcs_cutoff);

  for(i=0;i<=n_atoms;i++) {
    part->opt_lcs[i]=0;
  }

  j1=0;
  j2=0;
  strcpy(si1," ");
  strcpy(si2," ");
  strcpy(sj1," ");
  strcpy(sj2," ");

  max_n=0;
  k1=1;
  k2=3;
  k_opt=3;
  m_k1_opt[max_n]=0;
  m_k2_opt[max_n]=0;
  n=0;
  n1=0;
  n2=0;
  out=0;
  while(k2<=n_atoms) {
    n1=k1;
    n2=k2;
    next=0;
    while(n2<=n_atoms) {
      n=n2-n1+1;    
      for(i=n1;i<=n2;i++) {
        j=i-n1+1;
        calphas[0][j]=atoms[0][i];
        calphas[1][j]=atoms[1][i];
      }
      calphas[0][0]=n;
      calphas[1][0]=n;

      rms_local=rmsd(n, calphas, pdb, r, v);
      rms_all=rmsd_dist(n_atoms, atoms, pdb, r, v, part);

      if(rms_local<part->rms_lcs_cutoff) {
        if(n>k_opt) {
          max_n=0;
          k_opt=n;
          out=1;
        }
        if(n>=k_opt) {
          if(m_k1_opt[max_n]!=n1 || m_k2_opt[max_n]!=n2) {
            m_k1_opt[max_n]=n1;
            m_k2_opt[max_n]=n2;
            m_rms_local[max_n]=rms_local;
            m_rms_all[max_n]=rms_all;
            max_n++;
          }
        }
        for(i=n1;i<=n2;i++) {
          if(part->opt_lcs[i]<n) part->opt_lcs[i]=n;
        }
        k2=n2;
        next=1;
      }
      else {
        next=0;
      }
      n2++;
    }
    if(next==0) {
      k1++;
    }
    else {
      k2=n2;
    }
    if(k1>k2-2) k2=k1+2;
  }

  if(out==0) {
    if(part->error[1]==0) {
      part->error[1]=1;
      printf("\n# WARNING! The change of the parameter DIST cutoff may give you better result.\n");
    }
/*
    return 0.0;
*/
  }
  mn=0;
  while(mn<max_n) {
    k1=m_k1_opt[mn];
    k2=m_k2_opt[mn];
    j1=atoms[1][k1];
    j2=atoms[1][k2];
    strcpy(sj1,pdb[1].atom[j1].res_i);
    strcpy(sj2,pdb[1].atom[j2].res_i);
    if(part->lcs_gdt_print==1) 
      printf("\n  LONGEST_CONTINUOUS_SEGMENT:  %4ld   %7s - %-7s   %6.2f   %6.2f",
                k_opt,sj1,sj2,m_rms_local[mn],m_rms_all[mn]);
    mn++;
    for(i=k1;i<=k2;i++) {
      if(part->opt_lcs[i]<k_opt) part->opt_lcs[i]=k_opt;
    }
  }

  lcs=0.0;
  if(n_atoms>0) {
    for(i=1;i<=n_atoms;i++) {
      lcs=lcs+part->opt_lcs[i];
    }
    lcs=100.0*lcs/(n_atoms*part->t_n_aa);
    if(part->lcs_gdt_print==1) 
      printf("\n  LCS_AVERAGE:    %6.2f\n",lcs);
  }

  return lcs;
}

/*-----------------------------------------------------------
/
/   gdt_run - gdt calculations with window 3
/
/------------------------------------------------------------*/
void gdt_run(long n_atoms, long atoms[2][MAXRES+1], pdb_struct pdb[2],
            float r[3][3], float v[3], part_struct *part)
{
  long    i, j, k1, k2, k, n, print_test=0;
  float  rms_local=0.0, rms_all=0.0;
  /* Large array moved off stack to avoid stack overflow. */
  static long  calphas[2][MAXRES+1];

  k1=1;
  k2=3;
  while(k2<=n_atoms) {
    k=k2-k1+1;    
    for(i=k1;i<=k2;i++) {
      j=i-k1+1;
      calphas[0][j]=atoms[0][i];
      calphas[1][j]=atoms[1][i];
    }
    calphas[0][0]=k;
    calphas[1][0]=k;

    rms_local=rmsd(k, calphas, pdb, r, v);
    for(i=0;i<=n_atoms;i++) {
      calphas[0][i]=atoms[0][i];
      calphas[1][i]=atoms[1][i];
    }
    part->isp=1;
    rms_all=rmsd_dist(n_atoms, calphas, pdb, r, v, part);

    n=0;
    for(i=1;i<=n_atoms;i++) {
      if(calphas[0][i]>0) {
        n++;
        calphas[0][n]=atoms[0][i];
        calphas[1][n]=atoms[1][i];
      }
    }
    calphas[0][0]=n;
    calphas[1][0]=n;
    if(n>=3) {
      rms_local=rmsd(n, calphas, pdb, r, v);
      part->isp=0;
      rms_all=rmsd_dist(n_atoms, atoms, pdb, r, v, part);
      k1++;
    }
    else {
      k2++;
      if(k2>n_atoms) {
        k1++;
        k2=k1+2;
      }
    }
    if(k1>k2-2) k2=k1+2;
  }
  if(print_test==1) {
    printf("\nTEST_GDT_RMS: %7.3f  %7.3f \n",rms_local,rms_all); 
  }

  return;
}

/*-------------------------------------------------------------
/
/   lcs_gdt_analysis - LCS and GDT analysis
/
/------------------------------------------------------------*/
void lcs_gdt_analysis(long n_ca_all, long atoms[2][MAXRES+1], pdb_struct pdb[2],
                      float r[3][3], float v[3], part_struct *part)
{
  long    i, j, j1, j2, n;
  /* Large arrays moved off stack to avoid stack overflow. */
  static long  opt_lcs125[3][MAXRES+1], gdt_chars[21][MAXRES+1], l_atoms[2][MAXRES+1];
  float  lcs[3], rms=0.0, step, q_norm, r_norm, lcs_gdt_ts;
  char   ch1, ch2;

  q_norm=rms;
  r_norm=(float)n_ca_all/part->t_n_aa;

/*  Clean part structure
--------------------------*/
  for(i=0;i<=20;i++) {
    part->opt_sup1[i][0]=0;
    part->opt_sup2[i][0]=0;
    for(j=0;j<=n_ca_all;j++) {
      part->opt_set_nb[i][j]=0;
    }
  }
  for(j=0;j<=n_ca_all;j++) {
    part->opt_lcs[j]=0;
  }
  lcs[0]=0;
  lcs[1]=0;
  lcs[2]=0;

/*  LCS calculations: 5, 2, 1
-------------------------------*/

  if(part->lcs_gdt_print>0) {
    part->isp=0;
    part->rms_lcs_cutoff=5.0;
    lcs[2]=lcs_run(n_ca_all, atoms, pdb, r, v, part);
    for(i=0;i<=n_ca_all;i++) {
      opt_lcs125[2][i]=part->opt_lcs[i];
    }
    part->rms_lcs_cutoff=2.0;
    lcs[1]=lcs_run(n_ca_all, atoms, pdb, r, v, part);
    for(i=0;i<=n_ca_all;i++) {
      opt_lcs125[1][i]=part->opt_lcs[i];
    }
    part->rms_lcs_cutoff=1.0;
    lcs[0]=lcs_run(n_ca_all, atoms, pdb, r, v, part);
    for(i=0;i<=n_ca_all;i++) {
      opt_lcs125[0][i]=part->opt_lcs[i];
    }
  }
  else {
    part->rms_lcs_cutoff=1.4;
    lcs[1]=lcs_run(n_ca_all, atoms, pdb, r, v, part);
    for(i=0;i<=n_ca_all;i++) {
      opt_lcs125[1][i]=part->opt_lcs[i];
    }
  }

/*  GDT calculations: window 3, segments under distance cutoffs 0.1 - 5.6
---------------------------------------------------------------------------*/

  if(part->lcs_gdt_print>0) {
    n=part->accuracy_gdt_n;
    step=part->accuracy_gdt_step;
  }
  else {
    n=part->accuracy_lga_n;
    step=part->accuracy_lga_step;
  }

  part->isp_iter_start=11;
  for(i=0;i<n;i++) {
    part->dist_cutoff=step*(n-i);
    gdt_run(n_ca_all, atoms, pdb, r, v, part);
  }
  part->isp_iter_start=0;
  part->isp=1;
  for(i=1;i<=n;i++) {
    for(j=0;j<=n_ca_all;j++) {
      l_atoms[0][j]=atoms[0][j];
      l_atoms[1][j]=atoms[1][j];
    }
    part->dist_cutoff=step*(n-i);
    rms=rmsd_isp(n_ca_all, l_atoms, pdb, r, v, part);
  }
  
/*  GDT calculations: standard distance cutoffs  
------------------------------------------------*/

  part->isp_iter_start=0;
  part->isp=1;
  for(i=0;i<=20;i++) {
    part->dist_cutoff=part->gdt_cutoff[i];
    for(j=0;j<=n_ca_all;j++) {
      l_atoms[0][j]=atoms[0][j];
      l_atoms[1][j]=atoms[1][j];
    }
    rms=rmsd_isp(n_ca_all, l_atoms, pdb, r, v, part);
  }
  part->dist_cutoff=part->gdt_cutoff[0];

/*  GDT calculations: the best 
-------------------------------*/

  part->isp=0;
  for(j=20;j>=0;j--) {
    n=part->opt_sup1[j][0];
    if(n>=3) {
      for(i=0;i<=n;i++) {
        l_atoms[0][i]=part->opt_sup1[j][i];
        l_atoms[1][i]=part->opt_sup2[j][i];
      }
      part->opt_sup_rms[j]=rmsd(n, l_atoms, pdb, r, v);
      part->rmsd_isp=rmsd_dist(n_ca_all, atoms, pdb, r, v, part);
      part->opt_all_rms[j]=part->rmsd_isp;
    }
  }
  for(j=20;j>0;j--) {
    if(part->gdt_cutoff[0]>=part->gdt_cutoff[j]) {
      if(part->opt_sup1[j][0]>part->opt_sup1[0][0] ||
         (part->opt_sup1[j][0]==part->opt_sup1[0][0] &&
          part->opt_sup_rms[j]<part->opt_sup_rms[0])) {
        for(i=0;i<=part->opt_sup1[j][0];i++) {
          part->opt_sup1[0][i]=part->opt_sup1[j][i];
          part->opt_sup2[0][i]=part->opt_sup2[j][i];
        }
        part->opt_sup_rms[0]=part->opt_sup_rms[j];
        for(n=0;n<3;n++) {
          part->opt_r[0][0][n]=part->opt_r[j][0][n];
          part->opt_r[0][1][n]=part->opt_r[j][1][n];
          part->opt_r[0][2][n]=part->opt_r[j][2][n];
          part->opt_v[0][n]=part->opt_v[j][n];
        }
      }
    }
  }
  if(part->gdt==0) { // standard superposition on identified N residues
    n=part->opt_sup1[0][0];
    if(n>=3) {
      for(i=0;i<=n;i++) {
        l_atoms[0][i]=part->opt_sup1[0][i];
        l_atoms[1][i]=part->opt_sup2[0][i];
      }
      part->rmsd_isp=rmsd(n, l_atoms, pdb, r, v); // calculates rotation r and shift v
    }
    else {
      part->gdt=1;
    }
  }  
  if(part->gdt==1) {  // superposition that fits maximum number of residues
    for(i=0;i<3;i++) {
      r[0][i]=part->opt_r[0][0][i];
      r[1][i]=part->opt_r[0][1][i];
      r[2][i]=part->opt_r[0][2][i];
      v[i]=part->opt_v[0][i];
    }
  }
  part->rmsd_isp=rmsd_dist(n_ca_all, atoms, pdb, r, v, part);

  part->gdt_ts=0.0;
  for(i=2;i<=19;i++) {
    j1=19-i+1;
    part->gdt_ts=part->gdt_ts + j1*part->opt_sup2[i][0];
  }
  part->gdt_ts=part->gdt_ts*r_norm/(1.71*n_ca_all);
  part->lcs_ts=lcs[1];

// NOTE: GDT results (LGA_S) may differ when calculated with "-3" and "-4" modes
//       EVEN IF PERFORMED ON THE SAME SET OF RESIDUES.
//       "-3" and "-4" use different procedures to generate sets of residue-pairs
//       to calculate "optimal" superpositions (maximum number of residues under 
//       distance cutoffs: 0.5, 1.0, 1.5, ... 10.0) 
/*
  printf("\nTEST_LCS_GDT: %7.3f  %7.3f \n",part->lcs_ts,part->gdt_ts); 
  for(j=20;j>0;j--) {
    printf(" %4ld",part->opt_sup2[j][0]);
  }
  printf("\n"); 
  
  printf("\nUnitary ROTATION matrix and the shift VECTOR superimpose MOLECULES  (1=>2)");
  printf("\n  X_new = %14.10f * X  + %14.10f * Y  + %14.10f * Z  + %14.10f",
           r[0][0],r[1][0],r[2][0],v[0]);
  printf("\n  Y_new = %14.10f * X  + %14.10f * Y  + %14.10f * Z  + %14.10f",
           r[0][1],r[1][1],r[2][1],v[1]);
  printf("\n  Z_new = %14.10f * X  + %14.10f * Y  + %14.10f * Z  + %14.10f \n",
           r[0][2],r[1][2],r[2][2],v[2]);
*/

  if(part->lcs_gdt_print!=1) {
    return;
  }

/* Presentation of the results
-------------------------------*/

/*
    part->lcs_ts=(5*lcs[0]+4*lcs[1]+lcs[2])/10.0;
*/
  lcs_gdt_ts=part->lga_w*part->gdt_ts+(1.0-part->lga_w)*part->lcs_ts;
  part->lcs_ts=(lcs[0]+lcs[1]+lcs[2])/3.0;
  part->gdt_ts=(part->opt_sup2[2][0]+part->opt_sup2[4][0]+
        part->opt_sup2[8][0]+part->opt_sup2[16][0])*25.0/n_ca_all;

  if(n_ca_all>0) {
    for(n=1;n<=20;n++) {
      for(i=1;i<=n_ca_all;i++) {
        j=atoms[1][i];
        j=pdb[1].atom[j].res_n_local;
        gdt_chars[n-1][i]=part->opt_set_nb[n][j];
      }
    }
  }
  else {
    for(n=0;n<=20;n++) {
      for(i=1;i<=n_ca_all;i++) {
        gdt_chars[n][i]=0;
      }
    }
  }

  n=n_ca_all;

  printf("\nLCS_GDT    MOLECULE-1    MOLECULE-2     LCS_DETAILS     GDT_DETAILS                                                    TOTAL NUMBER OF RESIDUE PAIRS: %4ld",n_ca_all);
  printf("\nLCS_GDT     RESIDUE       RESIDUE       SEGMENT_SIZE    GLOBAL DISTANCE TEST COLUMNS: number of residues under the threshold assigned to each residue pair");
  printf("\nLCS_GDT   NAME NUMBER   NAME NUMBER    1.0  2.0  5.0    0.5  1.0  1.5  2.0  2.5  3.0  3.5  4.0  4.5  5.0  5.5  6.0  6.5  7.0  7.5  8.0  8.5  9.0  9.5 10.0");

  for(i=1;i<=n;i++) {
    j1=atoms[0][i];
    C_321(pdb[0].atom[j1].res_name,&ch1);
    j2=atoms[1][i];
    C_321(pdb[1].atom[j2].res_name,&ch2);
    printf("\nLCS_GDT     %1c %7s     %1c %7s   %4ld %4ld %4ld   %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld ",
            ch1,pdb[0].atom[j1].res_i,ch2,pdb[1].atom[j2].res_i,
            opt_lcs125[0][i],opt_lcs125[1][i],opt_lcs125[2][i],
            gdt_chars[0][i],gdt_chars[1][i],gdt_chars[2][i],gdt_chars[3][i],
            gdt_chars[4][i],gdt_chars[5][i],gdt_chars[6][i],gdt_chars[7][i],
            gdt_chars[8][i],gdt_chars[9][i],
            gdt_chars[10][i],gdt_chars[11][i],gdt_chars[12][i],gdt_chars[13][i],
            gdt_chars[14][i],gdt_chars[15][i],gdt_chars[16][i],gdt_chars[17][i],
            gdt_chars[18][i],gdt_chars[19][i]);
  }

  printf("\nLCS_AVERAGE  LCS_A: %6.2f  ( %6.2f  %6.2f  %6.2f )\n",
          part->lcs_ts,lcs[0],lcs[1],lcs[2]);
  printf("\nGLOBAL_DISTANCE_TEST (summary information about detected largest sets of residues (represented by selected AToms) that can fit under specified thresholds)");
  printf("\nGDT DIST_CUTOFF");
  printf("%6.2f",part->gdt_cutoff[1]);
  for(i=2;i<=20;i++) {
    printf(" %6.2f",part->gdt_cutoff[i]);
  }

  printf("\nGDT NUMBER_AT ");
  for(i=1;i<=20;i++) {
    printf("%6ld ",part->opt_sup2[i][0]);
  }

  printf("\nGDT PERCENT_AT");
  for(i=1;i<=20;i++) {
    printf(" %6.2f",part->opt_sup2[i][0]*100.0*r_norm/n);
  }

  printf("\nGDT RMS_LOCAL ");
  for(i=1;i<=20;i++) {
    printf(" %6.2f",part->opt_sup_rms[i]);
  }

  printf("\nGDT RMS_ALL_AT");
  for(i=1;i<=20;i++) {
    printf(" %6.2f",part->opt_all_rms[i]);
  }
  
  part->gdt_ts=part->gdt_ts*r_norm;

  n=part->opt_sup1[0][0];
  q_norm=part->opt_sup_rms[0];
  if(q_norm<0.0) {
    q_norm=0.0;
  }
  else {
    q_norm=n*0.1/(q_norm+0.1);
  }
  if(part->opt_sup_rms[0]>999 || part->opt_sup_rms[0]<0) {
    lcs_gdt_ts=0.0;
  }
  part->summary.rms_local=part->opt_sup_rms[0];
  part->summary.lcs_gdt_ts=lcs_gdt_ts;

  printf("\n");
  sprintf(part->summary.s[0].line,"\n#%-8s      N1   N2   DIST      N    RMSD    GDT_TS    LGA_S3     LGA_Q ",part->atoms);
  sprintf(part->summary.s[1].line,"\nSUMMARY(GDT) %4ld %4ld   %4.1f   %4ld  %6.2f   %7.3f   %7.3f  %8.3f",
          part->m_n_aa,part->t_n_aa,part->dist_cutoff,n,part->opt_sup_rms[0],part->gdt_ts,lcs_gdt_ts,q_norm);
  part->summary.n=2;

  return;
}
