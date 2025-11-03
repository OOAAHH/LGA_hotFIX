#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tools.h"

/*--------------------------------------------------------------------
/
/   align_search - search for the alignment under given superposition
/
/--------------------------------------------------------------------*/
long align_search(long anchor[2], long atoms[2][MAXRES+1], pdb_struct pdb[2],
                 float r[3][3], float v[3], part_struct *part, part_ear *ear)
{
  long    i,j,i0,j0,ok,n,na;
  lga_workspace *ws = part->ws;
  if (!ws) { /* defensive: ensure workspace exists */
    if (lga_ws_init(&part->ws) != 0) {
      fprintf(stderr, "\n# ERROR! Failed to allocate LGA workspace (align_search)\n");
      exit(1);
    }
    ws = part->ws;
  }
  long  *tmp = ws ? ws->tmp : NULL;
  float (*xa)[3] = ws ? (float (*)[3])ws->xa3 : NULL;
  float (*xb)[3] = ws ? (float (*)[3])ws->xb3 : NULL;
  float dist2, t, t10, t01, x, y, z;

  dist2=part->dist_cutoff*part->dist_cutoff;
  
  for(i=0;i<atoms[0][0];i++) {
    j=i+1;
    xa[i][0]=pdb[0].atom[atoms[0][j]].R.x;
    xa[i][1]=pdb[0].atom[atoms[0][j]].R.y;
    xa[i][2]=pdb[0].atom[atoms[0][j]].R.z;
  }

  for(i=0;i<atoms[1][0];i++) {
    j=i+1;
    xb[i][0]=r[0][0]*(pdb[1].atom[atoms[1][j]].R.x - v[0])+
             r[0][1]*(pdb[1].atom[atoms[1][j]].R.y - v[1])+
             r[0][2]*(pdb[1].atom[atoms[1][j]].R.z - v[2]);
    xb[i][1]=r[1][0]*(pdb[1].atom[atoms[1][j]].R.x - v[0])+
             r[1][1]*(pdb[1].atom[atoms[1][j]].R.y - v[1])+
             r[1][2]*(pdb[1].atom[atoms[1][j]].R.z - v[2]);
    xb[i][2]=r[2][0]*(pdb[1].atom[atoms[1][j]].R.x - v[0])+
             r[2][1]*(pdb[1].atom[atoms[1][j]].R.y - v[1])+
             r[2][2]*(pdb[1].atom[atoms[1][j]].R.z - v[2]);
  }                             

  n=0;
  na=0;
  i0=anchor[0]-1;
  j0=anchor[1]-1;
  if(anchor[0]<1 || anchor[1]<1 || 
     anchor[0]>pdb[0].n_aa || anchor[1]>pdb[1].n_aa) {
    printf("\n ERROR! %ld %ld %ld %ld",anchor[0],anchor[1],atoms[0][0],atoms[1][0]);
    return n;
  }

  ok=0;
  while((i0>0 || j0>0) && ok<2) {
    n=1;
    ok=0;
    i=i0;
    j=j0;
    while(ok==0 && n<=MAXGAP2) {
      i=i0-part->pairs[0][n];
      j=j0-part->pairs[1][n];
      if(i<0 && j<0) {
        ok=2;
      }
      else if(i>=0 && j>=0) {
        x=xa[i][0] - xb[j][0];
        y=xa[i][1] - xb[j][1];
        z=xa[i][2] - xb[j][2];
        t=x*x + y*y + z*z;
        t10=t;
        t01=t;
        if(t<dist2) {
          na++;

          if(i>0) {
            x=xa[i-1][0] - xb[j][0];
            y=xa[i-1][1] - xb[j][1];
            z=xa[i-1][2] - xb[j][2];
            t10=x*x + y*y + z*z;
          }
          if(j>0) {
            x=xa[i][0] - xb[j-1][0];
            y=xa[i][1] - xb[j-1][1];
            z=xa[i][2] - xb[j-1][2];
            t01=x*x + y*y + z*z;
          }
          if(t>t10 || t>t01) {
            if(t10>t01) {
              t=t01;
              j--;
            }
            else {
              t=t10;
              i--;
            }
          }

          ear->equiv[j+1]=i+1;
          ear->rmsd[j+1]=t;
          tmp[na]=j+1;
          ok=1;
          i0=i;
          j0=j;
        }
      }
      n++;
    }
    if(n>MAXGAP2) {
      i0=i;
      j0=j;
    }
  }
  if(na>0) for(i=1;i<=na;i++) ear->align[i]=tmp[na-i+1];
  i0=anchor[0]-1;
  j0=anchor[1]-1;
  i=i0;
  j=j0;
  x=xa[i][0] - xb[j][0];
  y=xa[i][1] - xb[j][1];
  z=xa[i][2] - xb[j][2];
  t=x*x + y*y + z*z;
  if(t<dist2) {
    na++;
    ear->equiv[j+1]=i+1;
    ear->align[na]=j+1;
    ear->rmsd[j+1]=t;
  }
  ok=0;
  while((i0<atoms[0][0]-1 || j0<atoms[1][0]-1) && ok<2) {
    n=1;
    ok=0;
    i=i0;
    j=j0;
    while(ok==0 && n<=MAXGAP2) {
      i=i0+part->pairs[0][n];
      j=j0+part->pairs[1][n];
      if(i>=atoms[0][0] && j>=atoms[1][0]) {
        ok=2;
      }
      else if(i<atoms[0][0] && j<atoms[1][0]) {
        x=xa[i][0] - xb[j][0];
        y=xa[i][1] - xb[j][1];
        z=xa[i][2] - xb[j][2];
        t=x*x + y*y + z*z;
        t10=t;
        t01=t;
        if(t<dist2) {
          na++;

          if(i<atoms[0][0]-1) {
            x=xa[i+1][0] - xb[j][0];
            y=xa[i+1][1] - xb[j][1];
            z=xa[i+1][2] - xb[j][2];
            t10=x*x + y*y + z*z;
          }
          if(j<atoms[1][0]-1) {
            x=xa[i][0] - xb[j+1][0];
            y=xa[i][1] - xb[j+1][1];
            z=xa[i][2] - xb[j+1][2];
            t01=x*x + y*y + z*z;
          }
          if(t>t10 || t>t01) {
            if(t10>t01) {
              t=t01;
              j++;
            }
            else {
              t=t10;
              i++;
            }
          }

          ear->equiv[j+1]=i+1;
          ear->align[na]=j+1;
          ear->rmsd[j+1]=t;
          ok=1;
          i0=i;
          j0=j;
        }
      }
      n++;
    }
    if(n>MAXGAP2) {
      i0=i;
      j0=j;
    }
  }
  ear->align[0]=na;

  return na;
}

/*-----------------------------------------------------------
/
/   alignment_run - alignment searching with window 'nw'
/
/------------------------------------------------------------*/
long alignment_run(long n_atoms, long atoms[2][MAXRES+1],
                  long all_atoms[2][MAXRES+1], pdb_struct pdb[2],
                  float r[3][3], float v[3], part_struct *part)
{
  long i, j, j1, j2, k1, k2, n, nw, nw2, n_max, n_step, min_pdb;
  long n1=0, n2=0;
  float rms_local;
  lga_workspace *ws = part->ws;
  if (!ws) { /* defensive */
    if (lga_ws_init(&part->ws) != 0) {
      fprintf(stderr, "\n# ERROR! Failed to allocate LGA workspace (alignment_run)\n");
      exit(1);
    }
    ws = part->ws;
  }
  long (*calphas)[MAXRES+1] = ws ? (long (*)[MAXRES+1])ws->calphas2 : NULL;
  long anchor[2];
  part_ear ear1, ear2;

  n=n1;
  n_max=n2;
  
  if(part->accuracy_opt>1) {
    nw=4;              //  4
  }
  else {
    nw=6;              //  5  // 6
  }
  n_step=3;            //  2  // 3
  nw2=(long)(nw/2);

  if(pdb[0].n_aa>pdb[1].n_aa) min_pdb=pdb[1].n_aa;
  else min_pdb=pdb[0].n_aa;

  calphas[0][0]=nw;
  calphas[1][0]=nw;
  for(k2=nw;k2<=n_atoms;k2=k2+n_step) {
    for(i=1;i<=nw;i++) {
      j=k2-nw+i;
      calphas[0][i]=atoms[0][j];
      calphas[1][i]=atoms[1][j];
    }
    rms_local=rmsd(nw, calphas, pdb, r, v);
    j1=atoms[0][k2-nw2];
    j2=atoms[1][k2-nw2];
    anchor[0]=pdb[0].atom[j1].res_n_local;
    anchor[1]=pdb[1].atom[j2].res_n_local;
    n=align_search(anchor, all_atoms, pdb, r, v, part, &part->ear[10]);
/*
      printf("\n TEST 0 rmsd = %6.2f  %ld  %ld  %ld  %ld ",rms_local,n,k2,anchor[0],anchor[1]);
*/      
    if(n>2) {
      for(i=1;i<=n;i++) {
        j2=part->ear[10].align[i];
        j1=part->ear[10].equiv[j2];
        calphas[0][i]=all_atoms[0][j1];
        calphas[1][i]=all_atoms[1][j2];
      }
      calphas[0][0]=n;
      calphas[1][0]=n;
      rms_local=rmsd(n, calphas, pdb, r, v);

      anchor[1]=part->ear[10].align[1];
      anchor[0]=part->ear[10].equiv[anchor[1]];
      n1=align_search(anchor, all_atoms, pdb, r, v, part, &ear1);

      anchor[1]=part->ear[10].align[n];
      anchor[0]=part->ear[10].equiv[anchor[1]];
      n2=align_search(anchor, all_atoms, pdb, r, v, part, &ear2);

      n=align_select(&part->ear[10], &ear1, &ear2);
        
      if(n>2) {
        for(i=1;i<=n;i++) {
          j2=part->ear[10].align[i];
          j1=part->ear[10].equiv[j2];
          calphas[0][i]=all_atoms[0][j1];
          calphas[1][i]=all_atoms[1][j2];
        }
        calphas[0][0]=n;
        calphas[1][0]=n;
        rms_local=rmsd(n, calphas, pdb, r, v);
      }
    }

//      printf("\n TEST 1 rmsd = %6.2f  %ld  %ld  %ld",rms_local,n,n1,n2);

    part->rmsd_align[10]=rms_local;

    j1=0;
    while(j1<10) {
      j=part->best_ind[j1];
      if(n==part->ear[j].align[0] && part->rmsd_align[j]==rms_local) {
        j1=10;
      }
      else if(n>=part->ear[j].align[0] || n>=min_pdb) {
        if(((part->ear[j].align[0]>=min_pdb || n==part->ear[j].align[0]) && 
                                          part->rmsd_align[j]>rms_local) ||
           (n>part->ear[j].align[0] && part->ear[j].align[0]<min_pdb)) {
          k1=part->best_ind[9];
          for(j2=9;j2>j1;j2--) {
            part->best_ind[j2]=part->best_ind[j2-1];
          }
          part->best_ind[j1]=k1;
          part->rmsd_align[k1]=rms_local;
          for(i=0;i<=n;i++) {
            part->ear[k1].align[i]=part->ear[10].align[i];
          }
          for(i=0;i<=all_atoms[1][0];i++) {
            part->ear[k1].equiv[i]=part->ear[10].equiv[i];
          }
          part->ear[k1].equiv[0]=all_atoms[1][0];
          part->ear[k1].align[0]=n;
          j1=10;
        }
      }
      j1++;
    }

    if(n>n_max) n_max=n;
  }

  return n_max;
}

/*-----------------------------------------------------------
/
/   alignment_best_match - alignment search for the best match
/
/------------------------------------------------------------*/
void alignment_best_match(long atoms[2][MAXRES+1], pdb_struct pdb[2],
                          float r[3][3], float v[3], part_struct *part)
{
  long i, j, j1, j2, j12, pdb_min, k1, k2, ms, si, n_max_ms;
  long n, n_min, n_max, n_diff, min_seq, best_match_1, best_match_2;
  lga_workspace *ws = part->ws;
  if (!ws) { /* defensive */
    if (lga_ws_init(&part->ws) != 0) {
      fprintf(stderr, "\n# ERROR! Failed to allocate LGA workspace (alignment_best_match)\n");
      exit(1);
    }
    ws = part->ws;
  }
  long (*l_atoms)[MAXRES+1] = ws ? (long (*)[MAXRES+1])ws->l_atoms2 : NULL;
  long (*ind)[MAXRES+1]     = ws ? (long (*)[MAXRES+1])ws->ind3     : NULL;
  long (*out)[MAXRES+1]     = ws ? (long (*)[MAXRES+1])ws->out2     : NULL;
  long anchor[2], anch_1_0, anch_1_1, anch_n_0, anch_n_1;
  float rms_local, lcs_gdt_ts, lcs_gdt_ts1, lcs_gdt_ts2, rms_lg1, rms_lg2;
  float q_norm, dist_cutoff, nref, n_nref, n_nref1, n_nref2, n_max_nref, dcs;
  part_ear ear[3];

  min_seq=5;  // minimum number of residues in a structure
  ms=15;
  si=0;
  dcs=0.08;   // 0.2  // 0.08

  dist_cutoff=part->dist_cutoff;

/*  Index_table construction to slide first structure along the second
-----------------------------------------------------------------------*/

  j1=pdb[0].n_aa;
  j2=pdb[1].n_aa;

  if(j1<min_seq || j2<min_seq) {
    printf("\n# ERROR! The number of residues in molecule is less then %ld\n",min_seq);
    return;
  }
  if(part->accuracy_opt>1) {
    nref=30.0/(float)j2;  // parameters tuned for the best performance: 20.0 - 40.0
  }
  else {
    nref=40.0/(float)j2;  // parameters tuned for the best performance: 30.0 - 50.0
  }

  atoms[0][0]=j1;
  atoms[1][0]=j2;
  if(j1<j2) {
    n_min=j1;
    n_max=j2;
  }
  else {
    n_min=j2;
    n_max=j1;
  }
  n_diff=n_max-n_min;
  j12=j1+j2;
  pdb_min=n_min;

  i=0;
  for(n=0;n<=n_diff;n++) {  // overlaps in the middle part
    if(j1>=j2) {
      k1=n;
      k2=0;
    }
    else {
      k1=0;
      k2=n;
    }
    if(n_min>=min_seq) {
      i++;
      ind[0][i]=n_min;
      ind[1][i]=k1;
      ind[2][i]=k2;
    }
  }
  for(n=1;n<n_min;n++) {  // overlaps in the terminal parts
    if(n_min-n>=min_seq) {
      i++;
      ind[0][i]=n_min-n;
      ind[1][i]=j1-n_min+n;
      ind[2][i]=0;
      i++;
      ind[0][i]=n_min-n;
      ind[1][i]=0;
      ind[2][i]=j2-n_min+n;
    }
  }
  j12=i;

/*  Using index_table quick search for matches by NB under DIST
-----------------------------------------------------------------*/

  part->dist_cutoff=dist_cutoff*1.25;  // 1.25

  if(part->full_print==1) printf("\nQuick search for maximum N_RES under the distance: %6.2f",part->dist_cutoff);
  n_max=0;
  n_max_ms=2*ms;
  for(i=1;i<=j12;i++) {
    n=ind[0][i];  
    k1=ind[1][i];
    k2=ind[2][i];
    for(j=1;j<=n;j++) {
      l_atoms[0][j]=atoms[0][k1+j];
      l_atoms[1][j]=atoms[1][k2+j];
    }
    l_atoms[0][0]=n;
    l_atoms[1][0]=n;
    n_min=alignment_run(n, l_atoms, atoms, pdb, r, v, part);
 
    j1=l_atoms[0][1];
    j2=l_atoms[1][1];
    if(part->full_print==1) printf("\n i =%4ld   of %4ld   Size =%4ld  from k1 =%4ld k2 =%4ld  N_RES =%4ld",i,j12,n,k1,k2,n_min);
    if(n_min>n_max) {
      n_max=n_min;
    }
    if(n_max>n+n_max_ms) i=j12;
  }

  if(part->full_print==1) printf("\n\nMaximum N_RES =%4ld\n",n_max);
  if(n_max>pdb_min) {
    if(part->error[1]==0) {
      part->error[1]=1;
      printf("\n# WARNING! Check your data. Change of the parameter DIST (-d:f.f) may give you better result.\n");
    }
  }

  part->dist_cutoff=dist_cutoff*1.05;  // 1.05
  best_match_1=1;
  best_match_2=0;
  lcs_gdt_ts1=0.0;
  lcs_gdt_ts2=0.0;
  j=part->best_ind[0];
  n_max=part->ear[j].align[0];
  if(n_max>pdb_min) { 
    n_max_nref=(pdb_min+(float)n_max/pdb_min)*nref; 
    n_max_ms=pdb_min-ms; 
  }
  else { 
    n_max_nref=n_max*nref; 
    n_max_ms=n_max-ms; 
  }
  j=part->best_ind[1];
  n=part->ear[j].align[0];
  
  if(n>=n_max_ms) {
    part->lcs_gdt_print=0;
    lcs_gdt_ts=0.0;
    for(k1=0;k1<10;k1++) {
      j=part->best_ind[k1];
      n=part->ear[j].align[0];
//      if(n>=n_max_ms) {
        for(i=1;i<=n;i++) {
          j2=part->ear[j].align[i];
          j1=part->ear[j].equiv[j2];
          l_atoms[0][i]=atoms[0][j1];
          l_atoms[1][i]=atoms[1][j2];
        }
        l_atoms[0][0]=n;
        l_atoms[1][0]=n;
        lcs_gdt_analysis(n, l_atoms, pdb, r, v, part);
//        lcs_gdt_ts=part->lga_w*part->gdt_ts+(1.0-part->lga_w)*part->lcs_ts-part->rmsd_isp;
        lcs_gdt_ts=part->lga_w*part->gdt_ts+(1.0-part->lga_w)*part->lcs_ts-part->rmsd_isp;
        if(n>pdb_min) { 
          n_nref=(pdb_min+(float)n/pdb_min)*nref; 
          lcs_gdt_ts=lcs_gdt_ts*pdb_min/n;
        }
        else { 
          n_nref=n*nref; 
        }
        if(part->full_print==1) printf("\nN_RES =%4ld   RMSD =%6.2f   GDT =%6.2f   LCS =%6.2f   LGA =%7.3f   NREF =%8.3f",n,part->rmsd_isp,part->gdt_ts,part->lcs_ts,lcs_gdt_ts,n_nref);
        if(n_nref+lcs_gdt_ts>=n_max_nref+lcs_gdt_ts1) {           
          lcs_gdt_ts2=lcs_gdt_ts1;
          lcs_gdt_ts1=lcs_gdt_ts;
          best_match_2=best_match_1;
          best_match_1=k1;
        }       
        else if(n_nref+lcs_gdt_ts>=n_max_nref+lcs_gdt_ts2) {
          lcs_gdt_ts2=lcs_gdt_ts;
          best_match_2=k1;
        }
//      }
    }
    if(part->full_print==1) printf("\n");
  }

  if(part->full_print==1) printf("\nBEST MATCH %ld %ld",best_match_1+1,best_match_2+1);

  part->dist_cutoff=dist_cutoff;

  if(best_match_1!=best_match_2) {
    j=part->best_ind[best_match_2];
    n=part->ear[j].align[0];
    for(i=1;i<=n;i++) {
      j2=part->ear[j].align[i];
      j1=part->ear[j].equiv[j2];
      ear[0].align[i]=j2;
      ear[0].equiv[j2]=j1;
    }
    ear[0].align[0]=n;
  }

  j=part->best_ind[best_match_1];
  n=part->ear[j].align[0];
  for(i=1;i<=n;i++) {
    j2=part->ear[j].align[i];
    j1=part->ear[j].equiv[j2];
    part->ear[10].align[i]=j2;
    part->ear[10].equiv[j2]=j1;
    l_atoms[0][i]=atoms[0][j1];
    l_atoms[1][i]=atoms[1][j2];
  }
  part->ear[10].align[0]=n;
  l_atoms[0][0]=n;
  l_atoms[1][0]=n;
  rms_local=rmsd(n, l_atoms, pdb, r, v);

  n=part->ear[j].align[0];
  anchor[1]=part->ear[j].align[1];
  anchor[0]=part->ear[j].equiv[anchor[1]];
  anch_1_0=anchor[0];
  anch_1_1=anchor[1];
  anchor[1]=part->ear[j].align[n];
  anchor[0]=part->ear[j].equiv[anchor[1]];
  anch_n_0=anchor[0];
  anch_n_1=anchor[1];
  for(j=1;j<10;j++) {
    part->dist_cutoff=dist_cutoff*(0.45+dcs*j);  // 0.0  // 0.45
    anchor[0]=anch_1_0;
    anchor[1]=anch_1_1;
    j1=align_search(anchor, atoms, pdb, r, v, part, &ear[1]);
    anchor[0]=anch_n_0;
    anchor[1]=anch_n_1;
    j2=align_search(anchor, atoms, pdb, r, v, part, &ear[2]);
    n=align_select(&part->ear[j], &ear[1], &ear[2]);
  }

  align_collect(&part->ear[0], &part->ear[1], &part->ear[2], atoms, pdb, r, v, dist_cutoff);
  for(j=1;j<=8;j++) {
    align_collect(&part->ear[j], &part->ear[j-1], &part->ear[j+2], atoms, pdb, r, v, dist_cutoff);
  }
  j=8;

  j1=part->ear[j].align[1];
  j2=part->ear[j].align[2];
  k1=1;
  if(j1==0 || j1>=j2 || part->ear[j].equiv[j1]>=part->ear[j].equiv[j2]) k1=2;
  n=part->ear[j].align[0];
  k2=n;
  for(i=k1;i<=n;i++) {
    k2=i-k1+1;
    j2=part->ear[j].align[i];
    j1=part->ear[j].equiv[j2];
    l_atoms[0][k2]=atoms[0][j1];
    l_atoms[1][k2]=atoms[1][j2];
  }
  n=k2;
  l_atoms[0][0]=n;
  l_atoms[1][0]=n;

  part->lcs_gdt_print=0;
  part->dist_cutoff=dist_cutoff;
  lcs_gdt_analysis(n, l_atoms, pdb, r, v, part);

  n=part->opt_sup1[0][0];
  for(i=0;i<=n;i++) {
    l_atoms[0][i]=part->opt_sup1[0][i];
    l_atoms[1][i]=part->opt_sup2[0][i];
  }
  l_atoms[0][0]=n;
  l_atoms[1][0]=n;

//  lcs_gdt_ts1=part->lga_w*part->gdt_ts+(1.0-part->lga_w)*part->lcs_ts-part->rmsd_isp;
  lcs_gdt_ts1=part->lga_w*part->gdt_ts+(1.0-part->lga_w)*part->lcs_ts;
  rms_lg1=rmsd(n, l_atoms, pdb, r, v);
  if(n>pdb_min) { 
    n_nref1=(pdb_min+(float)n/pdb_min)*nref; 
  }
  else { 
    n_nref1=n*nref; 
  }

  anch_1_0=pdb[0].atom[l_atoms[0][1]].res_n_local;
  anch_1_1=pdb[1].atom[l_atoms[1][1]].res_n_local;
  anch_n_0=pdb[0].atom[l_atoms[0][n]].res_n_local;
  anch_n_1=pdb[1].atom[l_atoms[1][n]].res_n_local;
  for(j=1;j<=10;j++) {
    part->dist_cutoff=dist_cutoff*(0.45+dcs*j);  // 0.0  // 0.45
    anchor[0]=anch_1_0;
    anchor[1]=anch_1_1;
    j1=align_search(anchor, atoms, pdb, r, v, part, &ear[1]);
    anchor[0]=anch_n_0;
    anchor[1]=anch_n_1;
    j2=align_search(anchor, atoms, pdb, r, v, part, &ear[2]);
    n=align_select(&part->ear[j], &ear[1], &ear[2]);
  }

  align_collect(&part->ear[0], &part->ear[1], &part->ear[2], atoms, pdb, r, v, dist_cutoff);
  for(j=1;j<=8;j++) {
    align_collect(&part->ear[j], &part->ear[j-1], &part->ear[j+2], atoms, pdb, r, v, dist_cutoff);
  }
  j=8;

  if(best_match_1!=best_match_2) {
    n=ear[0].align[0];
    for(i=1;i<=n;i++) {
      j2=ear[0].align[i];
      j1=ear[0].equiv[j2];
      part->ear[10].align[i]=j2;
      part->ear[10].equiv[j2]=j1;
      l_atoms[0][i]=atoms[0][j1];
      l_atoms[1][i]=atoms[1][j2];
    }
    part->ear[10].align[0]=n;
    l_atoms[0][0]=n;
    l_atoms[1][0]=n;
    rms_local=rmsd(n, l_atoms, pdb, r, v);
  
    j=8;
    n=part->ear[j].align[0];
    for(i=1;i<=n;i++) {
      j2=part->ear[j].align[i];
      j1=part->ear[j].equiv[j2];
      ear[0].align[i]=j2;
      ear[0].equiv[j2]=j1;
    }
    ear[0].align[0]=n;
  
    j=10;
    n=part->ear[j].align[0];
    anchor[1]=part->ear[j].align[1];
    anchor[0]=part->ear[j].equiv[anchor[1]];
    anch_1_0=anchor[0];
    anch_1_1=anchor[1];
    anchor[1]=part->ear[j].align[n];
    anchor[0]=part->ear[j].equiv[anchor[1]];
    anch_n_0=anchor[0];
    anch_n_1=anchor[1];
    for(j=1;j<10;j++) {
      part->dist_cutoff=dist_cutoff*(0.45+dcs*j);  // 0.0  // 0.45
      anchor[0]=anch_1_0;
      anchor[1]=anch_1_1;
      j1=align_search(anchor, atoms, pdb, r, v, part, &ear[1]);
      anchor[0]=anch_n_0;
      anchor[1]=anch_n_1;
      j2=align_search(anchor, atoms, pdb, r, v, part, &ear[2]);
      n=align_select(&part->ear[j], &ear[1], &ear[2]);
    }
  
    align_collect(&part->ear[0], &part->ear[1], &part->ear[2], atoms, pdb, r, v, dist_cutoff);
    for(j=1;j<=8;j++) {
      align_collect(&part->ear[j], &part->ear[j-1], &part->ear[j+2], atoms, pdb, r, v, dist_cutoff);
    }
    j=8;

    j1=part->ear[j].align[1];
    j2=part->ear[j].align[2];
    k1=1;
    if(j1==0 || j1>=j2 || part->ear[j].equiv[j1]>=part->ear[j].equiv[j2]) k1=2;
    n=part->ear[j].align[0];
    k2=n;
    for(i=k1;i<=n;i++) {
      k2=i-k1+1;
      j2=part->ear[j].align[i];
      j1=part->ear[j].equiv[j2];
      l_atoms[0][k2]=atoms[0][j1];
      l_atoms[1][k2]=atoms[1][j2];
    }
    n=k2;
    l_atoms[0][0]=n;
    l_atoms[1][0]=n;
  
    part->lcs_gdt_print=0;
    part->dist_cutoff=dist_cutoff;
    lcs_gdt_analysis(n, l_atoms, pdb, r, v, part);
  
    n=part->opt_sup1[0][0];
    for(i=0;i<=n;i++) {
      l_atoms[0][i]=part->opt_sup1[0][i];
      l_atoms[1][i]=part->opt_sup2[0][i];
    }
    l_atoms[0][0]=n;
    l_atoms[1][0]=n;
  
//    lcs_gdt_ts2=part->lga_w*part->gdt_ts+(1.0-part->lga_w)*part->lcs_ts-part->rmsd_isp;
    lcs_gdt_ts2=part->lga_w*part->gdt_ts+(1.0-part->lga_w)*part->lcs_ts;
    rms_lg2=rmsd(n, l_atoms, pdb, r, v);
    if(n>pdb_min) { 
      n_nref2=(pdb_min+(float)n/pdb_min)*nref; 
    }
    else { 
      n_nref2=n*nref; 
    }

    anch_1_0=pdb[0].atom[l_atoms[0][1]].res_n_local;
    anch_1_1=pdb[1].atom[l_atoms[1][1]].res_n_local;
    anch_n_0=pdb[0].atom[l_atoms[0][n]].res_n_local;
    anch_n_1=pdb[1].atom[l_atoms[1][n]].res_n_local;
    for(j=1;j<=10;j++) {
      part->dist_cutoff=dist_cutoff*(0.45+dcs*j);  // 0.0  // 0.45
      anchor[0]=anch_1_0;
      anchor[1]=anch_1_1;
      j1=align_search(anchor, atoms, pdb, r, v, part, &ear[1]);
      anchor[0]=anch_n_0;
      anchor[1]=anch_n_1;
      j2=align_search(anchor, atoms, pdb, r, v, part, &ear[2]);
      n=align_select(&part->ear[j], &ear[1], &ear[2]);
    }
  
    align_collect(&part->ear[0], &part->ear[1], &part->ear[2], atoms, pdb, r, v, dist_cutoff);
    for(j=1;j<=8;j++) {
      align_collect(&part->ear[j], &part->ear[j-1], &part->ear[j+2], atoms, pdb, r, v, dist_cutoff);
    }
    j=8;

    if(ear[0].align[0]+n_nref1+lcs_gdt_ts1-rms_lg1 >= part->ear[j].align[0]+n_nref2+lcs_gdt_ts2-rms_lg2) {
      if(part->full_print==1) printf("  Choice1: %ld : %ld %6.2f %6.2f : %ld %6.2f %6.2f",best_match_1+1,ear[0].align[0],rms_lg1,lcs_gdt_ts1,part->ear[j].align[0],rms_lg2,lcs_gdt_ts2);
      n=ear[0].align[0];
      for(i=1;i<=n;i++) {
        j2=ear[0].align[i];
        j1=ear[0].equiv[j2];
        part->ear[j].align[i]=j2;
        part->ear[j].equiv[j2]=j1;
      }
      part->ear[j].align[0]=n;
    }
    else {
      if(part->full_print==1) printf("  Choice2: %ld : %ld %6.2f %6.2f : %ld %6.2f %6.2f",best_match_2+1,ear[0].align[0],rms_lg1,lcs_gdt_ts1,part->ear[j].align[0],rms_lg2,lcs_gdt_ts2);
    }
  }
  else {
    if(part->full_print==1) printf("  Choice3: %ld : %ld %6.2f %6.2f",best_match_1+1,ear[0].align[0],rms_lg1,lcs_gdt_ts1);
  }
  if(part->full_print==1) printf("\n");

  j=8;
  j1=part->ear[j].align[1];
  j2=part->ear[j].align[2];
  k1=1;
  if(j1==0 || j1>=j2 || part->ear[j].equiv[j1]>=part->ear[j].equiv[j2]) k1=2;
  n=part->ear[j].align[0];
  k2=n;
  for(i=k1;i<=n;i++) {
    k2=i-k1+1;
    j2=part->ear[j].align[i];
    j1=part->ear[j].equiv[j2];
    l_atoms[0][k2]=atoms[0][j1];
    l_atoms[1][k2]=atoms[1][j2];
    out[0][k2]=atoms[0][j1];
    out[1][k2]=atoms[1][j2];
  }
  n=k2;
  l_atoms[0][0]=n;
  l_atoms[1][0]=n;

  part->lcs_gdt_print=2;
  part->dist_cutoff=dist_cutoff;

  lcs_gdt_analysis(n, l_atoms, pdb, r, v, part);
  lcs_gdt_ts=part->lga_w*part->gdt_ts+(1.0-part->lga_w)*part->lcs_ts;

  n=part->opt_sup1[0][0];
  for(i=1;i<=n;i++) {
    l_atoms[0][i]=part->opt_sup1[0][i];
    l_atoms[1][i]=part->opt_sup2[0][i];
  }
  l_atoms[0][0]=n;
  l_atoms[1][0]=n;

  rms_local=rmsd(n, l_atoms, pdb, r, v);
  si=check_align_lga(atoms, out, pdb, r, v, part, k2);
  k1=out[0][0];

/*
  printf("\nTEST1: k1, k2, si: %ld %ld %ld \n",k1,k2,si); 
*/

/*
  if(si>0) {
    n=k1;
    k2=n;
    for(i=1;i<=n;i++) {
      l_atoms[0][i]=out[0][i];
      l_atoms[1][i]=out[1][i];
    }
    l_atoms[0][0]=n;
    l_atoms[1][0]=n;

    lcs_gdt_analysis(n, l_atoms, pdb, r, v, part);
    lcs_gdt_ts=part->lga_w*part->gdt_ts+(1.0-part->lga_w)*part->lcs_ts;

    n=part->opt_sup1[0][0];
    for(i=1;i<=n;i++) {
      l_atoms[0][i]=part->opt_sup1[0][i];
      l_atoms[1][i]=part->opt_sup2[0][i];
    }
    l_atoms[0][0]=n;
    l_atoms[1][0]=n;
    rms_local=rmsd(n, l_atoms, pdb, r, v);
    si=check_align_lga(atoms, out, pdb, r, v, part, k2);
    k1=out[0][0];
  }
*/

  if(k1>k2 || si>0) {
    n=k1;
    k2=n;
    for(i=1;i<=n;i++) {
      l_atoms[0][i]=out[0][i];
      l_atoms[1][i]=out[1][i];
    }
    l_atoms[0][0]=n;
    l_atoms[1][0]=n;

    lcs_gdt_analysis(n, l_atoms, pdb, r, v, part);
    lcs_gdt_ts=part->lga_w*part->gdt_ts+(1.0-part->lga_w)*part->lcs_ts;

    n=part->opt_sup1[0][0];
    for(i=1;i<=n;i++) {
      l_atoms[0][i]=part->opt_sup1[0][i];
      l_atoms[1][i]=part->opt_sup2[0][i];

/*
      printf("\nTEST_3 LGA  %4ld atoms %4ld %4ld",i,l_atoms[0][i],l_atoms[1][i]);
*/
    }
    l_atoms[0][0]=n;
    l_atoms[1][0]=n;
    rms_local=rmsd(n, l_atoms, pdb, r, v);
  }

/*
  printf("\nTEST2: k1, k2, si: %ld %ld %ld \n",k1,k2,si);
*/

  if(rms_local<0.0) {
    q_norm=0.0;
  }
  else {
    q_norm=n*0.1/(rms_local+0.1);
  }
  if(rms_local>999 || rms_local<0) {
    lcs_gdt_ts=0.0;
  }
  if(1>n) {
    printf("\n ERROR! Check the number of aligned residues! \n");
    n=0;
    lcs_gdt_ts=0.0;
  }
  part->summary.rms_local=rms_local;
  part->summary.lcs_gdt_ts=lcs_gdt_ts;
  part->summary.q_norm=q_norm;
  part->summary.n=-4;
  
  si=print_lga(atoms, l_atoms, out, pdb, r, v, part, n, k2);
  
  return;
}

/*--------------------------------------------------------------------
/
/   align_select - select the alignment under given superposition
/
/--------------------------------------------------------------------*/
long align_select(part_ear *ear0, part_ear *ear1, part_ear *ear2)
{
  long ok1, ok2, i1, i2, j1, j2, n, n1, n2;

  ok1=0;
  ok2=0;
  n=0;
  n1=ear1->align[0];
  n2=ear2->align[0];
  i1=1;
  i2=1;
  while(ok1==0 || ok2==0) {
    if(ok1==0 && ok2==0) {
      if(ear1->align[i1]==ear2->align[i2]) {
        n++;                  
        j2=ear1->align[i1];
        ear0->align[n]=j2;
        if(ear1->rmsd[j2]>ear2->rmsd[j2]) {
          j1=ear2->equiv[j2];  
        }        
        else {
          j1=ear1->equiv[j2];
        }        
        ear0->equiv[j2]=j1;
        i1++;
        i2++;
        if(i1>n1) ok1=1;
        if(i2>n2) ok2=1;
      }
      else if(ear1->align[i1]<ear2->align[i2]) {
        n++;                 
        j2=ear1->align[i1];
        ear0->align[n]=j2;
        j1=ear1->equiv[j2];
        ear0->equiv[j2]=j1;
        i1++;
        if(i1>n1) ok1=1;
      }
      else {
        n++;                 
        j2=ear2->align[i2];
        ear0->align[n]=j2;
        j1=ear2->equiv[j2];
        ear0->equiv[j2]=j1;
        i2++;
        if(i2>n2) ok2=1;
      }
    }
    else if(ok1==0) {
      while(ok1==0) {
        n++;
        j2=ear1->align[i1];
        ear0->align[n]=j2;
        j1=ear1->equiv[j2];
        ear0->equiv[j2]=j1;
        i1++;
        if(i1>n1) ok1=1;
      }
    }
    else {
      while(ok2==0) {
        n++;
        j2=ear2->align[i2];
        ear0->align[n]=j2;
        j1=ear2->equiv[j2];
        ear0->equiv[j2]=j1;
        i2++;
        if(i2>n2) ok2=1;
      }
    }
    if(i1>n1 && i2>n2) {
      ok1=1;
      ok2=1;
    }
  }
  ear0->align[0]=n;
  return n;
}

/*--------------------------------------------------------------------
/
/   align_collect - collect the alignment from distance levels
/
/--------------------------------------------------------------------*/
void align_collect(part_ear *ear0, part_ear *ear1, part_ear *ear2,
                   long atoms[2][MAXRES+1], pdb_struct pdb[2], 
                   float r[3][3], float v[3], float dist_cutoff)
{
  long i, ok1, ok2, i1, i2, j1, j2, n, n1, n2, prev_1_j1, prev_2_j1;
  long present[MAXRES+1], fixed_1[MAXRES+1], fixed_i[MAXRES+1];

  for(i=0;i<=MAXRES;i++) {
    present[i]=0;
    fixed_1[i]=0;
    fixed_i[i]=0;
  }

  ok1=0;
  ok2=0;
  n=0;
  n1=ear1->align[0];
  n2=ear2->align[0];
  for(i=1;i<=n1;i++) {
    j2=ear1->align[i];
    j1=ear1->equiv[j2];
    fixed_1[j1]=j2;
    fixed_i[i]=j1;
  }
  i1=1;
  i2=1;
  prev_1_j1=0;
  prev_2_j1=0;
  while(ok1==0 || ok2==0) {
    if(ok1==0 && ok2==0) {
      if(ear1->align[i1]==ear2->align[i2]) {
        n++;                  
        j2=ear1->align[i1];
        ear0->align[n]=j2;
        j1=ear1->equiv[j2];
        if(present[j1]>0 || (j1<=prev_1_j1 && j1<=prev_2_j1)) {
          n--;
        }
        else if(j1<=prev_1_j1 && j1>prev_2_j1) {
          n--;
          ear0->align[n]=j2;
          present[j1]=1;
          ear0->equiv[j2]=j1;
          prev_1_j1=j1;
        }
        else {
          prev_2_j1=prev_1_j1;
          prev_1_j1=j1;
          present[j1]=1;
          ear0->equiv[j2]=j1;
        }     
        i1++;
        i2++;
        if(i1>n1) ok1=1;
        if(i2>n2) ok2=1;
      }
      else if(ear1->align[i1]<ear2->align[i2]) {
        n++;                 
        j2=ear1->align[i1];
        ear0->align[n]=j2;
        j1=ear1->equiv[j2];
        if(present[j1]>0 || (j1<=prev_1_j1 && j1<=prev_2_j1)) {
          n--;
        }
        else if(j1<=prev_1_j1 && j1>prev_2_j1) {
          n--;
          ear0->align[n]=j2;
          present[j1]=1;
          ear0->equiv[j2]=j1;
          prev_1_j1=j1;
        }
        else {
          prev_2_j1=prev_1_j1;
          prev_1_j1=j1;
          present[j1]=1;
          ear0->equiv[j2]=j1;
        }     
        i1++;
        if(i1>n1) ok1=1;
      }
      else {
        n++;                 
        j2=ear2->align[i2];
        ear0->align[n]=j2;
        j1=ear2->equiv[j2];
        if(present[j1]>0 || fixed_1[j1]>0 ||
           fixed_i[i1]<j1 || fixed_i[i1-1]>j1 ||
           (j1<=prev_1_j1 && j1<=prev_2_j1)) {
          n--;
        }
        else if(j1<=prev_1_j1 && j1>prev_2_j1) {
          n--;
          ear0->align[n]=j2;
          present[j1]=1;
          ear0->equiv[j2]=j1;
          prev_1_j1=j1;
        }
        else {
          prev_2_j1=prev_1_j1;
          prev_1_j1=j1;
          present[j1]=1;
          ear0->equiv[j2]=j1;
        }     
        i2++;
        if(i2>n2) ok2=1;
      }
    }
    else if(ok1==0) {
      while(ok1==0) {
        n++;
        j2=ear1->align[i1];
        ear0->align[n]=j2;
        j1=ear1->equiv[j2];
        if(present[j1]>0 || (j1<=prev_1_j1 && j1<=prev_2_j1)) {
          n--;
        }
        else if(j1<=prev_1_j1 && j1>prev_2_j1) {
          n--;
          ear0->align[n]=j2;
          present[j1]=1;
          ear0->equiv[j2]=j1;
          prev_1_j1=j1;
        }
        else {
          prev_2_j1=prev_1_j1;
          prev_1_j1=j1;
          present[j1]=1;
          ear0->equiv[j2]=j1;
        }     
        i1++;
        if(i1>n1) ok1=1;
      }
    }
    else {
      while(ok2==0) {
        n++;
        j2=ear2->align[i2];
        ear0->align[n]=j2;
        j1=ear2->equiv[j2];
        if(present[j1]>0 || fixed_1[j1]>0 || fixed_i[i1-1]>j1 ||
           (j1<=prev_1_j1 && j1<=prev_2_j1)) {
          n--;
        }
        else if(j1<=prev_1_j1 && j1>prev_2_j1) {
          n--;
          ear0->align[n]=j2;
          present[j1]=1;
          ear0->equiv[j2]=j1;
          prev_1_j1=j1;
        }
        else {
          prev_2_j1=prev_1_j1;
          prev_1_j1=j1;
          present[j1]=1;
          ear0->equiv[j2]=j1;
        }     
        i2++;
        if(i2>n2) ok2=1;
      }
    }
    if(i1>n1 && i2>n2) {
      ok1=1;
      ok2=1;
    }
  }
  ear0->align[0]=n;

  align_final(ear0, atoms, pdb, r, v, dist_cutoff);

  return;
}

/*--------------------------------------------------------------------
/
/   align_final - select the final alignment
/
/--------------------------------------------------------------------*/
void align_final(part_ear *ear0, long atoms[2][MAXRES+1], pdb_struct pdb[2],
                 float r[3][3], float v[3], float dist_cutoff)
{
  long i, j1, j2, k, k0, k1, k2, n, m, a;
  long align[MAXRES+1], equiv[MAXRES+1], fixed1[MAXRES+1], fixed2[MAXRES+1];
  float x2,y2,z2, d[4], t;
  float xa[MAXRES+1][3], xb[MAXRES+1][3];

  for(i=1;i<=atoms[0][0];i++) {
    xa[i][0]=pdb[0].atom[atoms[0][i]].R.x;
    xa[i][1]=pdb[0].atom[atoms[0][i]].R.y;
    xa[i][2]=pdb[0].atom[atoms[0][i]].R.z;
  }

  for(i=1;i<=atoms[1][0];i++) {
    xb[i][0]=r[0][0]*(pdb[1].atom[atoms[1][i]].R.x - v[0])+
             r[0][1]*(pdb[1].atom[atoms[1][i]].R.y - v[1])+
             r[0][2]*(pdb[1].atom[atoms[1][i]].R.z - v[2]);
    xb[i][1]=r[1][0]*(pdb[1].atom[atoms[1][i]].R.x - v[0])+
             r[1][1]*(pdb[1].atom[atoms[1][i]].R.y - v[1])+
             r[1][2]*(pdb[1].atom[atoms[1][i]].R.z - v[2]);
    xb[i][2]=r[2][0]*(pdb[1].atom[atoms[1][i]].R.x - v[0])+
             r[2][1]*(pdb[1].atom[atoms[1][i]].R.y - v[1])+
             r[2][2]*(pdb[1].atom[atoms[1][i]].R.z - v[2]);
  }                             

  for(i=0;i<=atoms[1][0];i++) {
    fixed1[i]=-999;
    fixed2[i]=-999;
  }

  a=0;
  n=ear0->align[0];
  for(i=1;i<=n;i++) {
    j2=ear0->align[i];
    j1=ear0->equiv[j2];
/*
printf("\n j1=%ld j2=%ld  aj1=%ld aj2=%ld ",j1,j2,atoms[0][j1],atoms[1][j2]);
*/
    for(k=1;k<=3;k++) {
      d[k]=0.0;
      m=0;
      for(k0=1;k0<=3;k0++) {
        k1=j1+k+k0-4;
        k2=j2+k0-2;
        if(k1>=1 && k1<=atoms[0][0] && k2>=1 && k2<=atoms[1][0]) {
          x2=xa[k1][0]-xb[k2][0];
          y2=xa[k1][1]-xb[k2][1];
          z2=xa[k1][2]-xb[k2][2];
          t=sqrt(x2*x2+y2*y2+z2*z2);
          if(k0==2 || t<=dist_cutoff*1.8) {
            m=m+1;
            d[k]=d[k]+t;
          }
        }
      }
      if(m>0) d[k]=d[k]/m;
      else d[k]=999.99;
/*
printf(" d=%8.2f ",d[k]);
*/
    }
    k=0;
    if(d[1]<d[2]) {
      k=-1;
      if(d[3]<d[1]) k=1;
    }
    if(d[3]<d[2]) {
      k=1;
      if(d[1]<d[3]) k=-1;
    }

    for(k0=1;k0<=3;k0++) {
      k1=j1+k+k0-2;
      k2=j2+k0-2;
      if(k1>=1 && k1<=atoms[0][0] && k2>=1 && k2<=atoms[1][0]) {
        x2=xa[k1][0]-xb[k2][0];
        y2=xa[k1][1]-xb[k2][1];
        z2=xa[k1][2]-xb[k2][2];
        t=sqrt(x2*x2+y2*y2+z2*z2);
        if(k0==2 || t<=dist_cutoff*0.8) {
          if(fixed1[a]<k1 && fixed2[a]<k2) {
            a++;
            align[a]=k2;
            equiv[k2]=k1;
            fixed1[a]=k1;
            fixed2[a]=k2;
          }
        }
      }
    }

  }
  ear0->align[0]=a;
  for(i=1;i<=a;i++) {
    j2=align[i];
    j1=equiv[j2];
    ear0->align[i]=j2;
    ear0->equiv[j2]=j1;
/*
printf("\n## j1=%ld j2=%ld  aj1=%ld aj2=%ld ",j1,j2,atoms[0][j1],atoms[1][j2]);
*/
  }
/*
printf("\n");
exit(0);
*/
  return;
}
