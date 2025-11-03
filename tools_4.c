#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "tools.h"

/*-----------------------------------------------------------
/
/   check_align_lga - check LGA records
/
/------------------------------------------------------------*/
long check_align_lga(long atoms[2][MAXRES+1],
               long out[2][MAXRES+1],
               pdb_struct pdb[2], float r[3][3], float v[3],
               part_struct *part, long k2)
{
  long i, j1, j2, j12, n1, n2, al1, al2, n_new, ok, ol;
  long n_max, i1, i2, si;
  long new_atoms[2][MAXRES+1];
  float t, t1, t2, t3, t4, t5, ct1, ct2;

/*
  for(i=1;i<=k2;i++) {
    printf("\nTEST_1 LGA  %4ld atoms %4ld %4ld",i,out[0][i],out[1][i]);
  }
*/

  si=0;
  ct1=1.25;
  ct2=part->dist_cutoff*1.50;
  n_max=atoms[0][0]+atoms[1][0];
  n_new=1;
  n1=1;
  n2=1;
  i=1;
  ok=0;
  ol=0;
  for(j12=1;j12<=n_max;j12++) {
    al1=0;
    al2=0;
    if(out[0][i]==atoms[0][n1] && n1<=atoms[0][0]) al1=1;
    if(out[1][i]==atoms[1][n2] && n2<=atoms[1][0]) al2=1;

    if(al1==1 && al2==1 && i<=k2) {
      j1=out[0][i];
      j2=out[1][i];
      new_atoms[0][n_new]=atoms[0][j1];
      new_atoms[1][n_new]=atoms[1][j2];
      ok=ok+10;

/*
      ol=10;
*/
      if(ol==11) {
        ol=0;
        i1=n1-1;
        i2=n2-2;
        t1=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-1;
        i2=n2-1;
        t2=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1;
        i2=n2-1;
        t3=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1;
        i2=n2;
        t4=distance_calc(atoms, pdb, r, v, i1, i2);

        if(t2*3<t1 || t3*3<t4) {
          if(t2<t3) {
            new_atoms[0][n_new-1]=atoms[0][n1-1];
            new_atoms[1][n_new-1]=atoms[1][n2-1];
          }
          else {
            new_atoms[0][n_new]=atoms[0][n1];
            new_atoms[1][n_new]=atoms[1][n2-1];
          }
          ok=0; ol=0;
          si++;
        }
        else {
          ol=10;
        }
      }
      else if(ol==12) {
        ol=0;
        i1=n1-2;
        i2=n2-1;
        t1=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-1;
        i2=n2-1;
        t2=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-1;
        i2=n2;
        t3=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1;
        i2=n2;
        t4=distance_calc(atoms, pdb, r, v, i1, i2);

        if(t2*3<t1 || t3*3<t4) {
          if(t2<t3) {
            new_atoms[0][n_new-1]=atoms[0][n1-1];
            new_atoms[1][n_new-1]=atoms[1][n2-1];
          }
          else {
            new_atoms[0][n_new]=atoms[0][n1-1];
            new_atoms[1][n_new]=atoms[1][n2];
          }
          ok=0; ol=0;
          si++;
        }
        else {
          ol=0;
        }
      }
      else {
        ol=10;
      }

      n_new++;
      i++;
      n1++;
      n2++;
    }
    else if(al1==1 && n2<=atoms[1][0]) {
      if(ok==11) {
        i1=n1-1;
        i2=n2-2;
        t1=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-1;
        i2=n2-1;
        t2=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-1;
        i2=n2;
        t3=distance_calc(atoms, pdb, r, v, i1, i2);

        if(t1<t2*ct1 || t3<t2*ct1) {
          if(t1<t3) {
            new_atoms[0][n_new-1]=atoms[0][n1-1];
            new_atoms[1][n_new-1]=atoms[1][n2-2];
          }
          else {
            new_atoms[0][n_new-1]=atoms[0][n1-1];
            new_atoms[1][n_new-1]=atoms[1][n2];
          }
          si++;
          ok=0; ol=0;
        }
        else {
          ok=1;
        }
      }
      else if(ok==12) {
        i1=n1-2;
        i2=n2-1;
        t1=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-1;
        i2=n2;
        t2=distance_calc(atoms, pdb, r, v, i1, i2);

        if((t1<=ct2 || n1==3) && (t2<=ct2 || n2==atoms[1][0])) {
          new_atoms[0][n_new-1]=atoms[0][n1-2];
          new_atoms[1][n_new-1]=atoms[1][n2-1];
          new_atoms[0][n_new]=atoms[0][n1-1];
          new_atoms[1][n_new]=atoms[1][n2];
          n_new++;
          si++;
          ok=0; ol=0;
        }
        else {
          ok=1;
        }
      }
      else if(ok==22) {
        i1=n1-3;
        i2=n2-2;
        t1=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-2;
        i2=n2-1;
        t2=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-1;
        i2=n2;
        t3=distance_calc(atoms, pdb, r, v, i1, i2);

        if((t1<=ct2 || n1==4) && t2<=ct2 && (t3<=ct2 || n2==atoms[1][0])) {
          new_atoms[0][n_new-2]=atoms[0][n1-3];
          new_atoms[1][n_new-2]=atoms[1][n2-2];
          new_atoms[0][n_new-1]=atoms[0][n1-2];
          new_atoms[1][n_new-1]=atoms[1][n2-1];
          new_atoms[0][n_new]=atoms[0][n1-1];
          new_atoms[1][n_new]=atoms[1][n2];
          n_new++;
          si++;
          ok=0; ol=0;
        }
        else {
          ok=1;
        }
      }
      else if(ok==32) {
        i1=n1-4;
        i2=n2-3;
        t1=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-3;
        i2=n2-2;
        t2=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-2;
        i2=n2-1;
        t3=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-1;
        i2=n2;
        t4=distance_calc(atoms, pdb, r, v, i1, i2);

        if((t1<=ct2 || n1==5) && t2<=ct2 && t3<=ct2 && (t4<=ct2 || n2==atoms[1][0])) {
          new_atoms[0][n_new-3]=atoms[0][n1-4];
          new_atoms[1][n_new-3]=atoms[1][n2-3];
          new_atoms[0][n_new-2]=atoms[0][n1-3];
          new_atoms[1][n_new-2]=atoms[1][n2-2];
          new_atoms[0][n_new-1]=atoms[0][n1-2];
          new_atoms[1][n_new-1]=atoms[1][n2-1];
          new_atoms[0][n_new]=atoms[0][n1-1];
          new_atoms[1][n_new]=atoms[1][n2];
          n_new++;
          si++;
          ok=0; ol=0;
        }
        else {
          ok=1;
        }
      }
      else if(ok==42) {
        i1=n1-5;
        i2=n2-4;
        t1=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-4;
        i2=n2-3;
        t2=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-3;
        i2=n2-2;
        t3=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-2;
        i2=n2-1;
        t4=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-1;
        i2=n2;
        t5=distance_calc(atoms, pdb, r, v, i1, i2);

        if((t1<=ct2 || n1==6) && t2<=ct2 && t3<=ct2 && t4<=ct2 && (t5<=ct2 || n2==atoms[1][0])) {
          new_atoms[0][n_new-4]=atoms[0][n1-5];
          new_atoms[1][n_new-4]=atoms[1][n2-4];
          new_atoms[0][n_new-3]=atoms[0][n1-4];
          new_atoms[1][n_new-3]=atoms[1][n2-3];
          new_atoms[0][n_new-2]=atoms[0][n1-3];
          new_atoms[1][n_new-2]=atoms[1][n2-2];
          new_atoms[0][n_new-1]=atoms[0][n1-2];
          new_atoms[1][n_new-1]=atoms[1][n2-1];
          new_atoms[0][n_new]=atoms[0][n1-1];
          new_atoms[1][n_new]=atoms[1][n2];
          n_new++;
          si++;
          ok=0; ol=0;
        }
        else {
          ok=1;
        }
      }
      else {
        ok=1;
        if(ol==10) ol=11;
      }
      if(ol!=11) ol=0;
      n2++;
    }
    else if(al2==1 && n1<=atoms[0][0]) {
      if(ok==12) {
        i1=n1-2;
        i2=n2-1;
        t1=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-1;
        i2=n2-1;
        t2=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1;
        i2=n2-1;
        t3=distance_calc(atoms, pdb, r, v, i1, i2);

        if(t1<t2*ct1 || t3<t2*ct1) {
          if(t1<t3) {
            new_atoms[0][n_new-1]=atoms[0][n1-2];
            new_atoms[1][n_new-1]=atoms[1][n2-1];
          }
          else {
            new_atoms[0][n_new-1]=atoms[0][n1];
            new_atoms[1][n_new-1]=atoms[1][n2-1];
          }
          si++;
          ok=0; ol=0;
        }
        else {
          ok=2;
        }
      }
      else if(ok==11) {
        i1=n1-1;
        i2=n2-2;
        t1=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1;
        i2=n2-1;
        t2=distance_calc(atoms, pdb, r, v, i1, i2);

        if((t1<=ct2 || n2==3) && (t2<=ct2 || n1==atoms[0][0])) {
          new_atoms[0][n_new-1]=atoms[0][n1-1];
          new_atoms[1][n_new-1]=atoms[1][n2-2];
          new_atoms[0][n_new]=atoms[0][n1];
          new_atoms[1][n_new]=atoms[1][n2-1];
          n_new++;
          si++;
          ok=0; ol=0;
        }
        else {
          ok=2;
        }
      }
      else if(ok==21) {
        i1=n1-2;
        i2=n2-3;
        t1=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-1;
        i2=n2-2;
        t2=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1;
        i2=n2-1;
        t3=distance_calc(atoms, pdb, r, v, i1, i2);

        if((t1<=ct2 || n2==4) && t2<=ct2 && (t3<=ct2 || n1==atoms[0][0])) {
          new_atoms[0][n_new-2]=atoms[0][n1-2];
          new_atoms[1][n_new-2]=atoms[1][n2-3];
          new_atoms[0][n_new-1]=atoms[0][n1-1];
          new_atoms[1][n_new-1]=atoms[1][n2-2];
          new_atoms[0][n_new]=atoms[0][n1];
          new_atoms[1][n_new]=atoms[1][n2-1];
          n_new++;
          si++;
          ok=0; ol=0;
        }
        else {
          ok=2;
        }
      }
      else if(ok==31) {
        i1=n1-3;
        i2=n2-4;
        t1=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-2;
        i2=n2-3;
        t2=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-1;
        i2=n2-2;
        t3=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1;
        i2=n2-1;
        t4=distance_calc(atoms, pdb, r, v, i1, i2);

        if((t1<=ct2 || n2==5) && t2<=ct2 && t3<=ct2 && (t4<=ct2 || n1==atoms[0][0])) {
          new_atoms[0][n_new-3]=atoms[0][n1-3];
          new_atoms[1][n_new-3]=atoms[1][n2-4];
          new_atoms[0][n_new-2]=atoms[0][n1-2];
          new_atoms[1][n_new-2]=atoms[1][n2-3];
          new_atoms[0][n_new-1]=atoms[0][n1-1];
          new_atoms[1][n_new-1]=atoms[1][n2-2];
          new_atoms[0][n_new]=atoms[0][n1];
          new_atoms[1][n_new]=atoms[1][n2-1];
          n_new++;
          si++;
          ok=0; ol=0;
        }
        else {
          ok=2;
        }
      }
      else if(ok==41) {
        i1=n1-4;
        i2=n2-5;
        t1=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-3;
        i2=n2-4;
        t2=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-2;
        i2=n2-3;
        t3=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1-1;
        i2=n2-2;
        t4=distance_calc(atoms, pdb, r, v, i1, i2);

        i1=n1;
        i2=n2-1;
        t5=distance_calc(atoms, pdb, r, v, i1, i2);

        if((t1<=ct2 || n2==6) && t2<=ct2 && t3<=ct2 && t4<=ct2 && (t5<=ct2 || n1==atoms[0][0])) {
          new_atoms[0][n_new-4]=atoms[0][n1-4];
          new_atoms[1][n_new-4]=atoms[1][n2-5];
          new_atoms[0][n_new-3]=atoms[0][n1-3];
          new_atoms[1][n_new-3]=atoms[1][n2-4];
          new_atoms[0][n_new-2]=atoms[0][n1-2];
          new_atoms[1][n_new-2]=atoms[1][n2-3];
          new_atoms[0][n_new-1]=atoms[0][n1-1];
          new_atoms[1][n_new-1]=atoms[1][n2-2];
          new_atoms[0][n_new]=atoms[0][n1];
          new_atoms[1][n_new]=atoms[1][n2-1];
          n_new++;
          si++;
          ok=0; ol=0;
        }
        else {
          ok=2;
        }
      }
      else {
        ok=2;
        if(ol==10) ol=12;
      }
      if(ol!=12) ol=0;
      n1++;
    }
    else {
      if(n1<=atoms[0][0] && n2<=atoms[1][0]) {
        t=distance_calc(atoms, pdb, r, v, n1, n2);
        if(t<=ct2) {
          new_atoms[0][n_new]=atoms[0][n1];
          new_atoms[1][n_new]=atoms[1][n2];
          n_new++;
          si++;
          ok=ok+10;
          ol=10;
        }
        else {
          ok=0;
        }
      }
      else {
        ok=0;
      }
      n1++;
      n2++;
    }
  }
  for(i=1;i<n_new;i++) {
    out[0][i]=new_atoms[0][i];
    out[1][i]=new_atoms[1][i];

/*
    printf("\nTEST_2 LGA  %4ld atoms %4ld %4ld",i,out[0][i],out[1][i]);
*/
  }
  out[0][0]=n_new-1;
  out[1][0]=n_new-1;

  return si;
}

/*-----------------------------------------------------------
/
/   print_lga - print LGA records
/
/------------------------------------------------------------*/
long print_lga(long atoms[2][MAXRES+1], long l_atoms[2][MAXRES+1],
               long out[2][MAXRES+1],
               pdb_struct pdb[2], float r[3][3], float v[3],
               part_struct *part, long k1, long k2)
{
  long i, j, j12, k, n1, n2, al1, al2, si, rw, n_pl, nb_sl, al_span;
  long n_max, s1, s2, s, s_cutoff, ok_s, si_s, print_s, nb_s, t_s, t_si_s;
  char ch1, ch2, al, span1[MAXRES+1], span2[MAXRES+1];
  char fname_lga[200], sub[MAXRES+1], span1_fix[MAXRES+1], span2_fix[MAXRES+1];
  float t, r_cutoff;
  lga_workspace *ws = part->ws;
  if (!ws) {
    if (lga_ws_init(&part->ws) != 0) {
      fprintf(stderr, "\n*** ERROR! Failed to allocate LGA workspace (print_lga) ***\n");
      exit(1);
    }
    ws = part->ws;
  }
  check_mol2 *mol = ws->mol;
  stext *straline = ws->straline;
  FILE *fp;

  s_cutoff=3;                  // minimum length of the SPAN segment cutoff
  r_cutoff=part->stral_r;      // maximum local RMSD of the lw:1 segment cutoff
  part->resrangaa[0].found=0;  // reset residue range calculation (option "5" or multiple runs)
    
  strcpy(fname_lga,part->fname_lga);
  if(part->stral==0) {
    strcat(fname_lga,".lga");
  }
  else {
    strcat(fname_lga,".stral");
  }
  
  if((fp = fopen(fname_lga,"w"))==NULL) {
    printf("\n*** ERROR! opening file %s for write ***\n",fname_lga);
    exit(0);
  }

  n_pl=0;
  if(part->rw_l>0) {
    rms_window(atoms, l_atoms, out, pdb, part, k1, k2);
    sprintf(part->printlga[n_pl++].line,"\n#      Molecule1      Molecule2  DISTANCE  RMSD(lw:%-ld)",part->rw_l);
    if(part->stral!=0) {
      fprintf(fp,"\n# STRAL results for: %s",part->mname_lga);
      fprintf(fp,"\n# Ord ResName Res    Ord ResName Res   DISTANCE  RMSD(lw:%-ld)",part->rw_l);
    } 
  } 
  else {
    sprintf(part->printlga[n_pl++].line,"\n#      Molecule1      Molecule2  DISTANCE");
  }
  n_max=atoms[0][0]+atoms[1][0];
  n1=1;
  n2=1;
  i=1;
  j=1;
  nb_s=0;
  nb_sl=0;
  si=0;
  si_s=0;
  t_s=0;      // total number of residues in spans
  t_si_s=0;   // total number of identical residues in spans
  print_s=0;
  rw=0;       // counts aligned residues
  ok_s=0;  
  s=0;        // length of the current SPAN
  s1=0;
  s2=0;
  al_span=0;
  for(j12=1;j12<=n_max;j12++) {
    al1=0;
    al2=0;
    if(out[0][i]==atoms[0][n1] && n1<=atoms[0][0]) al1=1;
    if(out[1][i]==atoms[1][n2] && n2<=atoms[1][0]) al2=1;

    al='-';
    if(al1==1 && al2==1 && i<=k2) {
      al_span++;
      if(j<=k1 && out[0][i]==l_atoms[0][j] && out[1][i]==l_atoms[1][j]) {
        j++;
        al='+';  // residues are aligned
      }
      else {
        al='#';  // residues are aligned, but above selected distance cutoff
      }

      if(out[0][i]!=n1) { printf("\nERROR (print_lga) check indexes!\n"); exit(0); }
      if(out[1][i]!=n2) { printf("\nERROR (print_lga) check indexes!\n"); exit(0); }
      C_321(pdb[0].atom[atoms[0][n1]].res_name,&ch1);
      C_321(pdb[1].atom[atoms[1][n2]].res_name,&ch2);

      if(part->resranges==1) {
        if(part->resrangaa[0].found==0) {
          strcpy(part->resrangaa[0].aa,pdb[0].atom[atoms[0][n1]].res_i);
          strcpy(part->resrangaa[2].aa,pdb[1].atom[atoms[1][n2]].res_i);
          part->resrangaa[0].found=1;
          part->resrangaa[2].found=1;
        }
        strcpy(part->resrangaa[1].aa,pdb[0].atom[atoms[0][n1]].res_i);
        strcpy(part->resrangaa[3].aa,pdb[1].atom[atoms[1][n2]].res_i);
        part->resrangaa[1].found=1;
        part->resrangaa[3].found=1;
      }
      
      if(al=='+') {  // residues are aligned
        rw++;
        strcpy(mol[0].equiv_aa[rw].aa,pdb[0].atom[atoms[0][n1]].res_i);
        strcpy(mol[1].equiv_aa[rw].aa,pdb[1].atom[atoms[1][n2]].res_i);
        mol[0].equiv_aa[rw].code1=ch1;
        mol[1].equiv_aa[rw].code1=ch2;
        mol[0].equiv_aa[rw].found=0;
        mol[1].equiv_aa[rw].found=0;
        mol[0].equiv_aa[rw].n_def=0;
        mol[1].equiv_aa[rw].n_def=0;

        t=distance_calc(atoms, pdb, r, v, n1, n2);
/*        
        if(part->sia==1 && !strcmp(pdb[0].atom[atoms[0][n1]].res_i,part->fit_r)) {
          for(k=0;k<3;k++) {
            r[0][k]=part->opt_r2[0][0][k];
            r[1][k]=part->opt_r2[0][1][k];
            r[2][k]=part->opt_r2[0][2][k];
            v[k]=part->opt_v2[0][k];
          }
        }
*/
        if(part->rw_l>0 && part->rw_rms[rw]>=0.0) {
          sprintf(part->printlga[n_pl++].line,"\nLGA    %c %7s      %c %7s   %7.3f   %7.3f",
                  ch1,pdb[0].atom[atoms[0][n1]].res_i,
                  ch2,pdb[1].atom[atoms[1][n2]].res_i,t,part->rw_rms[rw]);
          if(part->stral!=0) {
            fprintf(fp,"\n%4ld  %7s  %c   %4ld  %7s  %c    %7.3f  %7.3f",
                  n1,pdb[0].atom[atoms[0][n1]].res_i,ch1,
                  n2,pdb[1].atom[atoms[1][n2]].res_i,ch2,t,part->rw_rms[rw]);
            if(part->rw_rms[rw]<r_cutoff) {
              ok_s=1;
              span1[s]=ch1;
              span2[s]=ch2;
              s++;
              s1=n1;
              s2=n2;
              if(ch1==ch2) si_s++;
            }
            else if(ok_s==1 && s>=s_cutoff) {
              ok_s=0;
              print_s=1;
            }
            else {
              ok_s=0;
              si_s=0;
              s=0;
            }
          }
        }
        else if(part->rw_l>0) {
          sprintf(part->printlga[n_pl++].line,"\nLGA    %c %7s      %c %7s   %7.3f      -",
                  ch1,pdb[0].atom[atoms[0][n1]].res_i,
                  ch2,pdb[1].atom[atoms[1][n2]].res_i,t);
          if(part->stral!=0) {
            fprintf(fp,"\n%4ld  %7s  %c   %4ld  %7s  %c    %7.3f     -",
                  n1,pdb[0].atom[atoms[0][n1]].res_i,ch1,
                  n2,pdb[1].atom[atoms[1][n2]].res_i,ch2,t);
            if(ok_s==1 && s>=s_cutoff) {
              ok_s=0;
              print_s=1;
            }
            else {
              ok_s=0;
              si_s=0;
              s=0;
            }
          }
        }
        else {
          sprintf(part->printlga[n_pl++].line,"\nLGA    %c %7s      %c %7s   %7.3f",
                  ch1,pdb[0].atom[atoms[0][n1]].res_i,
                  ch2,pdb[1].atom[atoms[1][n2]].res_i,t);
        }
        if(ch1==ch2) si++;
      }
      if(al=='#') {  // residues are aligned, but above the distance cutoff
        if(part->rw_l>0) {
          sprintf(part->printlga[n_pl++].line,"\nLGA    %c %7s      %c %7s   %4c         -",
                  ch1,pdb[0].atom[atoms[0][n1]].res_i,
                  ch2,pdb[1].atom[atoms[1][n2]].res_i,al);
          if(part->stral!=0) {
            fprintf(fp,"\n%4ld  %7s  %c   %4ld  %7s  %c    %4c        -",
                  n1,pdb[0].atom[atoms[0][n1]].res_i,ch1,
                  n2,pdb[1].atom[atoms[1][n2]].res_i,ch2,al);
            if(ok_s==1 && s>=s_cutoff) {
              ok_s=0;
              print_s=1;
            }
            else {
              ok_s=0;
              si_s=0;
              s=0;
            }
          }
        }
        else {
          sprintf(part->printlga[n_pl++].line,"\nLGA    %c %7s      %c %7s   %4c",
                  ch1,pdb[0].atom[atoms[0][n1]].res_i,
                  ch2,pdb[1].atom[atoms[1][n2]].res_i,al);
        }
      }
      i++;
      n1++;
      n2++;
    }
    else if(al1==1 && n2<=atoms[1][0]) {  //  residues that are not aligned
      al_span=0;
      C_321(pdb[1].atom[atoms[1][n2]].res_name,&ch2);
      if(part->rw_l>0) {
        sprintf(part->printlga[n_pl++].line,"\nLGA    -       -      %c %7s      -         -",
                ch2,pdb[1].atom[atoms[1][n2]].res_i);
        if(part->stral!=0) {
          fprintf(fp,"\n   -       -   -   %4ld  %7s  %c       -        -",
                n2,pdb[1].atom[atoms[1][n2]].res_i,ch2);
          if(ok_s==1 && s>=s_cutoff) {
            ok_s=0;
            print_s=1;
          }
          else {
            ok_s=0;
            si_s=0;
            s=0;
          }
        }
      }
      else {
        sprintf(part->printlga[n_pl++].line,"\nLGA    -       -      %c %7s      -",
                ch2,pdb[1].atom[atoms[1][n2]].res_i);
      }
      n2++;
    }
    else if(al2==1 && n1<=atoms[0][0]) {
      al_span=0;
      C_321(pdb[0].atom[atoms[0][n1]].res_name,&ch1);
      if(part->rw_l>0) {
        sprintf(part->printlga[n_pl++].line,"\nLGA    %c %7s      -       -      -         -",
                ch1,pdb[0].atom[atoms[0][n1]].res_i);
        if(part->stral!=0) {
          fprintf(fp,"\n%4ld  %7s  %c      -       -   -       -        -",
                n1,pdb[0].atom[atoms[0][n1]].res_i,ch1);
          if(ok_s==1 && s>=s_cutoff) {
            ok_s=0;
            print_s=1;
          }
          else {
            ok_s=0;
            si_s=0;
            s=0;
          }
        }
      }
      else {
        sprintf(part->printlga[n_pl++].line,"\nLGA    %c %7s      -       -      -",
                ch1,pdb[0].atom[atoms[0][n1]].res_i);
      }
      n1++;
    }
    else {
      al_span=0;
      if(n1<=atoms[0][0] && n2<=atoms[1][0]) { // not aligned residues printed with "-" mark
        C_321(pdb[0].atom[atoms[0][n1]].res_name,&ch1);
        C_321(pdb[1].atom[atoms[1][n2]].res_name,&ch2);
        if(part->rw_l>0) {
          sprintf(part->printlga[n_pl++].line,"\nLGA    %c %7s      %c %7s      -         -",
                  ch1,pdb[0].atom[atoms[0][n1]].res_i,
                  ch2,pdb[1].atom[atoms[1][n2]].res_i);
          if(part->stral!=0) {
            fprintf(fp,"\n%4ld  %7s  %c   %4ld  %7s  %c       -        -",
                  n1,pdb[0].atom[atoms[0][n1]].res_i,ch1,
                  n2,pdb[1].atom[atoms[1][n2]].res_i,ch2);
            if(ok_s==1 && s>=s_cutoff) {
              ok_s=0;
              print_s=1;
            }
            else {
              ok_s=0;
              si_s=0;
              s=0;
            }
          }
        }
        else {
          sprintf(part->printlga[n_pl++].line,"\nLGA    %c %7s      %c %7s      -",
                  ch1,pdb[0].atom[atoms[0][n1]].res_i,
                  ch2,pdb[1].atom[atoms[1][n2]].res_i);
        }
      }
      else if(n1<=atoms[0][0]) {
        C_321(pdb[0].atom[atoms[0][n1]].res_name,&ch1);
        if(part->rw_l>0) {
          sprintf(part->printlga[n_pl++].line,"\nLGA    %c %7s      -       -      -         -",
                  ch1,pdb[0].atom[atoms[0][n1]].res_i);
          if(part->stral!=0) {
            fprintf(fp,"\n%4ld  %7s  %c      -       -   -       -        -",
                  n1,pdb[0].atom[atoms[0][n1]].res_i,ch1);
            if(ok_s==1 && s>=s_cutoff) {
              ok_s=0;
              print_s=1;
            }
            else {
              ok_s=0;
              si_s=0;
              s=0;
            }
          }
        }
        else {
          sprintf(part->printlga[n_pl++].line,"\nLGA    %c %7s      -       -      -",
                  ch1,pdb[0].atom[atoms[0][n1]].res_i);
        }
      }
      else if(n2<=atoms[1][0]) {
        C_321(pdb[1].atom[atoms[1][n2]].res_name,&ch2);
        if(part->rw_l>0) {
          sprintf(part->printlga[n_pl++].line,"\nLGA    -       -      %c %7s      -         -",
                  ch2,pdb[1].atom[atoms[1][n2]].res_i);
          if(part->stral!=0) {
            fprintf(fp,"\n   -       -   -   %4ld  %7s  %c       -        -",
                  n2,pdb[1].atom[atoms[1][n2]].res_i,ch2);
            if(ok_s==1 && s>=s_cutoff) {
              ok_s=0;
              print_s=1;
            }
            else {
              ok_s=0;
              si_s=0;
              s=0;
            }
          }
        }
        else {
          sprintf(part->printlga[n_pl++].line,"\nLGA    -       -      %c %7s      -",
                  ch2,pdb[1].atom[atoms[1][n2]].res_i);
        }
      }
      n1++;
      n2++;
    }
    if(print_s==1) {
// Span is expanded to the full sequence (both ends included) if molecules overlap entirely
      if(al_span==s+2 && 
         ((s1-s == 1 && s1+1 == atoms[0][0]) ||
         (s2-s == 1 && s2+1 == atoms[1][0]))) {
        C_321(pdb[0].atom[atoms[0][s1-s]].res_name,&ch1);
        C_321(pdb[1].atom[atoms[1][s2-s]].res_name,&ch2);
        if(ch1==ch2) si_s++;
        span1_fix[0]=ch1;
        span2_fix[0]=ch2;
        C_321(pdb[0].atom[atoms[0][s1+1]].res_name,&ch1);
        C_321(pdb[1].atom[atoms[1][s2+1]].res_name,&ch2);
        if(ch1==ch2) si_s++;
        span1_fix[s+1]=ch1;
        span2_fix[s+1]=ch2;
        for(k=0;k<s;k++) {
          span1_fix[k+1]=span1[k];
          span2_fix[k+1]=span2[k];
        }
        s=s+2;
        s1=s1+1;
        s2=s2+1;
        for(k=0;k<s;k++) {
          span1[k]=span1_fix[k];
          span2[k]=span2_fix[k];
        }
      }
      sprintf(straline[nb_sl++].line,"\n%-20s   %4ld-%-4ld  %4ld-%-4ld   %7s:%-7s    %7s:%-7s    %4ld   %6.2f",
              part->mname_lga,s1-s+1,s1,s2-s+1,s2,
              pdb[0].atom[atoms[0][s1-s+1]].res_i,pdb[0].atom[atoms[0][s1]].res_i,
              pdb[1].atom[atoms[1][s2-s+1]].res_i,pdb[1].atom[atoms[1][s2]].res_i,
              s,100.0*si_s/s);
      sprintf(straline[nb_sl].line,"\n#1  ");
      for(k=0;k<s;k++) {
        sprintf(sub,"%c",span1[k]);
        strcat(straline[nb_sl].line,sub);
      }
      sprintf(straline[++nb_sl].line,"\n#2  ");
      for(k=0;k<s;k++) {
        sprintf(sub,"%c",span2[k]);
        strcat(straline[nb_sl].line,sub);
      }
//    Calculate ranges from spans
      if(part->resranges==1) {
        if(part->resrangaa[4].found==0) {
          sprintf(part->resrangaa[4].aa,"%ld",s1-s+1);
          sprintf(part->resrangaa[6].aa,"%ld",s2-s+1);
          part->resrangaa[4].found=1;
          part->resrangaa[6].found=1;
        }
        sprintf(part->resrangaa[5].aa,"%ld",s1);
        sprintf(part->resrangaa[7].aa,"%ld",s2);
        part->resrangaa[5].found=1;
        part->resrangaa[7].found=1;
      }

      t_s=t_s+s;
      t_si_s=t_si_s+si_s;
      print_s=0;
      si_s=0;
      s=0;
      nb_sl++;
      nb_s++;
    }
  }
  
  mol[0].n_aa=rw;
  mol[1].n_aa=rw;
  
  if(part->summary.n==-4) {
    if(1>k1) {
      t=0.0;
    }
    else {
      t=100.0*si/k1;
    }
    part->summary.seq_id=t;  
    sprintf(part->summary.s[0].line,"\n#%-8s      N1   N2   DIST      N    RMSD   Seq_Id      LGA_S     LGA_Q ",part->atoms);
    sprintf(part->summary.s[1].line,"\nSUMMARY(LGA) %4ld %4ld   %4.1f   %4ld  %6.2f   %6.2f    %7.3f  %8.3f",
            part->m_n_aa, part->t_n_aa,part->dist_cutoff,k1,part->summary.rms_local,t,
            part->summary.lcs_gdt_ts,part->summary.q_norm);
    part->summary.n=2;  
    if(part->stral!=0) {
      if(t_s>0) t=100.0*t_si_s/t_s;
      else t=0.0;
      sprintf(sub," S_nb  S_N   S_Id ");
      strcat(part->summary.s[0].line,sub);
      sprintf(sub,"  %4ld %4ld %6.2f",nb_s,t_s,t);
      strcat(part->summary.s[1].line,sub);
    }
    if(part->lga_m==1) {
      if(part->t_n_aa>part->m_n_aa) t=part->summary.lcs_gdt_ts*part->t_n_aa/part->m_n_aa;
      else t=part->summary.lcs_gdt_ts;
      sprintf(sub,"   LGA_M ");
      strcat(part->summary.s[0].line,sub);
      sprintf(sub,"  %7.3f ",t);
      strcat(part->summary.s[1].line,sub);
    }
  } 
  sprintf(part->printlga[n_pl++].line,"\n");

// checking atoms within equivalent residues 
  if(part->all_rmsd == 1) {
    part->n_printlga=n_pl;
    t=check_all_atoms(pdb, part, mol, r, v);
    n_pl=part->n_printlga;
  }
  
  for(k=0;k<n_pl;k++) {
    printf("%s",part->printlga[k].line);
  }
  if(part->stral==0) {
    for(k=0;k<n_pl;k++) {
      fprintf(fp,"%s",part->printlga[k].line);
    }
  } 
  
  for(k=0;k<part->summary.n;k++) {
    fprintf(fp,"%s",part->summary.s[k].line);
  }
  if(part->stral==0) {
    fprintf(fp,"\n");
  } 
  else {
    fprintf(fp,"\n\n# STRAL results (SPANS) for: %s",part->mname_lga);
    fprintf(fp,"\n# Domain1.Domain2      Ord1-Ord1  Ord2-Ord2  ResName1:ResName1  ResName2:ResName2  Length  Seq_ID");
    for(k=0;k<nb_sl;k++) {
      fprintf(fp,"%s",straline[k].line);
    }
    fprintf(fp,"\n# END. NB of SPANS: %ld\n",nb_s);
  } 
  fclose(fp);

  return si;
}

/*-----------------------------------------------------------
/
/   check_all_atoms - checks identical residues, 
/                     calculates rmsd on main chain and side chain atoms,
/                     calculates GDC
/
/------------------------------------------------------------*/
float check_all_atoms(pdb_struct pdb[2], part_struct *part, check_mol2 mol[2],
                      float r[3][3], float v[3])
{
  typedef struct {
    char   aa[10];
  } eval_list;
  typedef struct {
    atom_coords coord[201];
  } eval_atoms;

  eval_list eset[MAXRES+1], esup[MAXRES+1];
  eval_atoms evalmol[2];
  long  i, j, k, k1, n, ok, n_atom, m, m1, m2, bug, eset_n, esup_n, open_set, open_sup;
  long i_at_aa, i_at_beg, i_at_end, prev, prev_mc, mis, found;
  long mc_atoms, mc_atoms_common, all_atoms, all_atoms_common, ca_atoms;
  long cd1_0, cd1_1, cd2_0, cd2_1, ce1_0, ce1_1, ce2_0, ce2_1;
  long od1_0, od1_1, od2_0, od2_1, oe1_0, oe1_1, oe2_0, oe2_1;
  long icd1, icd2, ice1, ice2, iod1, iod2, ioe1, ioe2, gdc_naa, tl_gdc[21], tl_gdc_mc[21];
  long gdc_nat, gdc_at_ok, tl_gdc_at[21], gdc_neval, tl_gdc_eval[21], gdcbin1;
  long list0[MAXRES+1], list1[MAXRES+1];
  float t, rl[3][3], vl[3], rg[3][3], vg[3];
  float rmsd_all, rmsd_mc, rmsd_ca, t_dists[21];
  float gdc_max, gdc_mc, gdc_sum, gdc_at, gdc_loc, gdc_loc_mc, gdc_eval, gdcbin_f, gdcbin_s;
  char keyword[200], line[500], curr_res[10], prev_res[10], chcurr[10], sub[10];
  char ch0[10], code1;
  float calc_cb_v, calc_m_v, calc_o_v;
  long calc_in_CB;
  FILE *fp_in;

  calc_cb_v=1.0;
  calc_m_v=0.0;
  calc_o_v=0.0;
  calc_in_CB=0;
  
  for(k=0;k<3;k++) {
    rg[0][k]=r[0][k];
    rg[1][k]=r[1][k];
    rg[2][k]=r[2][k];
    vg[k]=v[k];
  }
  n_atom=0;
  eset_n=0;
  esup_n=0;
  open_set=0;
  open_sup=0;
  ok=0;
  t=0.0;
  gdc_naa=part->t_n_aa;
  gdc_nat=0;
  gdc_neval=0;
  
  gdcbin1=part->gdc_bin + 1;
  gdcbin_f=(float)gdcbin1;
  gdcbin_s=200.0/(gdcbin_f*part->gdc_bin); 

// the lowest allowed distance error (first bin for GDC calculations): t_dists[1]=0.5 or 1.0
  gdc_sum=0.0;  // gdc_sum=0.0 or gdc_sum=0.5;
  for(i=0;i<21;i++) { t_dists[i]=gdc_sum + 0.5*i; } 

  gdc_mc=0.0;
  gdc_sum=0.0;
  gdc_at=-1.0;
  
  strcpy(curr_res,"XXXXX\0");
  strcpy(prev_res,"YYYYY\0");
  if((fp_in = fopen(part->fname_in,"r"))==NULL) {
    printf("\n# ERROR! opening file %s for read\n\n",part->fname_in);
    exit(0);
  }
  while(fgets(line,500,fp_in)!=NULL && ok<3) {
    for(i=0;i<500;i++) if(line[i]=='%') line[i]=' ';
    for(i=0;i<50;i++) keyword[i]=0;
    sscanf(line,"%s",keyword);
    if(!strncmp(keyword,"MOLECULE\0",9)) {
      ok++;
      if(ok<3) {
//        printf(line);
        mol[ok-1].n_atom=n_atom;
        for(n=1;n<=mol[ok-1].n_aa;n++) {
          mol[ok-1].equiv_aa[n].at_eq_cb=-1; // CB atom
        }
      }
      n_atom=0;
    }
    else if(!strncmp(keyword,"ATOM\0",5) ||
            !strncmp(keyword,"HETATM",6)) {

      for(j=0;j<10;j++) sub[j]=' ';
      strncpy(sub,&line[12],4);
      sscanf(sub,"%s",ch0);
      k1=strlen(ch0);       // atom name and length 
      for(j=0;j<10;j++) { sub[j]=' '; chcurr[j]=' '; curr_res[j]=' '; }
//      strncpy(sub,&line[21],1); sscanf(sub,"%s",chcurr); for(j=0;j<10;j++) sub[j]=' ';
      chcurr[0]=line[21]; chcurr[1]='\0';
      strncpy(sub,&line[22],5);
      sscanf(sub,"%s",curr_res);
      for(j=0;j<10;j++) sub[j]=' ';
      if(chcurr[0]!=' ' && strlen(chcurr)>0) {
        sub[0]='_';
        sub[1]=chcurr[0];
        sub[2]='\0';
        strcat(curr_res,sub);
      }
      k=strlen(curr_res);   // current amino-acid number and length (including chain ID)
// printf("# TEST (%ld): %s , %c , %c . \n",k,curr_res,chcurr[0],line[21]);

// BEGIN: Check for selected GDC_eat atoms.
      if(part->gdc_eval.ok>=0) {
        if(ok<3) {  // molecule 1
          for(m=0;m<part->gdc_eval.n;m++) {
            if(ok==2) { m++; } // molecule 2
            if(k == strlen(part->gdc_eval.b[m].aa) &&
               !strcmp(curr_res,part->gdc_eval.b[m].aa)) {
              if(k1 == strlen(part->gdc_eval.e[m].aa) &&
                 !strcmp(ch0,part->gdc_eval.e[m].aa)) {
// printf("# TEST GDC_eat (%ld): %s %s\n",ok,part->gdc_eval.b[m].aa,part->gdc_eval.e[m].aa);
                if(part->gdc_eval.b[m].found==0) {
                  for(j=0;j<10;j++) sub[j]=' ';
                  strncpy(sub,&line[30],8);
                  bug=buger(sub,line,30,37);
                  if(bug>1) {
                    printf("\n# ERROR! Check the molecule %ld and ATOM coordinates: %s\n",ok,line);
                    fclose(fp_in);
                    exit(0);
                  }
                  sscanf(sub,"%f",&evalmol[ok-1].coord[m].R.x);

                  for(j=0;j<10;j++) sub[j]=' ';
                  strncpy(sub,&line[38],8);
                  bug=buger(sub,line,38,45);
                  if(bug>1) {
                    printf("\n# ERROR! Check the molecule %ld and ATOM coordinates: %s\n",ok,line);
                    fclose(fp_in);
                    exit(0);
                  }
                  sscanf(sub,"%f",&evalmol[ok-1].coord[m].R.y);

                  for(j=0;j<10;j++) sub[j]=' ';
                  strncpy(sub,&line[46],8);
                  bug=buger(sub,line,46,53);
                  if(bug>1) {
                    printf("\n# ERROR! Check the molecule %ld and ATOM coordinates: %s\n",ok,line);
                    fclose(fp_in);
                    exit(0);
                  }
                  sscanf(sub,"%f",&evalmol[ok-1].coord[m].R.z);
                }
                part->gdc_eval.b[m].found=1;
                strncpy(part->gdc_eval.b[m].code3,&line[17],3);
              }
            }
            if(ok==1) { m++; } // molecule 1
          }
        }
      }

// BEGIN: Check for evaluation amino acids: gdc_set and gdc_sup. Creates lists.
      if(ok==2 && part->eval>0) {
        for(j=0;j<10;j++) sub[j]=0;
        strncpy(sub,&line[17],3);
        if((k != strlen(prev_res) || strcmp(prev_res,curr_res)) && strncmp(sub,"HOH",3)) {
          if(part->gdc_set.ok>=0) {
            for(j=0;j<part->gdc_set.n;j++) {
              if(k == strlen(part->gdc_set.b[j].aa) &&
                 !strcmp(curr_res,part->gdc_set.b[j].aa)) {
                if(part->gdc_set.b[j].found!=0) {
                  printf("\n# WARNING! Check gdc_set selection and molecule2 ATOM records. Only first range ( BEGIN: %s end: %s ) is used.\n",part->gdc_set.b[j].aa,part->gdc_set.e[j].aa);
                }
                else {
                  part->gdc_set.b[j].found=-1;
                  open_set=1;
// printf("# TEST1 open_set: %2ld %5ld %s %s\n",open_set,j,part->gdc_set.b[j].aa,curr_res);
                }
              }
              if(k == strlen(part->gdc_set.e[j].aa) &&
                 !strcmp(curr_res,part->gdc_set.e[j].aa)) {
// printf("# TEST2 open_set: %2ld %2ld %5ld %s %s %s %s\n",open_set,part->gdc_set.b[j].found,j,part->gdc_set.b[j].aa,part->gdc_set.e[j].aa,prev_res,curr_res);
                if(part->gdc_set.e[j].found!=0) {
                  printf("\n# WARNING! Check gdc_set selection and molecule2 ATOM records. Only first range ( begin: %s END: %s ) is used.\n",part->gdc_set.b[j].aa,part->gdc_set.e[j].aa);
                }
                else {
                  if(part->gdc_set.b[j].found!=-1) {
                    printf("\n# ERROR! Check gdc_set selection and molecule2 ATOM records. Range begin: %s END: %s \n\n",part->gdc_set.b[j].aa,part->gdc_set.e[j].aa);
                    exit(0);
                  }
                  part->gdc_set.b[j].found=1;
                  part->gdc_set.e[j].found=1;
                  open_set=2;
                }
              }
              else if(chcurr[0]!=' ') {
                if(strlen(part->gdc_set.b[j].aa)==1 && strlen(part->gdc_set.e[j].aa)==1 &&
                   part->gdc_set.b[j].aa[0] == part->gdc_set.e[j].aa[0] &&
                   part->gdc_set.b[j].aa[0] == chcurr[0]) {
                  part->gdc_set.b[j].found=1;
                  part->gdc_set.e[j].found=1;
                  open_set=2;
// printf("# TEST3 open_set: %2ld %5ld %s %s\n",open_set,j,part->gdc_set.b[j].aa,curr_res);
                }
              }
            }
            if(open_set>0) {
              if(open_set==2) open_set=0;
              strcpy(eset[eset_n].aa,curr_res);
              eset_n++;
            }
          }
          if(part->gdc_sup.ok>=0) {
            for(j=0;j<part->gdc_sup.n;j++) {
              if(k == strlen(part->gdc_sup.b[j].aa) &&
                 !strcmp(curr_res,part->gdc_sup.b[j].aa)) {
                if(part->gdc_sup.b[j].found!=0) {
                  printf("\n# WARNING! Check gdc_sup selection and molecule2 ATOM records. Only first range ( BEGIN: %s end: %s ) is used.\n",part->gdc_sup.b[j].aa,part->gdc_sup.e[j].aa);
                }
                else {
                  part->gdc_sup.b[j].found=-1;
                  open_sup=1;
// printf("# TEST1 open_sup: %2ld %5ld %s %s\n",open_sup,j,part->gdc_sup.b[j].aa,curr_res);
                }
              }
              if(k == strlen(part->gdc_sup.e[j].aa) &&
                 !strcmp(curr_res,part->gdc_sup.e[j].aa)) {
// printf("# TEST2 open_sup: %2ld %2ld %5ld %s %s %s %s\n",open_sup,part->gdc_sup.b[j].found,j,part->gdc_sup.b[j].aa,part->gdc_sup.e[j].aa,prev_res,curr_res);
                if(part->gdc_sup.e[j].found!=0) {
                  printf("\n# WARNING! Check gdc_sup selection and molecule2 ATOM records. Only first range ( begin: %s END: %s ) is used.\n",part->gdc_sup.b[j].aa,part->gdc_sup.e[j].aa);
                }
                else {
                  if(part->gdc_sup.b[j].found!=-1) {
                    printf("\n# ERROR! Check gdc_sup selection and molecule2 ATOM records. Range begin: %s end: %s \n\n",part->gdc_sup.b[j].aa,part->gdc_sup.e[j].aa);
                    exit(0);
                  }
                  part->gdc_sup.b[j].found=1;
                  part->gdc_sup.e[j].found=1;
                  open_sup=2;
                }
              }
              else if(chcurr[0]!=' ') {
                if(strlen(part->gdc_sup.b[j].aa)==1 && strlen(part->gdc_sup.e[j].aa)==1 &&
                   part->gdc_sup.b[j].aa[0] == part->gdc_sup.e[j].aa[0] &&
                   part->gdc_sup.b[j].aa[0] == chcurr[0]) {
                  part->gdc_sup.b[j].found=1;
                  part->gdc_sup.e[j].found=1;
                  open_sup=2;
//                  printf("# TEST3 open_sup: %2ld %5ld %s %s\n",open_sup,j,part->gdc_sup.b[j].aa,curr_res);
                }
              }
            }
            if(open_sup>0) {
              if(open_sup==2) open_sup=0;
              strcpy(esup[esup_n].aa,curr_res);
              esup_n++;
            }
          }
        }
      }
// END: Check for evaluation amino acids: gdc_set and gdc_sup (create lists)

// BEGIN: Check for equivalent aa and atoms
      for(n=1;n<=mol[0].n_aa;n++) {
// printf("CHECK: %10s : %10s : %c\n",mol[ok-1].equiv_aa[n].aa,curr_res,chcurr[0]);
        if(k == strlen(mol[ok-1].equiv_aa[n].aa) && !strcmp(mol[ok-1].equiv_aa[n].aa,curr_res) &&
           (mol[ok-1].equiv_aa[n].found<1 || (k == strlen(prev_res) && !strcmp(prev_res,curr_res)))) {
// printf("OK: %10s  %10s\n",mol[ok-1].equiv_aa[n].aa,curr_res);
          strcpy(prev_res,curr_res);
          n_atom++;
          
          mol[ok-1].coord[n_atom].id_aa=n;
          code1=mol[ok-1].equiv_aa[n].code1;
          mol[ok-1].coord[n_atom].id_atom=check_aa_atoms(code1, ch0, &j);
          if(mol[ok-1].equiv_aa[n].n_def<j) mol[ok-1].equiv_aa[n].n_def=j; // it also checks OXT atom
          if(mol[ok-1].coord[n_atom].id_atom == 4) mol[ok-1].equiv_aa[n].at_eq_cb=n_atom; // CB atom

          found=0;
          if(mol[ok-1].equiv_aa[n].found==1) {
            i_at_beg=mol[ok-1].equiv_aa[n].at_beg;
            i_at_end=mol[ok-1].equiv_aa[n].at_end;
            for(i=i_at_beg;i<=i_at_end;i++) {
              if(mol[ok-1].coord[i].id_atom==mol[ok-1].coord[n_atom].id_atom) {
                found=1;
                i=i_at_end;
              }
            }
          }
          
// printf("# TEST:  %c %4s %3ld %s",code1, ch0, j, line);
          if(mol[ok-1].coord[n_atom].id_atom<0 || found==1) {
// printf("# TEST-- %c %4s %3ld %s",code1, ch0, j, line);
            n_atom--;
          }
          else {
          
            if(n_atom>=MAXATOMS) {
              printf("\n# ERROR! Check the number of ATOM coordinates. (Limit MAX = %d) \n\n",
                      MAXATOMS);
              fclose(fp_in);
              exit(0);
            }

            if(mol[ok-1].equiv_aa[n].found==0) {
              mol[ok-1].equiv_aa[n].at_beg=n_atom;
              mol[ok-1].equiv_aa[n].found=1;
            }
            mol[ok-1].equiv_aa[n].at_end=n_atom;
            
            for(j=0;j<10;j++) sub[j]=' ';
            strncpy(sub,&line[30],8);
            bug=buger(sub,line,30,37);
            if(bug>1) {
              printf("\n# ERROR! Check the molecule %ld and ATOM coordinates: %s\n",ok,line);
              fclose(fp_in);
              exit(0);
            }
            sscanf(sub,"%f",&mol[ok-1].coord[n_atom].R.x);

            for(j=0;j<10;j++) sub[j]=' ';
            strncpy(sub,&line[38],8);
            bug=buger(sub,line,38,45);
            if(bug>1) {
              printf("\n# ERROR! Check the molecule %ld and ATOM coordinates: %s\n",ok,line);
              fclose(fp_in);
              exit(0);
            }
            sscanf(sub,"%f",&mol[ok-1].coord[n_atom].R.y);

            for(j=0;j<10;j++) sub[j]=' ';
            strncpy(sub,&line[46],8);
            bug=buger(sub,line,46,53);
            if(bug>1) {
              printf("\n# ERROR! Check the molecule %ld and ATOM coordinates: %s\n",ok,line);
              fclose(fp_in);
              exit(0);
            }
            sscanf(sub,"%f",&mol[ok-1].coord[n_atom].R.z);

//            printf(line);
            n=mol[0].n_aa;
          }
        }
        else {
//          printf("# %10s  %10s\n",mol[ok-1].equiv_aa[n].aa,curr_res);
        }
      }
// END: Check for equivalent aa and atoms
    }
    else if(!strncmp(keyword,"END\0",4) ||
            !strncmp(keyword,"ENDMDL",6)) {
      mol[ok-1].n_atom=n_atom;
//      printf(line);
    }
    if(ok>0 && ok<3) {
      mol[ok-1].n_atom=n_atom;
    }
    strcpy(prev_res,curr_res);
  }
  fclose(fp_in);

// checking equivalent atoms 
  prev=0;
  prev_mc=0;
  ca_atoms=0; 
  mc_atoms=0;
  mc_atoms_common=0;
  all_atoms=mol[1].n_atom; 
  all_atoms_common=0; 
  for(n=1;n<=all_atoms;n++) {
    if(mol[1].coord[n].id_atom<4) mc_atoms++;
    i_at_aa=mol[1].coord[n].id_aa;
    i_at_beg=mol[0].equiv_aa[i_at_aa].at_beg;
    i_at_end=mol[0].equiv_aa[i_at_aa].at_end;
    ok=5;    // for different amino-acids or for 'X' the mainchain + CB are considered as equivalent
    if(mol[0].equiv_aa[i_at_aa].code1 == '#' || mol[1].equiv_aa[i_at_aa].code1 == '#') {
      ok=4;  // for unknown amino-acids at least mainchain atoms (N, CA, C, O) are expected
    }
    else if(mol[0].equiv_aa[i_at_aa].code1 == mol[1].equiv_aa[i_at_aa].code1) {
      ok=15; // for identical amino-acids all atoms are taken longo account
    }
//    printf("\n# Check: %4ld %4ld %2ld %c ",i_at_aa,n,mol[1].coord[n].id_atom,mol[1].equiv_aa[i_at_aa].code1);
    if(mol[1].coord[n].id_atom<ok) {
      for(k=i_at_beg;k<=i_at_end;k++) {
        if(mol[1].coord[n].id_atom==mol[0].coord[k].id_atom) {
          i_at_aa=mol[0].coord[k].id_aa;
//          printf(" found: %4ld %4ld %2ld %c ",i_at_aa,k,mol[0].coord[k].id_atom,mol[0].equiv_aa[i_at_aa].code1);
          all_atoms_common++; 
          mol[0].equiv_atom[all_atoms_common]=k;
          mol[1].equiv_atom[all_atoms_common]=n;
          if(prev < i_at_aa) {
            prev_mc=1;
            mol[1].equiv_aa[i_at_aa].at_eq_al_beg=all_atoms_common;
          }
          mol[1].equiv_aa[i_at_aa].at_eq_al_end=all_atoms_common;
          if(mol[0].coord[k].id_atom<4) {
            mc_atoms_common++;
            mol[0].equiv_atom_mc[mc_atoms_common]=k;
            mol[1].equiv_atom_mc[mc_atoms_common]=n;
            if(prev_mc==1) {
              mol[1].equiv_aa[i_at_aa].at_eq_mc_beg=mc_atoms_common;
            }
            prev_mc=0;
            mol[1].equiv_aa[i_at_aa].at_eq_mc_end=mc_atoms_common;
            if(mol[0].coord[k].id_atom==1) {
              ca_atoms++;
              mol[0].equiv_atom_ca[ca_atoms]=k;
              mol[1].equiv_atom_ca[ca_atoms]=n;
              mol[1].equiv_aa[i_at_aa].at_eq_ca=ca_atoms;
            }
          }
          prev=i_at_aa;        // aa coordinates should begin from main chain atoms
        }
      }
    }
  }
  mol[0].equiv_atom[0]=all_atoms_common;
  mol[1].equiv_atom[0]=all_atoms_common;
  mol[0].equiv_atom_mc[0]=mc_atoms_common;
  mol[1].equiv_atom_mc[0]=mc_atoms_common;
  mol[0].equiv_atom_ca[0]=ca_atoms;
  mol[1].equiv_atom_ca[0]=ca_atoms;

/*
for(j=0;j<eset_n;j++) {
  printf("# TEST Eval_set: %5ld %10s\n",j+1,eset[j].aa);
}
for(j=0;j<esup_n;j++) {
  printf("# TEST Eval_sup: %5ld %10s\n",j+1,esup[j].aa);
}
*/
  if(part->gdc_set.ok>=0) {
    if(eset_n==0) {
      printf("\n# ERROR! Empty selection! Check gdc_set and ATOM records from molecule2.\n");
      exit(0);
    }
    else {
      gdc_naa=eset_n;
    }
    for(j=0;j<part->gdc_set.n;j++) {
      if(part->gdc_set.b[j].found!=1 || part->gdc_set.e[j].found!=1) {
        printf("\n# WARNING! Check gdc_set selection and molecule2 ATOM records. Range begin: %s end: %s \n",part->gdc_set.b[j].aa,part->gdc_set.e[j].aa);
      }
    }
  }
  if(part->gdc_sup.ok>=0) {
    for(j=0;j<part->gdc_sup.n;j++) {
      if(part->gdc_sup.b[j].found!=1 || part->gdc_sup.e[j].found!=1) {
        printf("\n# WARNING! Check gdc_sup selection and molecule2 ATOM records. Range begin: %s end: %s \n",part->gdc_sup.b[j].aa,part->gdc_sup.e[j].aa);
      }
    }
    n=0;
    for(k=1;k<=mol[0].n_aa;k++) {
      open_sup=0;
      for(j=0;j<esup_n;j++) {
        if(strlen(esup[j].aa) == strlen(mol[1].equiv_aa[k].aa) &&
           !strcmp(esup[j].aa,mol[1].equiv_aa[k].aa)) {
          open_sup=1;
          j=esup_n;
        }
      }
      if(open_sup==1) {
        n++;
        list0[n]=mol[0].equiv_atom_ca[k];
        list1[n]=mol[1].equiv_atom_ca[k];
      }
    }
    if(esup_n<3 || n<=3) {
      printf("\n# ERROR! Empty selection (%ld %ld ;n<3)! Check gdc_sup and ATOM records from molecule2.\n",esup_n,n);
      exit(0);
    }
// Additional global superposition is calculated using GDC_sup CA atoms
    rmsd_all=rmsd_any(n, list0, list1, mol[0].coord, mol[1].coord, rg, vg);
    if(part->gdc == 1) {
      for(k=0;k<3;k++) {
        r[0][k]=rg[0][k];
        r[1][k]=rg[1][k];
        r[2][k]=rg[2][k];
        v[k]=vg[k];
      }
    }
  }

// checking GDC_eat selection of atoms 
  if(part->gdc_eval.ok>=0) {
    for(i=0;i<21;i++) { tl_gdc_eval[i]=0; }
    ok=0;
    for(m=0;m<part->gdc_eval.n;m++) {
      ok++;
      if(part->gdc_eval.b[m].found==0) {
        printf("\n# WARNING! GDC_eat atom from the molecule%ld: %s.%s was not found\n",ok,part->gdc_eval.b[m].aa,part->gdc_eval.e[m].aa);
      }
      if(ok==2) { ok=0; }
    }
  }

// checking swapping atoms 
  if(part->swap == 1) {
    icd1=0; icd2=0; ice1=0; ice2=0; iod1=0; iod2=0; ioe1=0; ioe2=0;
    printf("\n# Checking swapping");
    for(n=1;n<=mol[0].n_aa;n++) {
      if(mol[0].equiv_aa[n].code1 == mol[1].equiv_aa[n].code1) {
        if(mol[0].equiv_aa[n].code1 == 'F' ||
           mol[0].equiv_aa[n].code1 == 'Y') {
          i_at_beg=mol[1].equiv_aa[n].at_eq_al_beg;
          i_at_end=mol[1].equiv_aa[n].at_eq_al_end;
          cd1_0=0;
          cd1_1=0;
          cd2_0=0;
          cd2_1=0;
          ce1_0=0;
          ce1_1=0;
          ce2_0=0;
          ce2_1=0;
          ok=0;
          for(i=i_at_beg;i<=i_at_end;i++) {
            k=mol[0].equiv_atom[i];
            if(mol[0].coord[k].id_atom == 6) {
              cd1_0=mol[0].equiv_atom[i];
              cd1_1=mol[1].equiv_atom[i];
              icd1=i;
              ok++;
            }
            if(mol[0].coord[k].id_atom == 7) {
              cd2_0=mol[0].equiv_atom[i];
              cd2_1=mol[1].equiv_atom[i];
              icd2=i;
              ok++;
            }
            if(mol[0].coord[k].id_atom == 8) {
              ce1_0=mol[0].equiv_atom[i];
              ce1_1=mol[1].equiv_atom[i];
              ice1=i;
              ok++;
            }
            if(mol[0].coord[k].id_atom == 9) {
              ce2_0=mol[0].equiv_atom[i];
              ce2_1=mol[1].equiv_atom[i];
              ice2=i;
              ok++;
            }
          }
          if(ok==4) {
            t=rot_sqr(mol[0].coord[cd1_0].R, mol[1].coord[cd1_1].R, r, v);
            t=t + rot_sqr(mol[0].coord[cd2_0].R, mol[1].coord[cd2_1].R, r, v);
            t=t - rot_sqr(mol[0].coord[cd2_0].R, mol[1].coord[cd1_1].R, r, v);
            t=t - rot_sqr(mol[0].coord[cd1_0].R, mol[1].coord[cd2_1].R, r, v);
            t=t + rot_sqr(mol[0].coord[ce1_0].R, mol[1].coord[ce1_1].R, r, v);
            t=t + rot_sqr(mol[0].coord[ce2_0].R, mol[1].coord[ce2_1].R, r, v);
            t=t - rot_sqr(mol[0].coord[ce2_0].R, mol[1].coord[ce1_1].R, r, v);
            t=t - rot_sqr(mol[0].coord[ce1_0].R, mol[1].coord[ce2_1].R, r, v);
            if(t>0) {
              printf("\n#   possible swapping detected:  %c %7s      %c %7s",
                        mol[0].equiv_aa[n].code1,mol[0].equiv_aa[n].aa,
                        mol[1].equiv_aa[n].code1,mol[1].equiv_aa[n].aa);
              mol[0].equiv_atom[icd1]=cd2_0;
              mol[0].equiv_atom[icd2]=cd1_0;
              mol[0].equiv_atom[ice1]=ce2_0;
              mol[0].equiv_atom[ice2]=ce1_0;
            }
          }
        }
        if(mol[0].equiv_aa[n].code1 == 'D') {
          i_at_beg=mol[1].equiv_aa[n].at_eq_al_beg;
          i_at_end=mol[1].equiv_aa[n].at_eq_al_end;
          od1_0=0;
          od1_1=0;
          od2_0=0;
          od2_1=0;
          ok=0;
          for(i=i_at_beg;i<=i_at_end;i++) {
            k=mol[0].equiv_atom[i];
            if(mol[0].coord[k].id_atom == 6) {
              od1_0=mol[0].equiv_atom[i];
              od1_1=mol[1].equiv_atom[i];
              iod1=i;
              ok++;
            }
            if(mol[0].coord[k].id_atom == 7) {
              od2_0=mol[0].equiv_atom[i];
              od2_1=mol[1].equiv_atom[i];
              iod2=i;
              ok++;
            }
          }
          if(ok==2) {
            t=rot_sqr(mol[0].coord[od1_0].R, mol[1].coord[od1_1].R, r, v);
            t=t + rot_sqr(mol[0].coord[od2_0].R, mol[1].coord[od2_1].R, r, v);
            t=t - rot_sqr(mol[0].coord[od2_0].R, mol[1].coord[od1_1].R, r, v);
            t=t - rot_sqr(mol[0].coord[od1_0].R, mol[1].coord[od2_1].R, r, v);
            if(t>0) {
              printf("\n#   possible swapping detected:  %c %7s      %c %7s",
                        mol[0].equiv_aa[n].code1,mol[0].equiv_aa[n].aa,
                        mol[1].equiv_aa[n].code1,mol[1].equiv_aa[n].aa);
              mol[0].equiv_atom[iod1]=od2_0;
              mol[0].equiv_atom[iod2]=od1_0;
            }
          }
        }
        if(mol[0].equiv_aa[n].code1 == 'E') {
          i_at_beg=mol[1].equiv_aa[n].at_eq_al_beg;
          i_at_end=mol[1].equiv_aa[n].at_eq_al_end;
          oe1_0=0;
          oe1_1=0;
          oe2_0=0;
          oe2_1=0;
          ok=0;
          for(i=i_at_beg;i<=i_at_end;i++) {
            k=mol[0].equiv_atom[i];
            if(mol[0].coord[k].id_atom == 7) {
              oe1_0=mol[0].equiv_atom[i];
              oe1_1=mol[1].equiv_atom[i];
              ioe1=i;
              ok++;
            }
            if(mol[0].coord[k].id_atom == 8) {
              oe2_0=mol[0].equiv_atom[i];
              oe2_1=mol[1].equiv_atom[i];
              ioe2=i;
              ok++;
            }
          }
          if(ok==2) {
            t=rot_sqr(mol[0].coord[oe1_0].R, mol[1].coord[oe1_1].R, r, v);
            t=t + rot_sqr(mol[0].coord[oe2_0].R, mol[1].coord[oe2_1].R, r, v);
            t=t - rot_sqr(mol[0].coord[oe2_0].R, mol[1].coord[oe1_1].R, r, v);
            t=t - rot_sqr(mol[0].coord[oe1_0].R, mol[1].coord[oe2_1].R, r, v);
            if(t>0) {
              t=sqrt(t/4.0);
              printf("\n#   possible swapping detected:  %c %7s      %c %7s",
                        mol[0].equiv_aa[n].code1,mol[0].equiv_aa[n].aa,
                        mol[1].equiv_aa[n].code1,mol[1].equiv_aa[n].aa);
              mol[0].equiv_atom[ioe1]=oe2_0;
              mol[0].equiv_atom[ioe2]=oe1_0;
            }
          }
        }
      }
    }
    printf("\n");
  }

// BEGIN: Local and Global calculations per amino acids
  k=0;
  for(j=0;j<part->n_printlga;j++) {
    if(strstr(part->printlga[j].line,"  Molecule1  ")) {
      if(part->rw_l==0) {
        strcat(part->printlga[j].line,"    Mis    MC     All ");
      }
      else {
        strcat(part->printlga[j].line,"  Mis    MC     All");
      }
      if(part->eval > 0) {
        strcat(part->printlga[j].line,"   Dist_max   GDC_mc  GDC_all");
      }
      if(part->gdc_at.ok>=0) {
        strcat(part->printlga[j].line,"  Dist_at");
// Checks GDC_at parameters
        for(i=0;i<21;i++) { tl_gdc_at[i]=0; }
        for(m=0;m<part->gdc_at.n;m++) {
          part->gdc_at.e[m].n_def=check_aa_atoms(part->gdc_at.b[m].aa[0], part->gdc_at.e[m].aa, &m1);
          if(part->gdc_at.e[m].n_def<0) {
            if(part->gdc_at.b[m].aa[0]=='G' && !strncmp(part->gdc_at.e[m].aa,"CB\0",3)) {
              part->gdc_at.e[m].n_def=5;
            }
            else if(part->gdc_at_mcb<0) {
              printf("\n# ERROR! Check the \"-gdc_at\" option and the aaname.atom definition: %s.%s\n",part->gdc_at.b[m].aa,part->gdc_at.e[m].aa);
              exit(0);
            }
          }
        }
      }
    }
    if(strstr(part->printlga[j].line,"LGA ")) {
// Local RMSDs calculated per aligned amino acids
      if(part->printlga[j].line[38] == '#' || part->printlga[j].line[38] == '-') {
        if(part->rw_l==0) strcat(part->printlga[j].line,"        -     -       -         -        -        -");
        else strcat(part->printlga[j].line,"        -     -       -        -        -        -");
        if(part->gdc_at.ok>=0) { strcat(part->printlga[j].line,"        -"); }
      }
      else {
        k++;
        strcpy(keyword,"  \0");
        if(part->gdc_ref>0) {
          mol[0].equiv_aa[k].n_def=0;
          for(i=mol[1].equiv_aa[k].at_beg;i<=mol[1].equiv_aa[k].at_end;i++) {
            if(mol[1].coord[i].id_atom<4) mol[0].equiv_aa[k].n_def++;
          }
          mol[1].equiv_aa[k].n_def=mol[1].equiv_aa[k].at_end-mol[1].equiv_aa[k].at_beg+1;
          if(part->gdc_ref>1 && mol[0].equiv_aa[k].code1 != mol[1].equiv_aa[k].code1) {
            mol[1].equiv_aa[k].n_def=0;
            if(mol[0].equiv_aa[k].code1 == 'G' || mol[1].equiv_aa[k].code1 == 'G' ||
               mol[0].equiv_aa[k].code1 == '#' || mol[1].equiv_aa[k].code1 == '#') m=4;
            else m=5;
            for(i=mol[1].equiv_aa[k].at_beg;i<=mol[1].equiv_aa[k].at_end;i++) {
              if(mol[1].coord[i].id_atom<m) mol[1].equiv_aa[k].n_def++;
            }
          }
          if(mol[0].equiv_aa[k].n_def<1) {
            mol[0].equiv_aa[k].n_def=mol[1].equiv_aa[k].n_def;
          }
        }
        else {
          if(mol[1].equiv_aa[k].n_def<=4) {
            mol[0].equiv_aa[k].n_def=mol[1].equiv_aa[k].n_def;
          }
          else {
            mol[0].equiv_aa[k].n_def=4;
          }
        }
        mis=mol[1].equiv_aa[k].n_def-
            (mol[1].equiv_aa[k].at_eq_al_end-mol[1].equiv_aa[k].at_eq_al_beg+1);
        if(part->printlga[j].line[48] == '-') strcat(part->printlga[j].line,"   ");
        n=mol[1].equiv_aa[k+part->rw_l].at_eq_mc_end-mol[1].equiv_aa[k-part->rw_l].at_eq_mc_beg+1;
        if(n >= 3 && k > part->rw_l && k <= ca_atoms - part->rw_l) {
          for(i=0;i<n;i++) {
            list0[i+1]=mol[0].equiv_atom_mc[mol[1].equiv_aa[k-part->rw_l].at_eq_mc_beg+i];
            list1[i+1]=mol[1].equiv_atom_mc[mol[1].equiv_aa[k-part->rw_l].at_eq_mc_beg+i];
// printf("# Test1: %4ld %4ld %4ld %4ld %2ld %c %8.3f \n",n,i+1,k,list0[i+1],mol[0].coord[list0[i+1]].id_atom,mol[0].equiv_aa[k].code1,mol[0].coord[list0[i+1]].R.x);
// printf("# Test2: %4ld %4ld %4ld %4ld %2ld %c %8.3f \n",n,i+1,k,list1[i+1],mol[1].coord[list1[i+1]].id_atom,mol[1].equiv_aa[k].code1,mol[1].coord[list1[i+1]].R.x);
          }
          rmsd_mc=rmsd_any(n, list0, list1, mol[0].coord, mol[1].coord, rl, vl);
          n=mol[1].equiv_aa[k+part->rw_l].at_eq_al_end-mol[1].equiv_aa[k-part->rw_l].at_eq_al_beg+1;
          for(i=0;i<n;i++) {
            list0[i+1]=mol[0].equiv_atom[mol[1].equiv_aa[k-part->rw_l].at_eq_al_beg+i];
            list1[i+1]=mol[1].equiv_atom[mol[1].equiv_aa[k-part->rw_l].at_eq_al_beg+i];
          }
          rmsd_all=rmsd_any(n, list0, list1, mol[0].coord, mol[1].coord, rl, vl);
          if(part->rw_l==0) sprintf(keyword,"    %2ld   %6.3f  %6.3f ",mis,rmsd_mc,rmsd_all);
          else sprintf(keyword,"    %2ld   %6.3f  %6.3f",mis,rmsd_mc,rmsd_all);
        }
        else {
          if(part->rw_l==0) sprintf(keyword,"    %2ld     -       -    ",mis);
          else sprintf(keyword,"    %2ld     -       -   ",mis);
        }
        strcat(part->printlga[j].line,keyword);
// GDC - global distance calculation: distance based evaluation of selected amino acids
        if(part->eval > 0) {
          gdc_at=-1.0;
          open_set=0;
          if(part->gdc_set.ok>=0 && eset_n>0) {
            for(i=0;i<eset_n;i++) {
              if(strlen(eset[i].aa) == strlen(mol[1].equiv_aa[k].aa) &&
                 !strcmp(eset[i].aa,mol[1].equiv_aa[k].aa)) {
                open_set=1;
                i=eset_n;
              }
            }
          }
          else {
            open_set=1;
          }
          if(open_set==1) {
            n=mol[1].equiv_aa[k].at_eq_al_end-mol[1].equiv_aa[k].at_eq_al_beg+1;
// GDC_mc (gdc_mc <- mainchain atoms), GDC_all (gdc <- all atoms), GDC_at calculations
            gdc_loc=0.0;
            gdc_loc_mc=0.0;
            gdc_max=0.0;
            gdc_at_ok=0;
            for(i=0;i<21;i++) { tl_gdc[i]=0; tl_gdc_mc[i]=0; }
            for(i=0;i<n;i++) {
              list0[i+1]=mol[0].equiv_atom[mol[1].equiv_aa[k].at_eq_al_beg+i];
              list1[i+1]=mol[1].equiv_atom[mol[1].equiv_aa[k].at_eq_al_beg+i];
              t=gdc_dist(list0[i+1], list1[i+1], mol[0].coord, mol[1].coord, rg, vg);
// GDC_at calculation: distances between corresponding aaname.atoms
              if(part->gdc_at.ok>=0 && gdc_at_ok==0) {
                if(part->gdc_at_mcb >= 0) {
                  if(part->gdc_at_mcb == mol[1].coord[list1[i+1]].id_atom) {
                    gdc_at_ok=1;
                  }
                  else if(part->gdc_at_mcb == 4 && i == 3 &&
                          (mol[0].equiv_aa[k].code1 == 'G' || mol[1].equiv_aa[k].code1 == 'G')) {
// printf("\n# TEST! CB requested  = %c  %c",mol[0].equiv_aa[k].code1,mol[1].equiv_aa[k].code1);
                    gdc_at_ok=2;
                  }
                }
                else if(mol[0].equiv_aa[k].code1 == mol[1].equiv_aa[k].code1) {
                  for(m=0;m<part->gdc_at.n;m++) {
                    if(mol[1].equiv_aa[k].code1 == part->gdc_at.b[m].aa[0]) {
                      if(mol[1].coord[list1[i+1]].id_atom == part->gdc_at.e[m].n_def) {
                        gdc_at_ok=1;
                      }
                      else if(part->gdc_at.b[m].aa[0] == 'G' && 
                              part->gdc_at.e[m].n_def == 5 && i == 3) {
                        gdc_at_ok=2;
// printf("# TEST GDC_at: %5ld %5ld %5ld %c %s %s %s %ld\n",m,n,gdc_nat+1,mol[0].equiv_aa[k].code1,mol[1].equiv_aa[k].aa,part->gdc_at.b[m].aa,part->gdc_at.e[m].aa,part->gdc_at.e[m].n_def);
                        m=part->gdc_at.n;
                      }
                    }
                  }
                }
                if(gdc_at_ok>=1) {
                  if(gdc_at_ok==2) {
                    if(mol[0].equiv_aa[k].code1 == 'G' || mol[0].equiv_aa[k].at_eq_cb < 0) {
// CB calculations for GLY from molecule1
                      for(m1=0;m1<4;m1++) {
                        m2=mol[0].coord[list0[m1+1]].id_atom;
                        pdb[1].atom[MAXRES-m2].R.x=mol[0].coord[list0[m1+1]].R.x;
                        pdb[1].atom[MAXRES-m2].R.y=mol[0].coord[list0[m1+1]].R.y;
                        pdb[1].atom[MAXRES-m2].R.z=mol[0].coord[list0[m1+1]].R.z;
                      }
                      calc_cb(calc_cb_v, calc_m_v, calc_o_v, pdb, calc_in_CB);
                      mol[0].coord[MAXATOMS].R.x=pdb[0].atom[MAXRES-1].R.x;
                      mol[0].coord[MAXATOMS].R.y=pdb[0].atom[MAXRES-1].R.y;
                      mol[0].coord[MAXATOMS].R.z=pdb[0].atom[MAXRES-1].R.z;
/*
                printf("\n# TEST! Completed: ATOM N  = %f",pdb[1].atom[MAXRES].R.x);
                printf("\n# TEST! Completed: ATOM CA = %f",pdb[1].atom[MAXRES-1].R.x);
                printf("\n# TEST! Completed: ATOM C  = %f",pdb[1].atom[MAXRES-2].R.x);
                printf("\n# TEST! Completed: ATOM O  = %f",pdb[1].atom[MAXRES-3].R.x);
                printf("\n# TEST! Completed: ATOM CB = %f\n",mol[0].coord[MAXATOMS].R.x);
*/
                    }
                    else {
                      m1=mol[0].equiv_aa[k].at_eq_cb;
                      mol[0].coord[MAXATOMS].R.x=mol[0].coord[m1].R.x;
                      mol[0].coord[MAXATOMS].R.y=mol[0].coord[m1].R.y;
                      mol[0].coord[MAXATOMS].R.z=mol[0].coord[m1].R.z;
                    }
                    if(mol[1].equiv_aa[k].code1 == 'G' || mol[1].equiv_aa[k].at_eq_cb < 0) {
// CB calculations for GLY from molecule2
                      for(m1=0;m1<4;m1++) {
                        m2=mol[1].coord[list1[m1+1]].id_atom;
                        pdb[1].atom[MAXRES-m2].R.x=mol[1].coord[list1[m1+1]].R.x;
                        pdb[1].atom[MAXRES-m2].R.y=mol[1].coord[list1[m1+1]].R.y;
                        pdb[1].atom[MAXRES-m2].R.z=mol[1].coord[list1[m1+1]].R.z;
                      }
                      calc_cb(calc_cb_v, calc_m_v, calc_o_v, pdb, calc_in_CB);
                      mol[1].coord[MAXATOMS].R.x=pdb[0].atom[MAXRES-1].R.x;
                      mol[1].coord[MAXATOMS].R.y=pdb[0].atom[MAXRES-1].R.y;
                      mol[1].coord[MAXATOMS].R.z=pdb[0].atom[MAXRES-1].R.z;
/*
                printf("\n# TEST! Completed: ATOM N  = %f",pdb[1].atom[MAXRES].R.x);
                printf("\n# TEST! Completed: ATOM CA = %f",pdb[1].atom[MAXRES-1].R.x);
                printf("\n# TEST! Completed: ATOM C  = %f",pdb[1].atom[MAXRES-2].R.x);
                printf("\n# TEST! Completed: ATOM O  = %f",pdb[1].atom[MAXRES-3].R.x);
                printf("\n# TEST! Completed: ATOM CB = %f\n",mol[1].coord[MAXATOMS].R.x);
*/                  
                    }
                    else {
                      m1=mol[1].equiv_aa[k].at_eq_cb;
                      mol[1].coord[MAXATOMS].R.x=mol[1].coord[m1].R.x;
                      mol[1].coord[MAXATOMS].R.y=mol[1].coord[m1].R.y;
                      mol[1].coord[MAXATOMS].R.z=mol[1].coord[m1].R.z;
                    }
                    
                    m1=MAXATOMS;
                    gdc_at=gdc_dist(m1, m1, mol[0].coord, mol[1].coord, rg, vg);
                    gdc_at_ok=1;
                  }
                  else {
                    gdc_at=t;
                  }
                  gdc_nat++;
                  for(m1=1;m1<21;m1++) {
                    if(gdc_at<=t_dists[m1]) {
                      tl_gdc_at[m1]++;
                    }
                  }
                }
              }
// GDC_max calculation: maximum distances between corresponding atoms within selected amino-acids
              if(t>gdc_max) {
                gdc_max=t;
              }
// GDC_mc (gdc_mc <- mainchain atoms) and GDC_all (gdc <- all atoms) calculations
              for(m=1;m<21;m++) {
                if(t<=t_dists[m]) {
                  tl_gdc[m]++;
                  if(mol[1].coord[list1[i+1]].id_atom<4) {
                    tl_gdc_mc[m]++;
                  }
                }
              }
            }
            sprintf(keyword," %8.3f",gdc_max);
            strcat(part->printlga[j].line,keyword);
            for(m=1;m<gdcbin1;m++) {               // counting local GDC_mc and GDC_all
              gdc_loc=gdc_loc + (float)(gdcbin_f-m)*tl_gdc[m];
              gdc_loc_mc=gdc_loc_mc + (float)(gdcbin_f-m)*tl_gdc_mc[m];
            }
            gdc_loc=gdc_loc*gdcbin_s/mol[1].equiv_aa[k].n_def;
            gdc_loc_mc=gdc_loc_mc*gdcbin_s/mol[0].equiv_aa[k].n_def;
            gdc_sum=gdc_sum+gdc_loc;
            gdc_mc=gdc_mc+gdc_loc_mc;
            sprintf(keyword," %8.3f %8.3f",gdc_loc_mc,gdc_loc);
            strcat(part->printlga[j].line,keyword);
          }
          else {
            sprintf(keyword,"     -        -        -   ");
            strcat(part->printlga[j].line,keyword);
          }
          if(part->gdc_at.ok>=0) {
            if(gdc_at>=0) {
              sprintf(keyword," %8.3f",gdc_at);
            }
            else {
              sprintf(keyword,"     -");
            }
            strcat(part->printlga[j].line,keyword);
          }
        }
      }
    }
  }
// END: Local and Global calculations per amino acids

  j=part->n_printlga;
// Calculating selected GDC_eat distances and a SUMMARY value 
  if(part->gdc_eval.ok>=0) {
    gdc_neval=0;  // number of evaluated pairs
    m1=0;
    ok=0;
    for(m=0;m<part->gdc_eval.n;m++) {
      ok++;
      if(part->gdc_eval.b[m].found==0) {
        m1=ok;
      }
      if(ok==2) {
        gdc_neval++;
        if(m1==0) {
          t=gdc_dist(m-1, m, evalmol[0].coord, evalmol[1].coord, rg, vg);
          for(i=1;i<21;i++) {
            if(t<=t_dists[i]) {
              tl_gdc_eval[i]++;
            }
          }
          sprintf(part->printlga[j++].line,"\nGDC_eat:   %3s %8s.%-4s    %3s %8s.%-4s  distance: %8.3f",part->gdc_eval.b[m-1].code3,part->gdc_eval.b[m-1].aa,part->gdc_eval.e[m-1].aa,part->gdc_eval.b[m].code3,part->gdc_eval.b[m].aa,part->gdc_eval.e[m].aa,t);
        }
        else {
          sprintf(part->printlga[j++].line,"\nGDC_eat:   %3s %8s.%-4s    %3s %8s.%-4s  distance:     - ",part->gdc_eval.b[m-1].code3,part->gdc_eval.b[m-1].aa,part->gdc_eval.e[m-1].aa,part->gdc_eval.b[m].code3,part->gdc_eval.b[m].aa,part->gdc_eval.e[m].aa);
        }
        ok=0; m1=0; 
      }
    }
    sprintf(part->printlga[j++].line,"\n");
  }

  gdc_mc=gdc_mc/gdc_naa;
  gdc_sum=gdc_sum/gdc_naa;
  if(mc_atoms==0) mc_atoms=1;
  if(all_atoms==0) all_atoms=1;
  sprintf(part->printlga[j++].line,"\n# RMSD_GDC results:       CA      MC common percent     ALL common percent   GDC_mc  GDC_all");
  sprintf(part->printlga[j++].line,"\nNUMBER_OF_ATOMS_AA:    %5ld   %5ld  %5ld %7.2f   %5ld  %5ld %7.2f             %5ld",
          ca_atoms,mc_atoms,mc_atoms_common,100.0*mc_atoms_common/mc_atoms,
          all_atoms,all_atoms_common,100.0*all_atoms_common/all_atoms,gdc_naa);
  
  rmsd_all=rmsd_any(all_atoms_common, mol[0].equiv_atom, mol[1].equiv_atom,
                    mol[0].coord, mol[1].coord, rl, vl);
  rmsd_mc=rmsd_any(mc_atoms_common, mol[0].equiv_atom_mc, mol[1].equiv_atom_mc,
                    mol[0].coord, mol[1].coord, rl, vl);
  rmsd_ca=rmsd_any(ca_atoms, mol[0].equiv_atom_ca, mol[1].equiv_atom_ca,
                    mol[0].coord, mol[1].coord, rl, vl);
  sprintf(part->printlga[j++].line,"\nSUMMARY(RMSD_GDC):  %8.3f       %8.3f               %8.3f         %8.3f %8.3f",rmsd_ca,rmsd_mc,rmsd_all,gdc_mc,gdc_sum);

  if(part->gdc_at.ok>=0) {
    if(part->gdc_at_mcb >= 0) { gdc_nat=gdc_naa; } // normalize to GDC_mc and GDC_all results
    gdc_at=0.0;
    if(gdc_nat>0) {
      for(m=1;m<gdcbin1;m++) {
        gdc_at=gdc_at + (float)(gdcbin_f-m)*tl_gdc_at[m];
      }
      gdc_at=gdc_at*gdcbin_s/gdc_nat;
    }
    strcat(part->printlga[j-3].line,"   GDC_at");
    sprintf(keyword,"    %5ld",gdc_nat);
    strcat(part->printlga[j-2].line,keyword);
    sprintf(keyword," %8.3f",gdc_at);
    strcat(part->printlga[j-1].line,keyword);
  }
  if(part->gdc_eval.ok>=0) {
    gdc_eval=0.0;
    if(gdc_neval>0) {
      for(m=1;m<gdcbin1;m++) {
        gdc_eval=gdc_eval + (float)(gdcbin_f-m)*tl_gdc_eval[m];
      }
      gdc_eval=gdc_eval*gdcbin_s/gdc_neval;
    }
    strcat(part->printlga[j-3].line,"  GDC_eat");
    sprintf(keyword,"    %5ld",gdc_neval);
    strcat(part->printlga[j-2].line,keyword);
    sprintf(keyword," %8.3f",gdc_eval);
    strcat(part->printlga[j-1].line,keyword);
  }
  strcat(part->printlga[j-1].line,"\n");
  part->n_printlga=j;
  
  return t;
}

/*-----------------------------------------------------------
/
/   rot_sqr - rotates first coordinates and calculates distance (sqr) residue window
/
/------------------------------------------------------------*/
float rot_sqr(vector c1, vector c2, float r[3][3], float v[3])
{
  float x, y, z, t;
  x=c2.x-(r[0][0]*c1.x + r[1][0]*c1.y + r[2][0]*c1.z + v[0]);
  y=c2.y-(r[0][1]*c1.x + r[1][1]*c1.y + r[2][1]*c1.z + v[1]);
  z=c2.z-(r[0][2]*c1.x + r[1][2]*c1.y + r[2][2]*c1.z + v[2]);
  t=x*x + y*y + z*z;
  return t;
}

/*-----------------------------------------------------------
/
/   rms_window - calculate rms value for a given residue window
/
/------------------------------------------------------------*/
void rms_window(long atoms[2][MAXRES+1], long l_atoms[2][MAXRES+1],
                long out[2][MAXRES+1], pdb_struct pdb[2],
                part_struct *part, long k1, long k2)
{
  long i, j, j1, j2, j12, n1, n2, al1, al2;
  long n_max, rw, rwi, rwj, rwn, rwn0, rwn1, rwn2;
  char al;
  float r[3][3], v[3];
  lga_workspace *ws = part->ws;
  if (!ws) { /* defensive */
    if (lga_ws_init(&part->ws) != 0) {
      fprintf(stderr, "\n# ERROR! Failed to allocate LGA workspace (rms_window)\n");
      exit(1);
    }
    ws = part->ws;
  }
  long (*calphas)[MAXRES+1] = ws ? (long (*)[MAXRES+1])ws->calphas2 : NULL;
  long *rwa1 = ws ? ws->rwa1 : NULL;
  long *rwa2 = ws ? ws->rwa2 : NULL;

  n_max=atoms[0][0]+atoms[1][0];
  n1=1;
  n2=1;
  i=1;
  j=1;
  rw=0;
  rwn=2*part->rw_l+1;
  for(j12=1;j12<=n_max;j12++) {
    part->rw_rms[j12]=-1.0;
    al1=0;
    al2=0;
    if(out[0][i]==atoms[0][n1] && n1<=atoms[0][0]) al1=1;
    if(out[1][i]==atoms[1][n2] && n2<=atoms[1][0]) al2=1;

    al='-';
    if(al1==1 && al2==1 && i<=k2) {
      if(j<=k1 && out[0][i]==l_atoms[0][j] && out[1][i]==l_atoms[1][j]) {
        j++;
        al='+';
      }
      else {
        al='#';
      }

      j1=out[0][i];
      j2=out[1][i];

      if(al=='+') {
        rw++;
        rwa1[rw]=j1;
        rwa2[rw]=j2;
        if(rw>=rwn) {
          rwn0=rw-part->rw_l;
          rwn1=rw-rwn+1;
          rwn2=rw;
          for(rwi=rwn1;rwi<=rwn2;rwi++) {
            rwj=rwi-rwn1+1;
            calphas[0][rwj]=atoms[0][rwa1[rwi]];
            calphas[1][rwj]=atoms[1][rwa2[rwi]];
          }
          calphas[0][0]=rwn;
          calphas[1][0]=rwn;

          part->rw_rms[rwn0]=rmsd(rwn, calphas, pdb, r, v);

        }
      }
      i++;
      n1++;
      n2++;
    }
    else if(al1==1 && n2<=atoms[1][0]) {
      n2++;
    }
    else if(al2==1 && n1<=atoms[0][0]) {
      n1++;
    }
    else {
      n1++;
      n2++;
    }
  }

  return;
}

/*-----------------------------------------------------------
/
/   gdc_dist - calculates distance between atoms
/
/------------------------------------------------------------*/
float gdc_dist(long j1, long j2,
               atom_coords coords1[MAXATOMS+1], atom_coords coords2[MAXATOMS+1],
               float r[3][3], float v[3])
{  
  float t, xa[3], xb[3];

        xa[0]=coords1[j1].R.x;
        xa[1]=coords1[j1].R.y;
        xa[2]=coords1[j1].R.z;

        xb[0]=r[0][0]*(coords2[j2].R.x - v[0])+
              r[0][1]*(coords2[j2].R.y - v[1])+
              r[0][2]*(coords2[j2].R.z - v[2]);
        xb[1]=r[1][0]*(coords2[j2].R.x - v[0])+
              r[1][1]*(coords2[j2].R.y - v[1])+
              r[1][2]*(coords2[j2].R.z - v[2]);
        xb[2]=r[2][0]*(coords2[j2].R.x - v[0])+
              r[2][1]*(coords2[j2].R.y - v[1])+
              r[2][2]*(coords2[j2].R.z - v[2]);

        t=(xa[0]-xb[0])*(xa[0]-xb[0])+
          (xa[1]-xb[1])*(xa[1]-xb[1])+
          (xa[2]-xb[2])*(xa[2]-xb[2]);
        t=sqrt(t);

  return t;
}

/*-----------------------------------------------------------
/
/   distance_calc - calculates distance between CA atoms
/
/------------------------------------------------------------*/
float distance_calc(long atoms[2][MAXRES+1],
               pdb_struct pdb[2], float r[3][3], float v[3],
               long j1, long j2)
{
  float t, xa[3], xb[3];

        xa[0]=pdb[0].atom[atoms[0][j1]].R.x;
        xa[1]=pdb[0].atom[atoms[0][j1]].R.y;
        xa[2]=pdb[0].atom[atoms[0][j1]].R.z;

        xb[0]=r[0][0]*(pdb[1].atom[atoms[1][j2]].R.x - v[0])+
              r[0][1]*(pdb[1].atom[atoms[1][j2]].R.y - v[1])+
              r[0][2]*(pdb[1].atom[atoms[1][j2]].R.z - v[2]);
        xb[1]=r[1][0]*(pdb[1].atom[atoms[1][j2]].R.x - v[0])+
              r[1][1]*(pdb[1].atom[atoms[1][j2]].R.y - v[1])+
              r[1][2]*(pdb[1].atom[atoms[1][j2]].R.z - v[2]);
        xb[2]=r[2][0]*(pdb[1].atom[atoms[1][j2]].R.x - v[0])+
              r[2][1]*(pdb[1].atom[atoms[1][j2]].R.y - v[1])+
              r[2][2]*(pdb[1].atom[atoms[1][j2]].R.z - v[2]);

        t=(xa[0]-xb[0])*(xa[0]-xb[0])+
          (xa[1]-xb[1])*(xa[1]-xb[1])+
          (xa[2]-xb[2])*(xa[2]-xb[2]);
        t=sqrt(t);

  return t;
}

/*-----------------------------------------------------------
/
/   best_fit - search for the best fit of the two structures
/
/------------------------------------------------------------*/
void best_fit(long atoms[2][MAXRES+1], pdb_struct pdb[2],
              float r[3][3], float v[3], part_struct *part, pdata_struct *pdata)
{
  long i, j, k, n, n1, n2, ok, ok1, ok11, ok2, ok21, nfits, nfits_max, number;
  float dist_min, dist_one, rms_local, rms_max, t, t_max;
  float l_r1[3][3], l_v1[3], l_r2[3][3], l_v2[3];
  lga_workspace *ws = part->ws;
  if (!ws) { /* defensive */
    if (lga_ws_init(&part->ws) != 0) {
      fprintf(stderr, "\n# ERROR! Failed to allocate LGA workspace (best_fit)\n");
      exit(1);
    }
    ws = part->ws;
  }
  long (*sup)[MAXRES+1]  = ws ? (long (*)[MAXRES+1])ws->sup  : NULL;
  long (*sup1)[MAXRES+1] = ws ? (long (*)[MAXRES+1])ws->sup1 : NULL;
  long (*sup2)[MAXRES+1] = ws ? (long (*)[MAXRES+1])ws->sup2 : NULL;
  char fname_lga[200], mname[200];

  strcpy(fname_lga,part->fname_lga);

  n=0;
  nfits=0;
  nfits_max=1;
  dist_min=0.45;  //  0.45
  dist_one=0.90;  //  0.90
  rms_max=999.99;
  if(part->sia==0) {
    number=atoms[0][0]-atoms[1][0]-part->fit_g;
    n=atoms[1][0];
    j=0;
    while(j<=number) {
      k=1+j;
      for(i=1;i<=n;i++) {
        sup1[0][i]=atoms[0][k];
        sup1[1][i]=atoms[1][i];
        if(part->fit_n==i) k=k+part->fit_g;
        k++;
      }
      sup1[0][0]=n;
      sup1[1][0]=n;
      rms_local=rmsd(n, sup1, pdb, l_r1, l_v1);
      if((rms_local<999 && nfits<20) || rms_local<=rms_max) {
        ok=1;
        ok1=0;
        ok2=0;
        if(part->fit_b>=part->fit_n) ok11=1;
        else ok11=0;
        if(part->fit_b>=n-part->fit_n) ok21=1;
        else ok21=0;
        t_max=-999.999;
        for(i=1;i<=n;i++) {
          t=distance_calc(sup1, pdb, l_r1, l_v1, i, i);
          if(t>t_max) t_max=t;
          if(t>part->dist_cutoff) {
            ok=0;
          }
          else {
            if(t<=dist_min) {
              if(part->fit_n>=i) ok1++;
              else ok2++;
            }
            if(t<=dist_one) {
              if(part->fit_n>=i) ok11++;
              else ok21++;
            }
          }
        }
        if(ok==1 && ok1>=part->fit_b && ok2>=part->fit_b &&
           ok11>part->fit_b && ok21>part->fit_b) {
          if(nfits==0 || rms_max<=rms_local) rms_max=rms_local;
          if(nfits<20) {
            nfits++;
            nfits_max=nfits;
          }
          for(i=0;i<3;i++) {
            part->opt_r[nfits_max][0][i]=l_r1[0][i];
            part->opt_r[nfits_max][1][i]=l_r1[1][i];
            part->opt_r[nfits_max][2][i]=l_r1[2][i];
            part->opt_v[nfits_max][i]=l_v1[i];
          }
          for(i=0;i<=n;i++) {
            part->opt_sup[nfits_max][0][i]=sup1[0][i];
            part->opt_sup[nfits_max][1][i]=sup1[1][i];
          }
          part->opt_sup_rms[nfits_max]=rms_local;
          part->opt_all_rms[nfits_max]=t_max;
          rms_max=0.0;
          for(i=1;i<=nfits;i++) {
            if(rms_max<=part->opt_sup_rms[i]) {
              nfits_max=i;
              rms_max=part->opt_sup_rms[nfits_max];
            }
          }
        }
      }
      j++;
    }
  }
  if(part->sia==1) {
    number=atoms[1][0]-atoms[0][0]-part->fit_g;
    n=atoms[0][0];
    j=0;
    while(j<=number) {
      k=1+j;
      ok=0;
      n1=0;
      n2=0;
      for(i=1;i<=n;i++) {
        if(ok==0) {
          n1++;
          sup1[0][n1]=atoms[0][i];
          sup1[1][n1]=atoms[1][k];
        }
        if(ok==1) {
          n2++;
          sup2[0][n2]=atoms[0][i];
          sup2[1][n2]=atoms[1][k];
        }
        if(part->fit_n==i) {
          ok=1;
          k=k+part->fit_g;
        }
        k++;
      }
      if((part->fit_r[0]!='#' && (n1<3 || n2<3)) || n!=n1+n2) {
        printf(" ERROR! CHECK the number of residues in Molecule1:\n");
        printf("        Molecule1_f1  %4ld \n",n1);
        printf("        Molecule1_f2  %4ld \n",n2);
        printf("        number of residues < 3  ERROR! \n");
        exit(0);
      }
      sup1[0][0]=n1;
      sup1[1][0]=n1;
      t=0.0;
      if(n1>=3) t=rmsd(n1, sup1, pdb, l_r1, l_v1);
      rms_local=t*t*n1;
      sup2[0][0]=n2;
      sup2[1][0]=n2;
      t=0.0;
      if(n2>=3) t=rmsd(n2, sup2, pdb, l_r2, l_v2);
      rms_local=rms_local+t*t*n2;
      rms_local=sqrt(rms_local/n);
      if((rms_local<999 && nfits<20) || rms_local<=rms_max) {
        ok=1;
        ok1=0;
        ok2=0;
        if(part->fit_b>=n1) ok11=1;
        else ok11=0;
        if(part->fit_b>=n2) ok21=1;
        else ok21=0;
        t_max=-999.999;
        for(i=1;i<=n1;i++) {
          t=distance_calc(sup1, pdb, l_r1, l_v1, i, i);
          if(t>t_max) t_max=t;
          if(t>part->dist_cutoff) {
            ok=0;
          }
          else {
            if(t<=dist_min) {
              ok1++;
            }
            if(t<=dist_one) {
              ok11++;
            }
          }
        }
        if(ok==1 && ok1>=part->fit_b && ok11>part->fit_b) {
          for(i=1;i<=n2;i++) {
            t=distance_calc(sup2, pdb, l_r2, l_v2, i, i);
            if(t>t_max) t_max=t;
            if(t>part->dist_cutoff) {
              ok=0;
            }
            else {
              if(t<=dist_min) {
                ok2++;
              }
              if(t<=dist_one) {
                ok21++;
              }
            }
          }
          if(ok==1 && ok2>=part->fit_b && ok21>part->fit_b) {
            if(nfits==0 || rms_max<=rms_local) rms_max=rms_local;
            if(nfits<20) {
              nfits++;
              nfits_max=nfits;
            }
            for(i=0;i<3;i++) {
              part->opt_r[nfits_max][0][i]=l_r1[0][i];
              part->opt_r[nfits_max][1][i]=l_r1[1][i];
              part->opt_r[nfits_max][2][i]=l_r1[2][i];
              part->opt_v[nfits_max][i]=l_v1[i];
              part->opt_r1[nfits_max][0][i]=l_r1[0][i];
              part->opt_r1[nfits_max][1][i]=l_r1[1][i];
              part->opt_r1[nfits_max][2][i]=l_r1[2][i];
              part->opt_v1[nfits_max][i]=l_v1[i];
              part->opt_r2[nfits_max][0][i]=l_r2[0][i];
              part->opt_r2[nfits_max][1][i]=l_r2[1][i];
              part->opt_r2[nfits_max][2][i]=l_r2[2][i];
              part->opt_v2[nfits_max][i]=l_v2[i];
            }
            part->opt_sup[nfits_max][0][0]=n;
            part->opt_sup[nfits_max][1][0]=n;
            for(i=1;i<=n1;i++) {
              part->opt_sup[nfits_max][0][i]=sup1[0][i];
              part->opt_sup[nfits_max][1][i]=sup1[1][i];
            }
            for(i=1;i<=n2;i++) {
              part->opt_sup[nfits_max][0][n1+i]=sup2[0][i];
              part->opt_sup[nfits_max][1][n1+i]=sup2[1][i];
            }
            part->opt_sup_rms[nfits_max]=rms_local;
            part->opt_all_rms[nfits_max]=t_max;
            rms_max=0.0;
            for(i=1;i<=nfits;i++) {
              if(rms_max<=part->opt_sup_rms[i]) {
                nfits_max=i;
                rms_max=part->opt_sup_rms[nfits_max];
              }
            }
          }
        }
      }
      j++;
    }
  }
  if(nfits>0) {
    for(k=1;k<=nfits;k++) {
      strcpy(part->fname_lga,fname_lga);
      strcat(part->fname_lga,part->model_index[k].aa);
      strcpy(mname,pdata->mname);
      strcat(mname,part->model_index[k].aa);
      part->opt_sup_rms[0]=part->opt_sup_rms[k];
      part->opt_all_rms[0]=part->opt_all_rms[k];
      for(i=0;i<3;i++) {
        r[0][i]=part->opt_r[k][0][i];
        r[1][i]=part->opt_r[k][1][i];
        r[2][i]=part->opt_r[k][2][i];
        v[i]=part->opt_v[k][i];
      }
      for(i=0;i<=n;i++) {
        sup[0][i]=part->opt_sup[k][0][i];
        sup[1][i]=part->opt_sup[k][1][i];
      }
      if(part->sia==1) {
        for(i=0;i<3;i++) {
          part->opt_r1[0][0][i]=part->opt_r1[k][0][i];
          part->opt_r1[0][1][i]=part->opt_r1[k][1][i];
          part->opt_r1[0][2][i]=part->opt_r1[k][2][i];
          part->opt_v1[0][i]=part->opt_v1[k][i];
          part->opt_r2[0][0][i]=part->opt_r2[k][0][i];
          part->opt_r2[0][1][i]=part->opt_r2[k][1][i];
          part->opt_r2[0][2][i]=part->opt_r2[k][2][i];
          part->opt_v2[0][i]=part->opt_v2[k][i];
        }
      }
      printf("\nFIT-MODEL: %3ld %3ld    %s \n",k,nfits,mname);
      part->fit_si=print_lga(atoms, sup, sup, pdb, r, v, part, n, n);
      write_output(part, pdata, atoms, r, v, part->fname_lga);
    }
    strcpy(part->fname_lga,fname_lga);
  }
  else {
    write_output(part, pdata, atoms, r, v, part->fname_lga);
  }

  return;
}
