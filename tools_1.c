#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "tools.h"

/*-----------------------------------------------------------
/
/   clean PDB structure
/
/------------------------------------------------------------*/
void clean_pdb(pdb_struct *molecule)
{
  long  i;
  char junk[10];

  strcpy(junk," \0");

  molecule->n_atoms_all=0;
  molecule->n_atoms=0;
  molecule->n_aa_all=0;
  molecule->n_aa=0;
  molecule->n_chains=0;
  for(i=0;i<200;i++) {
    molecule->chain_name[i]=' ';
  }
  for(i=0;i<=MAXRES;i++) {
    molecule->nb_ca[i]=0;
    molecule->atom[i].serial=0;
    strcpy(molecule->atom[i].name,junk);
    strcpy(molecule->atom[i].alt_loc,junk);
    strcpy(molecule->atom[i].res_name,junk);
    strcpy(molecule->atom[i].chain_id,junk);
    molecule->atom[i].res_seq=0;
    strcpy(molecule->atom[i].i_code,junk);
    molecule->atom[i].R.x=0.0; 
    molecule->atom[i].R.y=0.0; 
    molecule->atom[i].R.z=0.0;
    molecule->atom[i].occupancy=0.0;
    molecule->atom[i].temp_factor=0.0;
    strcpy(molecule->atom[i].seg_id,junk);
    strcpy(molecule->atom[i].element,junk);
    strcpy(molecule->atom[i].charge,junk);

    strcpy(molecule->atom[i].res_i,junk);
    molecule->atom[i].res_n_local=0;
  }
  return;
}

/*-----------------------------------------------------------
/
/   clean pdata structure
/
/------------------------------------------------------------*/
void clean_pdata(pdata_struct *pdata)
{
  long i;
  
  pdata->use_CB=0;
  pdata->cb_calc=0;
  pdata->cb_v=0.0;
  pdata->m_v=0.0;
  pdata->o_v=0.0;
  pdata->ignore_errors=0;
  pdata->one_output=1;
  pdata->n_atoms=0;
  pdata->s_analysis=1;
  pdata->sda=0;
  pdata->s_distance=5.0;
  pdata->ca_atoms=1;
  pdata->len_atom=3;
  pdata->ah_i=0;
  pdata->aa1_ch='*';
  pdata->aa2_ch='*';
  pdata->parameters_n=0;
  for(i=0;i<50;i++) {
    strcpy(pdata->parameters[i].line,"#");
  }
  strcpy(pdata->atoms,"CA\0");
  strcpy(pdata->mol1,"*");
  strcpy(pdata->mol2,"*");

  strcpy(pdata->molecule[0].molecule_name,"Molecule_1");
  strcpy(pdata->molecule[1].molecule_name,"Molecule_2");

  clean_pdb(&pdata->molecule[0]);
  clean_pdb(&pdata->molecule[1]);
  return;
}

/*-----------------------------------------------------------
/
/   clean part structure
/
/------------------------------------------------------------*/
void clean_part(part_struct *part)
{
  long i, j, n;
  char ch[10], ind[10];

  n=0;
  strncpy(ch,"0123456789",10);
  for(i=0;i<10;i++) {
    for(j=0;j<10;j++) {
      n++;
      ind[0]='_';
      ind[1]=ch[i];
      ind[2]=ch[j];
      ind[3]='\0';
      sscanf(ind,"%3s",part->model_index[n].aa);
    }
  }
  part->accuracy_opt=1;
  part->accuracy_gdt_n=56;
  part->accuracy_gdt_step=0.1;
  part->accuracy_lga_n=8;
  part->accuracy_lga_step=0.7;
  part->gdc=0;
  part->gdc_bin=20;
  part->gdc_ref=0;
  part->gdc_sup.ok=-1;
  part->gdc_sup.n=0;
  part->gdc_set.ok=-1;
  part->gdc_set.n=0;
  part->gdc_at.ok=-1;
  part->gdc_at_mcb=-1;
  part->gdc_at.n=0;
  part->gdc_eval.ok=-1;
  part->gdc_eval.n=0;
  for(j=0;j<2;j++) {
    part->error[j]=0;
    part->aar_g1[j]=9999;
    part->aar_g2[j]=-9999;
    part->aar_n1[j]=-9999;
    part->aar_n2[j]=9999;
    part->er[j].ok=-1;
    part->er[j].n=0;
  }
  for(i=0;i<8;i++) {
    strcpy(part->resrangaa[i].aa,"#");
    part->resrangaa[i].found=0;
  }
  part->all_rmsd=0;
  part->eval=0;
  part->swap=0;
  part->gdt=0;
  part->resranges=0;
  part->summary.rms_local=999.99;
  part->summary.seq_id=0.0;
  part->summary.lcs_gdt_ts=0.0;
  part->summary.n=0;
  part->fit_g=0;
  part->fit_b=0;
  part->fit_n=0;
  part->rw_l=0;
  part->lN_n=1;
  part->sia=0;
  strcpy(part->fit_r,"#");
  part->lcs_gdt_print=1;
  part->m_n_aa=0;
  part->t_n_aa=0;
  part->full_print=0;
  part->check=0;
  part->stral=0;
  part->stral_r=0.5;
  part->isp=0;
  part->isp_iter=0;
  part->isp_iter_start=0;
  part->isp_dist=999.999;
  part->lga_m=0;
  part->lga_w=0.75;
  part->rmsd=0.0;
  part->rmsd_isp=0.0;
  part->isp_cutoff[0]=50.0;
  part->isp_cutoff[1]=25.0;
  part->isp_cutoff[2]=15.0;
  part->isp_cutoff[3]=10.0;
  part->isp_cutoff[4]=7.00;
  part->isp_cutoff[5]=5.00;
  part->isp_cutoff[6]=3.00;
  part->isp_cutoff[7]=2.00;
  part->isp_cutoff[8]=1.50;
  part->isp_cutoff[9]=1.25;
  part->isp_cutoff[10]=1.1;
  for(i=0;i<3;i++) {
    part->opt_r[0][0][i]=0.0;
    part->opt_r[0][1][i]=0.0;
    part->opt_r[0][2][i]=0.0;
    part->opt_v[0][i]=0.0;
  }
  for(i=0;i<11;i++) {
    part->rmsd_align[i]=999.999;
    part->ear[i].align[0]=-1;
    part->best_ind[i]=i;
  }
  for(i=0;i<=20;i++) {
    part->gdt_cutoff[i]=0.5*i;
    part->opt_sup1[i][0]=0;
    part->opt_sup2[i][0]=0;
    part->opt_sup_rms[i]=999.999;
    part->opt_all_rms[i]=999.999;
    for(j=0;j<=MAXRES;j++) {
      part->opt_set_nb[i][j]=0;
    }
  }
  for(j=0;j<=MAXRES;j++) {
    part->opt_lcs[j]=0;
  }

  n=0;
  for(i=1;i<=MAXGAP;i++) {
    for(j=1;j<i;j++) {
      n++;
      part->pairs[0][n]=i;
      part->pairs[1][n]=j; 
      n++;
      part->pairs[0][n]=j; 
      part->pairs[1][n]=i;
    } 
    n++;
    part->pairs[0][n]=i; 
    part->pairs[1][n]=i;
  }
  part->pairs[0][0]=n;
  part->pairs[1][0]=n;

  return;
}

/*-----------------------------------------------------------
/
/   read_aamol2 - read AAMOL* and LGA records
/
/------------------------------------------------------------*/
void read_aamol2(mol2_aa aa_mol2[2], FILE* fp)
{
  long  i, j, n_mol1, n_mol2;
  char keyword[200], line[200], sub[10], sub1[10], sub2[10], sub3[10], sub4[10];

  n_mol1=0;
  n_mol2=0;
  while(fgets(line,200,fp)!=NULL) {
    for(i=0;i<200;i++) keyword[i]=0;
    sscanf(line,"%s",keyword);
    if(!strncmp(keyword,"AAMOL1\0",7)) {
      n_mol1++;
//
        for(j=0;j<10;j++) sub[j]=0;
        strncpy(sub,&line[21],1);
        sscanf(sub,"%s",sub1);          // chain
        for(j=0;j<10;j++) sub[j]=0;
        strncpy(sub,&line[22],5);
        sscanf(sub,"%s",sub2);          // residue
        for(j=0;j<10;j++) sub[j]=0;
        if(sub1[0]!=' ') {
          sub[0]='_';
          sub[1]=sub1[0];
          sub[2]='\0';
        }
        strcat(sub2,sub);               // residue_chain
        for(j=0;j<10;j++) sub[j]=0;
        strncpy(sub,&line[70],1);
        sscanf(sub,"%s",sub3);          // aa
//
      sscanf(sub2,"%s",aa_mol2[0].num[n_mol1].aa);
//      printf("# TEST! MOL1: %ld Residue: %6s \n",n_mol1,aa_mol2[0].num[n_mol1].aa);
    }
    if(!strncmp(keyword,"AAMOL2\0",7)) {
      n_mol2++;
//
        for(j=0;j<10;j++) sub[j]=0;
        strncpy(sub,&line[21],1);
        sscanf(sub,"%s",sub1);          // chain
        for(j=0;j<10;j++) sub[j]=0;
        strncpy(sub,&line[22],5);
        sscanf(sub,"%s",sub2);          // residue
        for(j=0;j<10;j++) sub[j]=0;
        if(sub1[0]!=' ') {
          sub[0]='_';
          sub[1]=sub1[0];
          sub[2]='\0';
        }
        strcat(sub2,sub);               // residue_chain
        for(j=0;j<10;j++) sub[j]=0;
        strncpy(sub,&line[70],1);
        sscanf(sub,"%s",sub3);          // aa
//
      sscanf(sub2,"%s",aa_mol2[1].num[n_mol2].aa);
    }
    if(!strncmp(keyword,"LGA\0",4)) {
      sscanf(line,"%s %s %s %s %s %s ",keyword,sub1,sub2,sub3,sub4,sub);
      if(sub[0]!='-') {
        n_mol1++;
        n_mol2++;
        strcpy(aa_mol2[0].num[n_mol1].aa,sub2);
        strcpy(aa_mol2[1].num[n_mol2].aa,sub4);
      }
    }
  }
  aa_mol2[0].n_aa=n_mol1;
  aa_mol2[1].n_aa=n_mol2;

  return;
}

/*-----------------------------------------------------------
/
/   read_pdb_2sets - read PDB format file (2 molecules;
/                    sets of coordinates)
/
/------------------------------------------------------------*/
void read_pdb_2sets(mol2_aa aa_mol2[2],
                    pdata_struct *pdata, part_struct *part, FILE* fp)
{
  long  i, j, k, n_aa, n_atoms, n_atoms_all, n_chains, nm, endflag, hetatm;
  long  n_aamol1, n_aamol2, alt_loc_ok, ok, selected, completed;
  long  in_C, in_CA, in_N, in_O, in_CB, cb_calculated, bug;
  char keyword[200], line[200], sub[10], sub_res[10], atom_name[10];
  char alt_loc, ch, code1;
  char chcurr[10], chprev[10], curr_res[10], last_res[10], prev_res[10], bmo_res[10];

  /* Read in the molecule
  -------------------------------------------*/
  nm=-1;
  endflag=0;
  selected=0;
  completed=0;
  n_atoms_all=0;
  n_atoms=0;
  n_aa=0;
  n_chains=0;
  ch='#';
  alt_loc='#';
  alt_loc_ok=0;
  strcpy(last_res,"###\0");
  strcpy(prev_res,"###\0");
  strcpy(bmo_res,"###\0");
  strcpy(chprev,"#\0");
  ok=1;
  in_C=0;
  in_CA=0;
  in_N=0;
  in_O=0;
  in_CB=0;
  cb_calculated=0;
  hetatm=0;
  n_aamol1=aa_mol2[0].n_aa;
  n_aamol2=aa_mol2[1].n_aa;
  while(fgets(line,200,fp)!=NULL && endflag<2) {
    for(i=0;i<200;i++) keyword[i]=0;
    sscanf(line,"%s",keyword);
    if(!strncmp(keyword,"MOLECULE\0",9))  {
      nm++;
      if(nm>1) {
        printf("\n# ERROR! Check the number of MOLECULE records.");
        printf("\n#        Each molecule should begin from MOLECULE record");
        printf("\n#        and end with END record \n");
        fclose(fp);
        exit(0);
      }
      endflag=0;
      n_atoms_all=0;
      n_atoms=0;
      n_aa=0;
      n_chains=0;
      ch='#';
      alt_loc='#';
      alt_loc_ok=0;
      strcpy(last_res,"###\0");
      strcpy(prev_res,"###\0");
      ok=1;
      in_C=0;
      in_CA=0;
      in_N=0;
      in_O=0;
      in_CB=0;
      cb_calculated=1;
      completed=1;
      i=0;
      while(line[i]!='\n' && i<200) i++;
      if(i>=10) {
        sscanf(line,"%s %s",keyword,pdata->molinput[nm].molecule_name);
      }
    }
    else if(!strncmp(keyword,"END",3))    {
      if(nm<0) {
        printf("\n# ERROR! There is no MOLECULE record before END record \n");
        fclose(fp);
        exit(0);
      }
      endflag++;
      if(completed!=1) {
        if(part->full_print==1 || part->check==1) {
          printf("\n# WARNING! (IN) No %s atom inside the residue: %s (%s - %s molecule: %ld)\n",
                  pdata->atoms,last_res,prev_res,curr_res,nm+1);
        }
        if(part->er[nm].ok>=0) {
          for(j=0;j<part->er[nm].n;j++) {
            if(part->er[nm].b[j].found==0) {
              if(strlen(last_res) == strlen(part->er[nm].b[j].aa) &&
                 !strcmp(last_res,part->er[nm].b[j].aa)) {
                if(part->full_print==1 || part->check==1) {
                  printf("\n# WARNING! Selected residue has been replaced %s (%s - %s molecule: %ld)\n",
                          part->er[nm].b[j].aa,curr_res,bmo_res,nm+1);
                }
                strcpy(part->er[nm].b[j].aa,curr_res);
              }
            }
            if(part->er[nm].e[j].found==0) {
              if(strlen(last_res) == strlen(part->er[nm].e[j].aa) &&
                 !strcmp(last_res,part->er[nm].e[j].aa)) {
                if(part->full_print==1 || part->check==1) {
                  printf("\n# WARNING! Selected residue has been replaced %s (%s - %s molecule: %ld)\n",
                          part->er[nm].e[j].aa,prev_res,bmo_res,nm+1);
                }
                if(pdata->cb_calc==1) {
                  strcpy(part->er[nm].e[j].aa,bmo_res);
                }
                else {
                  strcpy(part->er[nm].e[j].aa,prev_res);
                }
              }
            }
          }
        }
      }
    }
    else if((!strncmp(keyword,"ATOM\0",5) && pdata->ah_i<2) ||
            (!strncmp(keyword,"HETATM",6) && pdata->ah_i!=1)) {
      if(nm<0) {
        printf("\n# ERROR! Each molecule should begin from MOLECULE record");
        printf("\n#        and end with END record \n");
        fclose(fp);
        exit(0);
      }
      i=0;
      selected=0;
      for(j=0;j<10;j++) sub_res[j]=0;
      strncpy(sub_res,&line[17],3);

      while(line[i]!='\n' && i<120) i++;
      if(i>=54 && strncmp(sub_res,"HOH",3)) {
        n_atoms_all++;
        for(j=0;j<10;j++) {
          sub[j]=0;
          atom_name[j]=0;
        }
        strncpy(sub,&line[12],4);
        sscanf(sub,"%s",atom_name);

        if(pdata->cb_calc==1) {
          if(!strncmp(atom_name,"N\0",2) || 
             !strncmp(atom_name,"C\0",2) || 
             !strncmp(atom_name,"CA\0",3) ||
             !strncmp(atom_name,"O\0",2) ||
             !strncmp(atom_name,"CB\0",3)) {
            if(pdata->use_CB == 1 || strncmp(atom_name,"CB\0",3)) {
              selected=1;
            }
          }
        }
        else if(!strncmp(atom_name,pdata->atoms,pdata->len_atom)) {
//          strcpy(atom_name,"CA\0");
          selected=1;
        }
        strcpy(chcurr," \0");
        for(j=0;j<10;j++) sub[j]=0;
        strncpy(sub,&line[21],1);
        sscanf(sub,"%s",chcurr);
        if(nm==0 && pdata->aa1_ch!='*') {
          if(pdata->aa1_ch!=chcurr[0]) selected=0;
        }
        if(nm==1 && pdata->aa2_ch!='*') {
          if(pdata->aa2_ch!=chcurr[0]) selected=0;
        }
        
        for(j=0;j<10;j++) sub[j]=0;
        strncpy(sub,&line[22],5);
        sscanf(sub,"%s",curr_res);
        for(j=0;j<10;j++) sub[j]=0;

        if(chcurr[0]!=' ') {
          sub[0]='_';
          sub[1]=chcurr[0];
          sub[2]='\0';
        }
        strcat(curr_res,sub);
        if(strcmp(curr_res,last_res) && 
           ((pdata->cb_calc==0) || (selected!=0 && pdata->cb_calc==1))) {
          if(completed!=1) {
            if(part->full_print==1 || part->check==1) {
              printf("\n# WARNING! (IN) No %s atom inside the residue: %s (%s - %s molecule: %ld)\n",
                      pdata->atoms,last_res,prev_res,curr_res,nm+1);
            }
            if(part->er[nm].ok>=0) {
              for(j=0;j<part->er[nm].n;j++) {
                if(part->er[nm].b[j].found==0) {
                  if(strlen(last_res) == strlen(part->er[nm].b[j].aa) &&
                     !strcmp(last_res,part->er[nm].b[j].aa)) {
                    if(part->full_print==1 || part->check==1) {
                      printf("\n# WARNING! Selected residue has been replaced %s (%s - %s molecule: %ld)\n",
                                part->er[nm].b[j].aa,curr_res,last_res,nm+1);
                    }
                    strcpy(part->er[nm].b[j].aa,curr_res);
                  }
                }
                if(part->er[nm].e[j].found==0) {
                  if(strlen(last_res) == strlen(part->er[nm].e[j].aa) &&
                     !strcmp(last_res,part->er[nm].e[j].aa)) {
                    if(part->full_print==1 || part->check==1) {
                      printf("\n# WARNING! Selected residue has been replaced %s (%s - %s molecule: %ld)\n",
                                part->er[nm].e[j].aa,prev_res,last_res,nm+1);
                    }
                    strcpy(part->er[nm].e[j].aa,prev_res);
                  }
                }
              }
            }
            if(completed>1) {
              n_atoms=n_atoms-5-in_CB+completed;
              n_aa--;
            }
          }
          else if(part->er[nm].ok>=0) {
            for(j=0;j<part->er[nm].n;j++) {
              if(strlen(last_res) == strlen(part->er[nm].b[j].aa) &&
                 !strcmp(last_res,part->er[nm].b[j].aa)) {
                part->er[nm].b[j].found=-1;
              }
              if(strlen(last_res) == strlen(part->er[nm].e[j].aa) &&
                 !strcmp(last_res,part->er[nm].e[j].aa)) {
                part->er[nm].e[j].found=-1;
              }
              else if(chcurr[0]!=' ') {
                if(strlen(part->er[nm].b[j].aa)==1 && strlen(part->er[nm].e[j].aa)==1 &&
                   part->er[nm].b[j].aa[0] == part->er[nm].e[j].aa[0] &&
                   part->er[nm].b[j].aa[0] == chcurr[0]) {
                  part->er[nm].b[j].found=-1;
                  part->er[nm].e[j].found=-1;
                }
              }
            }
          }
          completed=0;
          strcpy(last_res,curr_res);
          in_C=0;
          in_CA=0;
          in_N=0;
          in_O=0;
          in_CB=0;
          cb_calculated=0;
        }
      }
      if(selected==1) {
        n_atoms++;
        if(n_atoms>=MAXRES) {
          printf("\n# ERROR! Check the number of ATOM records. (Limit MAX = %d) \n\n",
                  MAXRES);
          fclose(fp);
          exit(0);
        }

        for(j=0;j<10;j++) sub[j]=0;
        strncpy(sub,&line[6],5);
        sscanf(sub,"%ld",&pdata->molinput[nm].atom[n_atoms].serial);

        strcpy(pdata->molinput[nm].atom[n_atoms].name,atom_name);

// Modified residue name: http://ligand-expo.rcsb.org/index.html

        if(!strncmp(keyword,"HETATM",6)) {
          if(pdata->ca_atoms==0) {
            strcpy(sub_res,"XXX");
            hetatm=0;
          }
          else if(!strncmp(sub_res,"ABA",3) || // modified ALA
             !strncmp(sub_res,"BFD",3) ||      // modified ASP
             !strncmp(sub_res,"CAS",3) ||      // modified CYS
             !strncmp(sub_res,"CCS",3) ||      // modified CYS
             !strncmp(sub_res,"CME",3) ||      // modified CYS
             !strncmp(sub_res,"CSD",3) ||      // modified CYS
             !strncmp(sub_res,"CSE",3) ||      // modified CYS
             !strncmp(sub_res,"CSO",3) ||      // modified CYS
             !strncmp(sub_res,"CSS",3) ||      // modified CYS
             !strncmp(sub_res,"CSU",3) ||      // modified CYS
             !strncmp(sub_res,"CSW",3) ||      // modified CYS
             !strncmp(sub_res,"CSX",3) ||      // modified CYS
             !strncmp(sub_res,"CXM",3) ||      // modified MET
             !strncmp(sub_res,"CY3",3) ||      // modified CYS
             !strncmp(sub_res,"DAR",3) ||      // modified ARG
             !strncmp(sub_res,"DGL",3) ||      // modified GLU
             !strncmp(sub_res,"DHA",3) ||      // modified ALA dehydroalanine
             !strncmp(sub_res,"DIL",3) ||      // modified ILE
             !strncmp(sub_res,"DTR",3) ||      // modified TRP
             !strncmp(sub_res,"DXX",3) ||
             !strncmp(sub_res,"ETA",3) ||
             !strncmp(sub_res,"FME",3) ||      // modified MET
             !strncmp(sub_res,"GLZ",3) ||      // modified GLY
             !strncmp(sub_res,"HYP",3) ||      // modified PRO Hydroxyproline
             !strncmp(sub_res,"KCX",3) ||      // modified LYS
             !strncmp(sub_res,"LLP",3) ||      // modified LYS
             !strncmp(sub_res,"MGG",3) ||      // modified ARG
             !strncmp(sub_res,"MLY",3) ||      // modified LYS
             !strncmp(sub_res,"MSE",3) ||      // modified MET Selenomethionine
             !strncmp(sub_res,"M3L",3) ||      // modified LYS
             !strncmp(sub_res,"OCS",3) ||      // modified CYS
             !strncmp(sub_res,"PCA",3) ||      // modified GLU
             !strncmp(sub_res,"PTR",3) ||      // modified TYR
             !strncmp(sub_res,"PYT",3) ||
             !strncmp(sub_res,"ROP",3) ||
             !strncmp(sub_res,"SAC",3) ||      // modified SER
             !strncmp(sub_res,"SEC",3) ||      // modified CYS Selenocysteine, or modified ALA
             !strncmp(sub_res,"SEP",3) ||      // modified SER
             !strncmp(sub_res,"TPO",3) ||      // modified THR
             !strncmp(sub_res,"TPQ",3) ||      // modified TYR
             !strncmp(sub_res,"YCM",3)) {      // modified CYS
            strcpy(sub_res,"XXX");
            hetatm=0;
            ok=1;
          }
          else {
            hetatm=1;
          }
        }
        if(!strncmp(keyword,"ATOM\0",5)) {

          C_321(sub_res,&code1);
          if(pdata->ca_atoms==0 && code1=='#') {
            strcpy(sub_res,"XXX");
          }
          else if(!strncmp(sub_res,"ACE",3) || // Acetyl beginning group
             !strncmp(sub_res,"ADE",3) ||
             !strncmp(sub_res,"ASH",3) ||      // Neutral ASP
             !strncmp(sub_res,"ASX",3) ||      // ASN (N) or ASP (D)
             !strncmp(sub_res,"CYM",3) ||      // modified CYS Negative CYS
             !strncmp(sub_res,"CYT",3) ||
             !strncmp(sub_res,"CYX",3) ||      // Cystine, S--S crosslink, SS-bonded CYS
             !strncmp(sub_res,"GLH",3) ||      // modified GLN or Neutral GLU
             !strncmp(sub_res,"GLP",3) ||      // pyroglutamic acid
             !strncmp(sub_res,"GLX",3) ||      // GLN (Q) or GLU (E)
             !strncmp(sub_res,"GUA",3) ||
             !strncmp(sub_res,"HID",3) ||      // Histidine, delta H, neutral HIS, proton HD1 present
             !strncmp(sub_res,"HIE",3) ||      // Histidine, epsilon H, neutral HIS, proton HE2 present
             !strncmp(sub_res,"HIP",3) ||      // modified HIS Histidine, protonated, positive HIS
             !strncmp(sub_res,"HSD",3) ||
             !strncmp(sub_res,"HSE",3) ||      // modified SER homoserine
             !strncmp(sub_res,"HSL",3) ||      // modified SER homoserine lactone
             !strncmp(sub_res,"HSP",3) ||
             !strncmp(sub_res,"HYL",3) ||      // Hydroxylysine
             !strncmp(sub_res,"HYP",3) ||      // modified PRO Hydroxyproline
             !strncmp(sub_res,"IVA",3) ||      // Isovaline
             !strncmp(sub_res,"LYN",3) ||      // modified LYS Neutral LYS
             !strncmp(sub_res,"MSE",3) ||      // modified MET Selenomethionine
             !strncmp(sub_res,"NHE",3) ||      // Amine ending group
             !strncmp(sub_res,"NLE",3) ||      // modified LEU Norleucine
             !strncmp(sub_res,"NME",3) ||      // N-methylamine ending group 
             !strncmp(sub_res,"SEC",3) ||      // modified CYS Selenocysteine, or modified ALA
             !strncmp(sub_res,"THY",3) ||
             !strncmp(sub_res,"TYM",3) ||      // Negative TYR
             !strncmp(sub_res,"XXX",3) ||
             !strncmp(sub_res,"UNK",3)) {
            strcpy(sub_res,"XXX");
          }
          if(hetatm==1) ok=1;
          hetatm=0;
        }
        sscanf(sub_res,"%s",pdata->molinput[nm].atom[n_atoms].res_name);

        strcpy(pdata->molinput[nm].atom[n_atoms].chain_id,chcurr);
        if(ch!=pdata->molinput[nm].atom[n_atoms].chain_id[0]) {
          ch=pdata->molinput[nm].atom[n_atoms].chain_id[0];
          pdata->molinput[nm].chain_name[n_chains]=ch;
          if(hetatm==0) n_chains++;
          alt_loc='#';
          alt_loc_ok=0;
          strcpy(prev_res,"###\0");
        }

        for(j=0;j<10;j++) sub[j]=0;
        strncpy(sub,&line[22],4);
        sscanf(sub,"%ld",&pdata->molinput[nm].atom[n_atoms].res_seq);

        for(j=0;j<10;j++) sub[j]=0;
        strncpy(sub,&line[26],1);
        sscanf(sub,"%s",pdata->molinput[nm].atom[n_atoms].i_code);

        sscanf(curr_res,"%s",pdata->molinput[nm].atom[n_atoms].res_i);
//      printf("\n# TEST: %ld  %ld %s ",n_atoms,pdata->molinput[nm].atom[n_atoms].res_seq,pdata->molinput[nm].atom[n_atoms].res_i);
        if(nm==0) {
          if(hetatm==0 && (n_aamol1>0 || pdata->sda==2)) {
            hetatm=1;
            for(j=1;j<=n_aamol1;j++) {
              if(!strcmp(pdata->molinput[nm].atom[n_atoms].res_i,aa_mol2[0].num[j].aa))
                hetatm=0;
            }
          }
        }
        if(nm==1) {
          if(hetatm==0 && (n_aamol2>0 || pdata->sda==2)) {
            hetatm=1;
            for(j=1;j<=n_aamol2;j++) {
              if(!strcmp(pdata->molinput[nm].atom[n_atoms].res_i,aa_mol2[1].num[j].aa))
                hetatm=0;
            }
          }
        }
        if(strcmp(pdata->molinput[nm].atom[n_atoms].res_i,prev_res)) {
          n_aa++;
          alt_loc='#';
          alt_loc_ok=0;
          if(n_aa>=MAXRES) {
            printf("\n# ERROR! Check the number of residues. (Limit MAX = %d) \n\n",
                    MAXRES);
            fclose(fp);
            exit(0);
          }
          if(ok==0 && hetatm==0 && pdata->ca_atoms==1) {
            if(pdata->ignore_errors==0) {
              printf("\n# ERROR! There is no %s atom inside the residue: %s , %s (molecule: %ld)\n\n",
                      pdata->atoms,prev_res,pdata->molinput[nm].atom[n_atoms].res_i,nm+1);
              fclose(fp);
              exit(0);
            }
            else {
              n_aa--;
            }
            completed=1;
          }
          if(hetatm!=0) {
            n_aa--;
            completed=1;
          }
          ok=0;
          strcpy(prev_res,pdata->molinput[nm].atom[n_atoms].res_i);
        }

        for(j=0;j<10;j++) sub[j]=0;
        strncpy(sub,&line[16],1);
        sscanf(sub,"%s",pdata->molinput[nm].atom[n_atoms].alt_loc);
        if(alt_loc_ok==1 && line[16]==' ') {
          alt_loc=' ';
          alt_loc_ok=0;
        }
        if(alt_loc_ok==1 && line[16]!=alt_loc) {
          alt_loc_ok=2;
        }
        if(alt_loc_ok==0 && line[16]!=alt_loc) {
          alt_loc_ok=1;
          alt_loc=line[16];
        }

        for(j=0;j<10;j++) sub[j]=0;
        strncpy(sub,&line[30],8);
        bug=buger(sub,line,30,37);
        if(bug>1) {
          printf("\n# ERROR! Check the coordinates (X) of the ATOM number %ld (molecule: %ld)\n\n",
                  pdata->molinput[nm].atom[n_atoms].serial,nm+1);
          fclose(fp);
          exit(0);
        }
        sscanf(sub,"%f",&pdata->molinput[nm].atom[n_atoms].R.x);

        for(j=0;j<10;j++) sub[j]=0;
        strncpy(sub,&line[38],8);
        bug=buger(sub,line,38,45);
        if(bug>1) {
          printf("\n# ERROR! Check the coordinates (Y) of the ATOM number %ld (molecule: %ld)\n\n",
                  pdata->molinput[nm].atom[n_atoms].serial,nm+1);
          fclose(fp);
          exit(0);
        }
        sscanf(sub,"%f",&pdata->molinput[nm].atom[n_atoms].R.y);

        for(j=0;j<10;j++) sub[j]=0;
        strncpy(sub,&line[46],8);
        bug=buger(sub,line,46,53);
        if(bug>1) {
          printf("\n# ERROR! Check the coordinates (Z) of the ATOM number %ld (molecule: %ld)\n\n",
                  pdata->molinput[nm].atom[n_atoms].serial,nm+1);
          fclose(fp);
          exit(0);
        }
        sscanf(sub,"%f",&pdata->molinput[nm].atom[n_atoms].R.z);

        if(alt_loc_ok!=2) {
// Check for duplicated atoms. BEGIN
          k=strlen(pdata->molinput[nm].atom[n_atoms].res_i);
          for(j=1;j<n_atoms;j++) {
            if(k==strlen(pdata->molinput[nm].atom[j].res_i) && pdata->cb_calc!=1 &&
               !strcmp(pdata->molinput[nm].atom[n_atoms].res_i,pdata->molinput[nm].atom[j].res_i)) {
//              printf("# TEST: Resname %s length %ld\n",sub_res,(long)strlen(sub_res)); 
              if(pdata->ignore_errors<2 && sub_res[0] != ' ') { // e.g. hetatm ' CA' are ignored  
                printf("\n# ERROR! Check molecule: %ld and the amino acid number %s and name %s", 
                        nm+1,pdata->molinput[nm].atom[n_atoms].res_i,sub_res);
                printf("\n#        This number or some atom positions within this amino acid are duplicated.\n");
                if(pdata->ignore_errors==0) {
                  if(hetatm==1) {
                    printf("#        Only one atom (first) is used for calculations %d. \n",strncmp(sub_res,"CA\0",3));
                  }
                  else {
                    fclose(fp);
                    exit(0);
                  }
                }
              }
            }
            if(fabs(pdata->molinput[nm].atom[j].R.x-pdata->molinput[nm].atom[n_atoms].R.x)<0.01 &&
               fabs(pdata->molinput[nm].atom[j].R.y-pdata->molinput[nm].atom[n_atoms].R.y)<0.01 &&
               fabs(pdata->molinput[nm].atom[j].R.z-pdata->molinput[nm].atom[n_atoms].R.z)<0.01) {
              if(completed != 1 || strncmp(atom_name,"CB\0",3)) {
                if(pdata->ignore_errors<2) {
                  printf("\n# ERROR! Check the ATOM numbers %ld and %ld (molecule: %ld)",
                          pdata->molinput[nm].atom[j].serial,
                          pdata->molinput[nm].atom[n_atoms].serial,nm+1);
                  printf("\n#        The distance between atoms is too small! ");
                  printf("\n#        Check if some ATOM coordinates are not repeated.\n");
                }
                alt_loc_ok=2;
                if(pdata->ignore_errors==0) {
                  fclose(fp);
                  exit(0);
                }
              }
            }
          }
// Check for duplicated atoms. END
        }
        if(i>54) {
          for(j=0;j<10;j++) sub[j]=0;
          strncpy(sub,&line[54],6);
          sscanf(sub,"%f",&pdata->molinput[nm].atom[n_atoms].occupancy);
        }
        if(i>60) {
          for(j=0;j<10;j++) sub[j]=0;
          strncpy(sub,&line[60],6);
          sscanf(sub,"%f",&pdata->molinput[nm].atom[n_atoms].temp_factor);
        }
        if(i>72) {
          for(j=0;j<10;j++) sub[j]=0;
          strncpy(sub,&line[72],4);
          sscanf(sub,"%s",pdata->molinput[nm].atom[n_atoms].seg_id);
        }
        if(i>76) {
          for(j=0;j<10;j++) sub[j]=0;
          strncpy(sub,&line[76],2);
          sscanf(sub,"%s",pdata->molinput[nm].atom[n_atoms].element);
        }
        if(i>78) {
          for(j=0;j<10;j++) sub[j]=0;
          strncpy(sub,&line[78],2);
          sscanf(sub,"%s",pdata->molinput[nm].atom[n_atoms].charge);
        }
      }
      if(i<30) {
        printf("\n# ERROR! Check ATOM record lines below the ATOM number %ld (molecule: %ld)\n",
                pdata->molinput[nm].atom[n_atoms].serial,nm+1);
        fclose(fp);
        exit(0);
      }
      if(i<54) {
        printf("\n# ERROR! Check ATOM records and the format of XYZ coordinates.");
        printf("\n#        Check lines below the ATOM number %ld (molecule: %ld)\n",
                pdata->molinput[nm].atom[n_atoms].serial,nm+1);
        fclose(fp);
        exit(0);
      }
      if(selected==1) {
//        printf("\n# TEST! MOL: %ld Residue: %6s ATOM Name: %4s  Coord: X = %f",nm+1,pdata->molinput[nm].atom[n_atoms].res_i,pdata->molinput[nm].atom[n_atoms].name,pdata->molinput[nm].atom[n_atoms].R.x);
        if(hetatm==1) {
          n_atoms--;
        }
        else if(alt_loc_ok==2 && pdata->cb_calc==0) {
          alt_loc_ok=1;
// THERE is a problem with CA alternate positions when -bmo option is used!!!
//          printf("\n# TEST: Alternate ATOM records %ld (molecule: %ld)\n",
//                  pdata->molinput[nm].atom[n_atoms].serial,nm+1);
          n_atoms--;
        }
        else {
          if(pdata->cb_calc!=1 || completed!=1 || strncmp(atom_name,"CB\0",3)) {
          if(!strncmp(pdata->molinput[nm].atom[n_atoms].name,pdata->atoms,pdata->len_atom) ||
             completed==1) {
            pdata->molinput[nm].nb_ca[n_aa]=n_atoms;
            pdata->molinput[nm].atom[n_atoms].res_n_local=n_aa;
            if(ok==1 || completed==1) {
              if(alt_loc_ok!=2 || line[16]==' ' || pdata->cb_calc!=1) {
                printf("\n# WARNING! There is %s atom duplicated inside the residue %s (molecule: %ld)",
                      pdata->molinput[nm].atom[n_atoms].name,prev_res,nm+1);
                printf("\n#          Only one atom (first) is used for calculations. \n");
              }
              n_atoms--;
            }
            else {
              strcpy(chprev,chcurr);
            }
            ok=1;
          }
          }
          if(pdata->cb_calc==1 && (completed!=1 || !strncmp(atom_name,"CB\0",3))) {
            if(completed==0) {
              completed=5;
              strcpy(pdata->molinput[nm].atom[n_atoms].name,"BMO\0");
            }
            j=0;
            if(!strncmp(atom_name,"CB\0",3)) {
              if(in_CB==0) j=MAXRES-4;
              in_CB=1;
            }
            else if(!strncmp(atom_name,"O\0",2)) {
              if(in_O==0) j=MAXRES-3;
              in_O=1;
            }
            else if(!strncmp(atom_name,"C\0",2)) {
              if(in_C==0) j=MAXRES-2;
              in_C=1;
            }
            else if(!strncmp(atom_name,"CA\0",3)) {
              if(in_CA==0) j=MAXRES-1;
              in_CA=1;
            }
            else if(!strncmp(atom_name,"N\0",2)) {
              if(in_N==0) j=MAXRES;
              in_N=1;
            }
            else {
              printf("\n# ERROR! Check ATOM records %s (molecule: %ld)\n",curr_res,nm+1);
              fclose(fp);
              exit(0);
            }
            if(j>0) {
              if(strncmp(atom_name,"CB\0",3)) completed--;
              pdata->molinput[1].atom[j].R.x=pdata->molinput[nm].atom[n_atoms].R.x;
              pdata->molinput[1].atom[j].R.y=pdata->molinput[nm].atom[n_atoms].R.y;
              pdata->molinput[1].atom[j].R.z=pdata->molinput[nm].atom[n_atoms].R.z;
              if(completed==1) {
/*
                printf("\n# TEST! Completed: ATOM N  = %f",pdata->molinput[1].atom[MAXRES].R.x);
                printf("\n# TEST! Completed: ATOM CA = %f",pdata->molinput[1].atom[MAXRES-1].R.x);
                printf("\n# TEST! Completed: ATOM C  = %f",pdata->molinput[1].atom[MAXRES-2].R.x);
                printf("\n# TEST! Completed: ATOM O  = %f",pdata->molinput[1].atom[MAXRES-3].R.x);
                printf("\n# TEST! Completed: ATOM CB = %f",pdata->molinput[1].atom[MAXRES-4].R.x);
*/
                calc_cb(pdata->cb_v, pdata->m_v, pdata->o_v, pdata->molinput, in_CB);
                if(cb_calculated==0) {
                  cb_calculated=1;
                  n_atoms=n_atoms-3-in_CB;
                }
                else {
                  n_atoms--;
                }
                j=MAXRES-1;
                pdata->molinput[nm].atom[n_atoms].R.x=pdata->molinput[0].atom[j].R.x;
                pdata->molinput[nm].atom[n_atoms].R.y=pdata->molinput[0].atom[j].R.y;
                pdata->molinput[nm].atom[n_atoms].R.z=pdata->molinput[0].atom[j].R.z;
                strcpy(bmo_res,curr_res);
//                printf("\n# TEST! Completed: ATOM X  = %f\n",pdata->molinput[nm].atom[n_atoms].R.x);
              }
            }
            else {
              if(alt_loc_ok!=2 || line[16]==' ') {
                printf("\n# WARNING! There is %s atom duplicated inside the residue %s (molecule: %ld)",
                      pdata->molinput[nm].atom[n_atoms].name,prev_res,nm+1);
                printf("\n#          Only first one atom is used for CB calculations. \n");
              }
              n_atoms--;
            }
          }
          else {
            completed=1;
          }
        }
      }
    }      
    else if(!strncmp(keyword,"HET\0",4))    {}
    else {}
    if(endflag==1) {
      if(completed!=1) {
//        printf("\n# WARNING! (END) There is no %s atom inside the residue: %s , %s : %ld , %ld\n",
//                  pdata->atoms,last_res,curr_res,n_atoms,n_aa);
        if(completed>1) {
          n_atoms=n_atoms-5-in_CB+completed;
          n_aa--;
        }
      }
      if(ok==0 && hetatm==0 && pdata->ca_atoms==1 && pdata->ignore_errors==0) {
        printf("\n# ERROR! (IN) There is no %s atom inside the residue %s (molecule: %ld)\n",
                pdata->atoms,prev_res,nm+1);
        fclose(fp);
        exit(0);
      }
      pdata->molinput[nm].n_atoms_all=n_atoms_all;
      pdata->molinput[nm].n_atoms=n_atoms;
      pdata->molinput[nm].n_aa=n_aa;
      pdata->molinput[nm].n_chains=n_chains;
      if(nm<1) {endflag=0;}
      if(n_atoms<3 || n_aa<3) {
        printf("\n# ERROR! Check the number of selected atoms in molecule %ld",nm+1);
        printf("\n#        Selected number of residues: %ld",n_aa);
        printf("\n#        Use if necessary the options to specify chains: -ch1:A , -ch2:B ,...");
        printf("\n#        Check residue numbering and chain ID (see selected parameters)");
        printf("\n#        ");
        for(j=0;j<pdata->parameters_n;j++) {
          printf("%s ",pdata->parameters[j].line);
        }
        printf("\n");
        fclose(fp);
        exit(0);
      }
    }
  }
  if(nm!=1) {
    printf("\n# ERROR! Check the number of molecules begining from MOLECULE record \n");
    fclose(fp);
    exit(0);
  }
  if(endflag==0 && hetatm==0) {
    if(ok==0 && pdata->ca_atoms==1 && pdata->ignore_errors==0) {
      printf("\n# ERROR! (END) There is no %s atom inside the residue %s (molecule: %ld)\n",
              pdata->atoms,prev_res,nm+1);
      fclose(fp);
      exit(0);
    }
    pdata->molinput[nm].n_atoms_all=n_atoms_all;
    pdata->molinput[nm].n_atoms=n_atoms;
    pdata->molinput[nm].n_aa=n_aa;
    pdata->molinput[nm].n_chains=n_chains;
    if(n_aa<3) {
      printf("\n# ERROR! Check the number of residues in molecule %ld",nm+1);
      printf("\n#        Check residue numbering and chain ID (see selected parameters)");
      printf("\n#        ");
      for(j=0;j<pdata->parameters_n;j++) {
        printf("%s ",pdata->parameters[j].line);
      }
      printf("\n");
      fclose(fp);
      exit(0);
    }
  }
  if(nm>0 && endflag==0) {
    printf("\n# WARNING! Check the END record in molecule 2 \n");
  }
  
  return;
}

/*-----------------------------------------------------------
/
/   C_321 - converts an amino acid code3 to code1
/
/------------------------------------------------------------*/
void C_321(char* code3, char* code1)
{
  *code1='#';
  if(!strncmp(code3, "ALA", 3) || !strncmp(code3, "ala", 3)) *code1='A';
  else if(!strncmp(code3, "VAL", 3) || !strncmp(code3, "val", 3)) *code1='V';
  else if(!strncmp(code3, "LEU", 3) || !strncmp(code3, "leu", 3)) *code1='L';
  else if(!strncmp(code3, "ILE", 3) || !strncmp(code3, "ile", 3)) *code1='I';
  else if(!strncmp(code3, "PRO", 3) || !strncmp(code3, "pro", 3)) *code1='P';
  else if(!strncmp(code3, "MET", 3) || !strncmp(code3, "met", 3)) *code1='M';
  else if(!strncmp(code3, "PHE", 3) || !strncmp(code3, "phe", 3)) *code1='F';
  else if(!strncmp(code3, "TRP", 3) || !strncmp(code3, "trp", 3)) *code1='W';
  else if(!strncmp(code3, "GLY", 3) || !strncmp(code3, "gly", 3)) *code1='G';
  else if(!strncmp(code3, "SER", 3) || !strncmp(code3, "ser", 3)) *code1='S';
  else if(!strncmp(code3, "THR", 3) || !strncmp(code3, "thr", 3)) *code1='T';
  else if(!strncmp(code3, "CYS", 3) || !strncmp(code3, "cys", 3)) *code1='C';
  else if(!strncmp(code3, "TYR", 3) || !strncmp(code3, "tyr", 3)) *code1='Y';
  else if(!strncmp(code3, "ASN", 3) || !strncmp(code3, "asn", 3)) *code1='N';
  else if(!strncmp(code3, "GLN", 3) || !strncmp(code3, "gln", 3)) *code1='Q';
  else if(!strncmp(code3, "ASP", 3) || !strncmp(code3, "asp", 3)) *code1='D';
  else if(!strncmp(code3, "GLU", 3) || !strncmp(code3, "glu", 3)) *code1='E';
  else if(!strncmp(code3, "LYS", 3) || !strncmp(code3, "lys", 3)) *code1='K';
  else if(!strncmp(code3, "ARG", 3) || !strncmp(code3, "arg", 3)) *code1='R';
  else if(!strncmp(code3, "HIS", 3) || !strncmp(code3, "his", 3)) *code1='H';
  else if(!strncmp(code3, "UNK", 3) || !strncmp(code3, "unk", 3)) *code1='X';
  else if(!strncmp(code3, "XXX", 3) || !strncmp(code3, "xxx", 3)) *code1='X';

  return;
}

/*-----------------------------------------------------------
/
/   copy_atom - copy atom2 data to atom1
/
/------------------------------------------------------------*/
void copy_atom(atom_struct *atom1, atom_struct *atom2)
{
  atom1->serial=atom2->serial;
  strcpy(atom1->name,atom2->name);
  strcpy(atom1->alt_loc,atom2->alt_loc);
  strcpy(atom1->res_name,atom2->res_name);
  strcpy(atom1->chain_id,atom2->chain_id);
  atom1->res_seq=atom2->res_seq;
  strcpy(atom1->i_code,atom2->i_code);
  atom1->R.x=atom2->R.x;
  atom1->R.y=atom2->R.y;
  atom1->R.z=atom2->R.z;
  atom1->occupancy=atom2->occupancy;
  atom1->temp_factor=atom2->temp_factor;
  strcpy(atom1->seg_id,atom2->seg_id);
  strcpy(atom1->element,atom2->element);
  strcpy(atom1->charge,atom2->charge);
  strcpy(atom1->res_i,atom2->res_i);
  atom1->res_n_local=atom2->res_n_local;
  
  return;
}

/*-----------------------------------------------------------
/
/   select_atoms - select atoms for processing
/
/------------------------------------------------------------*/
long select_atoms(long ca_atoms[2][MAXRES+1],
                 pdata_struct *pdata, 
                 part_struct *part)
{
  long i, j, k, m, n, ok, fit_rn;
  char chcurr[10], chprev[10];

  fit_rn=0;

  for(i=0;i<2;i++) {
//    printf("\n# TEST_1 %ld %ld %ld\n",part->er[i].ok,i,pdata->sda);
    strcpy(chcurr,"#\0");
    strcpy(chprev,"#\0");
    k=0;
    for(j=1;j<=pdata->molinput[i].n_aa;j++) {
      ok=1;
      if(chcurr[0]!=pdata->molinput[i].atom[j].chain_id[0]) {
        strcpy(chprev,chcurr);
        strcpy(chcurr,pdata->molinput[i].atom[j].chain_id);
      }
//  BEGIN aa and gap selection
      if(part->aar_n1[i]>pdata->molinput[i].atom[j].res_seq ||
         part->aar_n2[i]<pdata->molinput[i].atom[j].res_seq) {
        ok=0;
      }
      if(part->aar_g1[i]<=pdata->molinput[i].atom[j].res_seq &&
         part->aar_g2[i]>=pdata->molinput[i].atom[j].res_seq) {
        ok=0;
      }
//  END aa and gap selection
//  BEGIN er selection
      if(part->er[i].ok>=0 && pdata->sda!=2) {
        for(n=0;n<part->er[i].n;n++) {
          if(part->er[i].b[n].aa[0]!='#') {
//  er: begin residue section 
            if(strlen(pdata->molinput[i].atom[j].res_i) ==
               strlen(part->er[i].b[n].aa) &&
               !strcmp(pdata->molinput[i].atom[j].res_i,part->er[i].b[n].aa)) {
              if(part->er[i].e[n].found!=1) part->er[i].ok=1;
              for(m=0;m<part->er[i].n;m++) {
                if(part->er[i].b[m].aa[0]=='+') part->er[i].b[m].aa[0]='#';
              }
              part->er[i].b[n].aa[0]='+';
              part->er[i].b[n].found=1;
            }
//  er: begin chain section
            else if(strlen(part->er[i].b[n].aa) == 1 && 
               strlen(part->er[i].e[n].aa) == 1 &&
               chcurr[0] != ' ' &&
               (part->er[i].b[n].aa[0] == part->er[i].e[n].aa[0] ||
               part->er[i].b[n].aa[0]=='+') && part->er[i].e[n].aa[0]==chcurr[0]) {
              part->er[i].ok=1;
              part->er[i].b[n].found=1;
              part->er[i].e[n].found=1;
              for(m=0;m<part->er[i].n;m++) {
                if(part->er[i].b[m].aa[0]=='+') part->er[i].b[m].aa[0]='#';
              }
              part->er[i].b[n].aa[0]='+';
            }
          }
        }
        if(part->er[i].ok==0) {
          ok=0;
        }
        for(n=0;n<part->er[i].n;n++) {
          if(part->er[i].b[n].aa[0]!='#') {
//  er: end residue section
            if(strlen(pdata->molinput[i].atom[j].res_i) ==
               strlen(part->er[i].e[n].aa) &&
               !strcmp(pdata->molinput[i].atom[j].res_i,part->er[i].e[n].aa)) {
              if(part->er[i].b[n].found==1) part->er[i].b[n].aa[0]='#';
              part->er[i].e[n].found=1;
              part->er[i].ok=0;
            }
//  er: end chain section
            else if(strlen(part->er[i].b[n].aa) == 1 && 
               strlen(part->er[i].e[n].aa) == 1 &&
               chcurr[0] != ' ' &&
               part->er[i].b[n].aa[0]=='+' && part->er[i].e[n].aa[0]!=chcurr[0] &&
               (chcurr[0]!=' ' || chcurr[0]!=chprev[0])) {
              ok=0;
              part->er[i].ok=0;
              part->er[i].b[n].aa[0]='#';
            }
          }
        }
      }
//  END er selection
      if(ok==1) {
        k++;
        copy_atom(&pdata->molecule[i].atom[k],&pdata->molinput[i].atom[j]);
        pdata->molecule[i].atom[k].res_n_local=k;
        pdata->molecule[i].nb_ca[k]=k;
        ca_atoms[i][k]=pdata->molecule[i].nb_ca[k];
        if(pdata->s_analysis==5 && i==1-part->sia && fit_rn==0) {
          if(!strcmp(pdata->molinput[i].atom[j].res_i,part->fit_r)) {
            fit_rn=k;
          }
        }
      }
    }
    ca_atoms[i][0]=k;
    pdata->molecule[i].n_aa=k;
    pdata->molecule[i].n_atoms=pdata->molinput[i].n_atoms;
    pdata->molecule[i].n_atoms_all=pdata->molinput[i].n_atoms_all;
    strcpy(pdata->molecule[i].molecule_name,pdata->molinput[i].molecule_name);

   if(pdata->sda!=2) {
    for(n=0;n<part->er[i].n;n++) {
      if(part->er[i].b[n].found!=1) {
        printf("\n# ERROR! Selected ATOM from residue %s not found",part->er[i].b[n].aa);
        printf("\n#        Check parameters -er%ld (check the description, NO overlaps allowed)",i+1);
        printf("\n#        ");
        for(j=0;j<pdata->parameters_n;j++) {
          printf("%s ",pdata->parameters[j].line);
        }
        printf("\n");
        exit(0);
      }
      if(part->er[i].e[n].found!=1) {
        printf("\n# ERROR! Selected ATOM from residue %s not found",part->er[i].e[n].aa);
        printf("\n#        Check parameters -er%ld (check the description, NO overlaps allowed)",i+1);
        printf("\n#        ");
        for(j=0;j<pdata->parameters_n;j++) {
          printf("%s ",pdata->parameters[j].line);
        }
        printf("\n");
        exit(0);
      }
    }
   }
  }
  if(ca_atoms[0][0]<ca_atoms[1][0])
    pdata->n_atoms=ca_atoms[0][0];
  else
    pdata->n_atoms=ca_atoms[1][0];
  part->fit_n=fit_rn;

  if(pdata->n_atoms<3) {
    printf("\n# ERROR! Check the number of selected atoms: %ld",pdata->n_atoms);
    printf("\n#        Check residue numbering and chain ID (see selected parameters)\n");
  }
  
  return pdata->n_atoms;
}

/*-----------------------------------------------------------
/
/   sda_list - list of the sequence identical residues
/
/------------------------------------------------------------*/
void sda_list(mol2_aa aa_mol2[2], long ca_atoms[2][MAXRES+1],
             pdata_struct *pdata)
{
  long i, j, k, n_aamol, ok;
  char chain_id[2][MAXRES+1];
  
  n_aamol=0;
  for(j=1;j<=pdata->molecule[1].n_aa;j++) {
    for(i=1;i<=pdata->molecule[0].n_aa;i++) {
      ok=0;
      if(pdata->aa1_ch != '*' && pdata->aa2_ch != '*') {
        if(pdata->aa1_ch == pdata->molecule[0].atom[ca_atoms[0][i]].chain_id[0] &&
           pdata->aa2_ch == pdata->molecule[1].atom[ca_atoms[1][j]].chain_id[0]) {
          if(pdata->molecule[0].atom[ca_atoms[0][i]].res_seq ==
             pdata->molecule[1].atom[ca_atoms[1][j]].res_seq &&
             pdata->molecule[0].atom[ca_atoms[0][i]].i_code[0] ==
             pdata->molecule[1].atom[ca_atoms[1][j]].i_code[0]) {
            ok=1;
          }
        }
      }
      if(ok==1 || !strcmp(pdata->molecule[0].atom[ca_atoms[0][i]].res_i,
                          pdata->molecule[1].atom[ca_atoms[1][j]].res_i)) {
        n_aamol++;
        strcpy(aa_mol2[0].num[n_aamol].aa,pdata->molecule[0].atom[ca_atoms[0][i]].res_i);
        strcpy(aa_mol2[1].num[n_aamol].aa,pdata->molecule[1].atom[ca_atoms[1][j]].res_i);
        chain_id[0][n_aamol]=pdata->molecule[0].atom[ca_atoms[0][i]].chain_id[0];
        chain_id[1][n_aamol]=pdata->molecule[1].atom[ca_atoms[1][j]].chain_id[0];
        if(pdata->ignore_errors==2) {
          for(k=1;k<n_aamol;k++) {
            if(!strcmp(aa_mol2[0].num[k].aa,pdata->molecule[0].atom[ca_atoms[0][i]].res_i) &&
               chain_id[0][k]==pdata->molecule[0].atom[ca_atoms[0][i]].chain_id[0]) {
              printf("\n# ERROR! There is a duplicated residue %s in the first molecule \n",aa_mol2[0].num[k].aa);
              exit(0);
            }
            if(!strcmp(aa_mol2[1].num[k].aa,pdata->molecule[1].atom[ca_atoms[1][j]].res_i) &&
               chain_id[1][k]==pdata->molecule[1].atom[ca_atoms[1][j]].chain_id[0]) {
              printf("\n# ERROR! There is a duplicated residue %s in the second molecule \n",aa_mol2[1].num[k].aa);
              exit(0);
            }
          }
        }
      }
    }
  }
  
  aa_mol2[0].n_aa=n_aamol;
  aa_mol2[1].n_aa=n_aamol;

  return;
}

/*-----------------------------------------------------------
/
/   print_aamol2 - list of the residues from the molecules 1 and 2
/
/------------------------------------------------------------*/
void print_aamol2(long ca_atoms[2][MAXRES+1], pdata_struct *pdata)
{
  long i,j,k;
  char ch;

  printf("\n");
  for(i=0;i<=1;i++) {
    printf("\n# List of residues and atom coordinates for LGA processing from Molecule%1ld: %s ",
              i+1,pdata->molecule[i].molecule_name);
    k=0;
    printf("\n\n>%s  FASTA sequence of residues from Molecule%1ld\n",
              pdata->molecule[i].molecule_name,i+1);
    for(j=1;j<=pdata->molecule[i].n_aa;j++) {
      k=k+1;
      if(k>100) {
        k=1;
        printf("\n");
      }
      C_321(pdata->molecule[i].atom[ca_atoms[i][j]].res_name,&ch);
      printf("%c",ch);
    }
    printf("\n");
    for(j=1;j<=pdata->molecule[i].n_aa;j++) {
      C_321(pdata->molecule[i].atom[ca_atoms[i][j]].res_name,&ch);
      printf("\nAAMOL%1ld %4ld %4s %3s %1s%4ld%1s",
                i+1,j,
                pdata->molecule[i].atom[ca_atoms[i][j]].name,
                pdata->molecule[i].atom[ca_atoms[i][j]].res_name,
                pdata->molecule[i].atom[ca_atoms[i][j]].chain_id,
                pdata->molecule[i].atom[ca_atoms[i][j]].res_seq,
                pdata->molecule[i].atom[ca_atoms[i][j]].i_code);
      printf("   %8.3f%8.3f%8.3f%6.2f%6.2f   %c",
                pdata->molecule[i].atom[ca_atoms[i][j]].R.x, 
                pdata->molecule[i].atom[ca_atoms[i][j]].R.y, 
                pdata->molecule[i].atom[ca_atoms[i][j]].R.z, 
                pdata->molecule[i].atom[ca_atoms[i][j]].occupancy, 
                pdata->molecule[i].atom[ca_atoms[i][j]].temp_factor, 
                ch); 
    }
    printf("\n");
  }

  return;
}

/*-----------------------------------------------------------
/
/   euler_angles - computing Euler angles from a rotation matrix
/
/------------------------------------------------------------*/
void euler_angles(float r[3][3], float euler[3][4])
{
  long i;
  double c1, s1, c2, del, pi, th1, th2, ps1, ps2, ph1, ph2, err;
  pi=(double)3.1415926535897932384626433;
  err=(double)0.00001;

// XYZ convention
//   phi is about x-axis
//   theta is about y-axis
//   psi is about z-axis
  if(fabs(r[0][2]) < 1.0-err) {
    th1=-asin(r[0][2]);
    th2=pi - th1;
    c1=cos(th1);
    c2=cos(th2);
    ps1=atan2(r[1][2]/c1,r[2][2]/c1);
    ps2=atan2(r[1][2]/c2,r[2][2]/c2);
    ph1=atan2(r[0][1]/c1,r[0][0]/c1);
    ph2=atan2(r[0][1]/c2,r[0][0]/c2);
  }
  else {
    ph1=(double)0.0;
    ph2=(double)pi*0.5;
    del=atan2(r[1][0],r[2][0]);
    if(r[0][2] < 0.0) {
      th1=ph2;
      th2=ph2;
      ps1=ph1 + del;
      ps2=ph2 + del;
    }
    else {
      th1=-ph2;
      th2=-ph2;
      ps1=-ph1 + del;
      ps2=-ph2 + del;
    }
  }
  if(th1 < -pi) th1=2.0*pi + th1;
  if(th2 > pi) th2=th2 - 2.0*pi;
  if(ps1 < -pi) ps1=2.0*pi + ps1;
  if(ps2 > pi) ps2=ps2 - 2.0*pi;
  if(ph1 < -pi) ph1=2.0*pi + ph1;
  if(ph2 > pi) ph2=ph2 - 2.0*pi;

  euler[0][0]=th1; euler[0][2]=(double)180.0*th1/pi;
  euler[0][1]=th2; euler[0][3]=(double)180.0*th2/pi;
  euler[1][0]=ps1; euler[1][2]=(double)180.0*ps1/pi;
  euler[1][1]=ps2; euler[1][3]=(double)180.0*ps2/pi;
  euler[2][0]=ph1; euler[2][2]=(double)180.0*ph1/pi;
  euler[2][1]=ph2; euler[2][3]=(double)180.0*ph2/pi;
  
  if(fabs(th1)<fabs(th2)) i=0;
  else i=1;
  printf("\nEuler angles from the ROTATION matrix. Conventions XYZ and ZXZ:");
  printf("\n           Phi     Theta       Psi   [DEG:       Phi     Theta       Psi ]");
  printf("\nXYZ: %9.6f %9.6f %9.6f   [DEG: %9.4f %9.4f %9.4f ]",
         euler[2][i],euler[0][i],euler[1][i],euler[2][2],euler[0][2],euler[1][2]);

// TEST XYZ:
/*
  long n1, n2;
  float m[3][3];
  double c3, s2, s3, z;
    c1 = cos(euler[2][i]); // Phi
    s1 = sin(euler[2][i]);
    c2 = cos(euler[0][i]); // Theta
    s2 = sin(euler[0][i]);
    c3 = cos(euler[1][i]); // Psi
    s3 = sin(euler[1][i]);
    m[0][0] =  c1 * c2;
    m[1][0] =  c1 * s2 * s3 - s1 * c3;
    m[2][0] =  c1 * s2 * c3 + s1 * s3;
    m[0][1] =  s1 * c2;
    m[1][1] =  s1 * s2 * s3 + c1 * c3;
    m[2][1] =  s1 * s2 * c3 - c1 * s3;
    m[0][2] = -s2;
    m[1][2] =  c2 * s3;
    m[2][2] =  c2 * c3;
  z=0.0;
  for(n1=0;n1<3;n1++)
    for(n2=0;n2<3;n2++)
      z=z+fabs(m[n1][n2]-r[n1][n2]);
  printf("\nTEST XYZ: error = %10.6f",z);
   
    printf("\nTEST XYZ: Unitary ROTATION matrix calculated from euler angles:");
    printf("\n  %10.6f   %10.6f   %10.6f ",
              m[0][0],m[1][0],m[2][0]);
    printf("\n  %10.6f   %10.6f   %10.6f ",
              m[0][1],m[1][1],m[2][1]);
    printf("\n  %10.6f   %10.6f   %10.6f \n",
              m[0][2],m[1][2],m[2][2]);
*/

// ZXZ convention
  if(fabs(r[2][2]) < 1.0-err) {
    th1=acos(r[2][2]);
    s1=sin(th1);
    ph1=atan2(r[2][0]/s1,-r[2][1]/s1);
    ps1=atan2(r[0][2]/s1,r[1][2]/s1);
  }
  else {
    th1=0.0;
    ph1=0.0;
    ps1=atan2(r[0][1],r[0][0]);
  }
  euler[0][0]=ph1; euler[0][2]=(double)180.0*ph1/pi;
  euler[1][0]=th1; euler[1][2]=(double)180.0*th1/pi;
  euler[2][0]=ps1; euler[2][2]=(double)180.0*ps1/pi;
  
  printf("\nZXZ: %9.6f %9.6f %9.6f   [DEG: %9.4f %9.4f %9.4f ]",
         euler[0][0],euler[1][0],euler[2][0],euler[0][2],euler[1][2],euler[2][2]);

// TEST ZXZ:
//  long n1, n2;
//  float m[3][3];
//  double c3, s2, s3, z;
/*
    c1 = cos(euler[0][0]); // Phi
    s1 = sin(euler[0][0]);
    c2 = cos(euler[1][0]); // Theta
    s2 = sin(euler[1][0]);
    c3 = cos(euler[2][0]); // Psi
    s3 = sin(euler[2][0]);
    m[0][0] =  c1 * c3 - s1 * c2 * s3;
    m[0][1] =  s1 * c3 + c1 * c2 * s3;
    m[0][2] =  s2 * s3; 
    m[1][0] = -c1 * s3 - s1 * c2 * c3;
    m[1][1] = -s1 * s3 + c1 * c2 * c3;
    m[1][2] =  s2 * c3;
    m[2][0] =  s1 * s2;
    m[2][1] = -c1 * s2;
    m[2][2] =  c2;
  z=0.0;
  for(n1=0;n1<3;n1++)
    for(n2=0;n2<3;n2++)
      z=z+fabs(m[n1][n2]-r[n1][n2]);
  printf("\nTEST ZXZ: error = %10.6f",z);

    printf("\nTEST ZXZ: Unitary ROTATION matrix calculated from euler angles:");
    printf("\n  %10.6f   %10.6f   %10.6f ",
              m[0][0],m[1][0],m[2][0]);
    printf("\n  %10.6f   %10.6f   %10.6f ",
              m[0][1],m[1][1],m[2][1]);
    printf("\n  %10.6f   %10.6f   %10.6f \n",
              m[0][2],m[1][2],m[2][2]);
*/
  printf("\n ");
  return;
}

/*-----------------------------------------------------------
/
/   write_output - write an output; PDB format file
/
/------------------------------------------------------------*/
void write_output(part_struct *part, pdata_struct *pdata,
                  long ca_atoms[2][MAXRES+1],
                  float r[3][3], float v[3], char* fname)
{
  long i, j, k, n, ok, second;
  float x, y, z, Rx, Ry, Rz, euler[3][4];
  char ch, fname_out[200], keyword[200], line[500], l_line[200], r_line[200];
  char sub[10], sub_res[200], res_i[200];
  FILE *fp, *fp_in;

  if(pdata->s_analysis==5 && part->opt_sup_rms[0]>999) {
    printf(" ERROR! The DISTANCE cutoff %7.3f A is too small.\n",part->dist_cutoff);
    return;
  }
  
  sprintf(sub,"%s",pdata->atoms);
  n=ca_atoms[0][0];
  if(pdata->s_analysis==1) {
    sprintf(part->summary.s[5].line,"REMARK   #%-8s  N1   N2  RMSD_global\n",sub);
    sprintf(part->summary.s[6].line,"REMARK   SUMMARY: %4ld %4ld   %6.2f\n",
              part->m_n_aa, part->t_n_aa, part->rmsd);              
    sprintf(part->summary.s[0].line,"\n# Standard rmsd calculations on all %ld assigned %s atoms:", n, pdata->atoms);
    sprintf(part->summary.s[1].line,"\nSUMMARY(RMS): N1 = %4ld  N2 = %4ld  RMSD =%7.3f\n",
              part->m_n_aa, part->t_n_aa, part->rmsd);
    part->summary.n=2;
  }
  if(pdata->s_analysis==2) {
    sprintf(part->summary.s[5].line,"REMARK   #%-8s  N1   N2  DIST    N   RMSD  RMSD_global\n",sub);
    sprintf(part->summary.s[6].line,"REMARK   SUMMARY: %4ld %4ld  %4.1f %4ld %6.2f   %6.2f\n",
              part->m_n_aa, part->t_n_aa, part->dist_cutoff, part->opt_sup1[0][0], part->opt_sup_rms[0], part->rmsd);
    sprintf(part->summary.s[0].line,"\n# Iterative Superposition Procedure on %s atoms:",pdata->atoms);
    sprintf(part->summary.s[1].line,"\nSUMMARY(ISP): N1 = %4ld  N2 = %4ld  RMSD =%7.3f  N = %4ld  DIST = %6.2f",
              part->m_n_aa, part->t_n_aa, part->opt_sup_rms[0], part->opt_sup1[0][0], part->dist_cutoff);
    sprintf(part->summary.s[2].line,"\n\nISP_ASGN_ATOMS RMSD: %7.3f  Number of assigned atoms: %4ld  Nb iter: %3ld ",
              part->rmsd_isp, pdata->n_atoms, part->isp_iter);
    sprintf(part->summary.s[3].line,"\nStd_ASGN_ATOMS RMSD: %7.3f  Standard rmsd on all %ld assigned %s atoms \n",
              part->rmsd, n, pdata->atoms);
    part->summary.n=4;
  }
  if(pdata->s_analysis==5) {
    sprintf(part->summary.s[5].line,"REMARK   #%-8s  N1   N2  DIST  GAP  R_NUM  Seq_Id  RMSD  D_MAX\n",sub);
    sprintf(part->summary.s[6].line,"REMARK   SUMMARY: %4ld %4ld  %4.1f %4ld %6s %7.2f %5.2f %6.2f\n",
              pdata->molecule[0].n_aa, pdata->molecule[1].n_aa, part->dist_cutoff, 
              part->fit_g, part->fit_r, 100.0*part->fit_si/pdata->molecule[1-part->sia].n_aa,
              part->opt_sup_rms[0], part->opt_all_rms[0]);
    sprintf(part->summary.s[0].line,"\n#              N1   N2   DIST    GAP   RES_NUM   Seq_Id     RMSD     D_MAX ");
    sprintf(part->summary.s[1].line,"\nSUMMARY(FIT) %4ld %4ld   %4.1f   %4ld   %7s   %6.2f  %7.2f   %7.2f",
              pdata->molecule[0].n_aa, pdata->molecule[1].n_aa, part->dist_cutoff,
              part->fit_g, part->fit_r, 100.0*part->fit_si/pdata->molecule[1-part->sia].n_aa,
              part->opt_sup_rms[0], part->opt_all_rms[0]);
    part->summary.n=2;
  }
  
  if(pdata->s_analysis<4) {
    ok=print_lga(ca_atoms, ca_atoms, ca_atoms, pdata->molecule, r, v, part, n, n);
  }

  for(k=0;k<part->summary.n;k++) {
    printf("%s",part->summary.s[k].line);
  }
  part->summary.n=0;
  if(pdata->s_analysis>=3 && pdata->s_analysis<=5) {
    if(part->resranges==1) {
      printf("  (%s:%s:%s:%s)",
            part->resrangaa[0].aa,part->resrangaa[1].aa,
            part->resrangaa[2].aa,part->resrangaa[3].aa);
      if(part->stral>0) {
        printf("[%s:%s:%s:%s]",
            part->resrangaa[4].aa,part->resrangaa[5].aa,
            part->resrangaa[6].aa,part->resrangaa[7].aa);
      }
    }
    printf("\n");
  }
  
  if(pdata->s_analysis==3) {
    sprintf(part->summary.s[5].line,"REMARK   #%-8s  N1   N2  DIST    N   RMSD  LGA_S3  RMSD_global\n",sub);
    sprintf(part->summary.s[6].line,"REMARK   SUMMARY: %4ld %4ld  %4.1f %4ld %6.2f %7.3f   %6.2f\n",
              part->m_n_aa, part->t_n_aa, part->dist_cutoff, part->opt_sup1[0][0], part->summary.rms_local, part->summary.lcs_gdt_ts, part->rmsd);
    printf("\nLGA_LOCAL      RMSD: %7.3f  Number of atoms: %4ld  under DIST: %6.2f",
              part->opt_sup_rms[0],part->opt_sup1[0][0],part->dist_cutoff);
    printf("\nLGA_ASGN_ATOMS RMSD: %7.3f  Number of assigned atoms: %4ld ",
              part->rmsd_isp, pdata->n_atoms);
    printf("\nStd_ASGN_ATOMS RMSD: %7.3f  Standard rmsd on all %ld assigned %s atoms \n",
              part->rmsd,n,pdata->atoms);
  }
  if(pdata->s_analysis==4) {
    sprintf(part->summary.s[5].line,"REMARK   #%-8s  N1   N2  DIST    N   RMSD  Seq_Id   LGA_S\n",sub);
    sprintf(part->summary.s[6].line,"REMARK   SUMMARY: %4ld %4ld  %4.1f %4ld %6.2f  %6.2f %7.3f\n",
              part->m_n_aa, part->t_n_aa, part->dist_cutoff, part->opt_sup1[0][0], part->summary.rms_local, part->summary.seq_id, part->summary.lcs_gdt_ts);
  }
  
  printf("\nUnitary ROTATION matrix and the SHIFT vector superimpose molecules  (1=>2)");
  printf("\n  X_new = %10.6f * X  + %10.6f * Y  + %10.6f * Z  + %10.6f",
           r[0][0],r[1][0],r[2][0],v[0]);
  printf("\n  Y_new = %10.6f * X  + %10.6f * Y  + %10.6f * Z  + %10.6f",
           r[0][1],r[1][1],r[2][1],v[1]);
  printf("\n  Z_new = %10.6f * X  + %10.6f * Y  + %10.6f * Z  + %10.6f \n",
           r[0][2],r[1][2],r[2][2],v[2]);

  euler_angles(r,euler);

  if(pdata->s_analysis>=2 && pdata->s_analysis<=4) {
    if(part->lN_n>part->opt_sup1[0][0]) {
      printf("\n# NOTE! The number N = %ld of calculated residues is smaller then cutoff lN:n = %ld \n",
             part->opt_sup1[0][0],part->lN_n);
      pdata->one_output=0;
    }
  }

  if((fp_in = fopen(pdata->fname_in,"r"))==NULL) {
    printf("\n# ERROR! opening file %s for read\n\n",pdata->fname_in);
    exit(0);
  }
  
  if(pdata->one_output==0) return;
/*-----------------------------------------------------------
/
/   write rotated structure in the PDB format
/
/------------------------------------------------------------*/
  
  ok=0;
  second=0;
  x=0.0;y=0.0;z=0.0;
  if(part->sia==1) {
    for(k=0;k<3;k++) {
      r[0][k]=part->opt_r1[0][0][k];
      r[1][k]=part->opt_r1[0][1][k];
      r[2][k]=part->opt_r1[0][2][k];
      v[k]=part->opt_v1[0][k];
    }
  }

  strcpy(fname_out,fname);
  strcat(fname_out,".pdb");
  if((fp = fopen(fname_out,"w"))==NULL) {
    printf("\n*** ERROR! opening file %s for write ***\n",fname_out);
    exit(0);
  }
  fprintf(fp,"REMARK  ---------------------------------------------------------- \n");
  fprintf(fp,"REMARK   Citing LGA: \n");
  fprintf(fp,"REMARK   Zemla A., LGA - a Method for Finding 3D Similarities in  \n");
  fprintf(fp,"REMARK   Protein Structures, Nucleic Acids Research, 2003, V. 31, \n");
  fprintf(fp,"REMARK   No. 13, pp. 3370-3374. \n");
  fprintf(fp,"REMARK  ---------------------------------------------------------- \n");
  if(pdata->one_output==0) {
    fprintf(fp,"REMARK   Superimposed MOLECULES (1=>2)   Output: no coordinates  \n");
  }
  if(pdata->one_output==1) {
    fprintf(fp,"REMARK   Superimposed MOLECULES (1=>2)   Output: 1 molecule      \n");
  }
  if(pdata->one_output==2) {
    fprintf(fp,"REMARK   Superimposed MOLECULES (1=>2)   Output: 2 molecules     \n");
  }
  fprintf(fp,"REMARK     1: %-32s              \n",pdata->molecule[0].molecule_name);
  fprintf(fp,"REMARK     2: %-32s              \n",pdata->molecule[1].molecule_name);
  if(pdata->s_analysis==1) {
    fprintf(fp,"REMARK   Standard RMSD calculation\n");
  }
  else if(pdata->s_analysis==2) {
    fprintf(fp,"REMARK   Iterative Superposition Procedure\n");
  }
  else if(pdata->s_analysis==3) {
    fprintf(fp,"REMARK   GDT and LCS analysis\n");
  }
  else if(pdata->s_analysis==4) {
    fprintf(fp,"REMARK   Structure alignment analysis\n");
  }
  else if(pdata->s_analysis==5) {
    fprintf(fp,"REMARK   Structure best fit analysis\n");
  }
  if(pdata->s_analysis==4) {
    fprintf(fp,"REMARK   Search for Atom-Atom correspondence\n");
  }
  else {
    fprintf(fp,"REMARK   FIXED Atom-Atom correspondence\n");
  }
  fprintf(fp,"REMARK   ");
  for(k=0;k<pdata->parameters_n;k++) {
    fprintf(fp,"%s",pdata->parameters[k].line);
  }
  fprintf(fp,"\n%s",part->summary.s[5].line);
  fprintf(fp,"%s",part->summary.s[6].line);
  fprintf(fp,"REMARK  ---------------------------------------------------------- \n");

  while(fgets(line,500,fp_in)!=NULL && ok<3) {
    for(i=0;i<500;i++) if(line[i]=='%') line[i]=' ';
    for(i=0;i<50;i++) keyword[i]=0;
    sscanf(line,"%s",keyword);
    if(!strncmp(keyword,"MOLECULE\0",9)) {
      ok++;
      if(ok>1 && pdata->one_output==1) {
        ok=3;
      }
      else {
        fprintf(fp,"%s",line);
      }
    }
    else if((!strncmp(keyword,"ATOM\0",5) ||
            !strncmp(keyword,"HETATM",6)) && ok==1) {
      i=0;
      while(line[i]!='\n' && i<500) i++;
      if(i>=54) {
        if(part->sia==1 && second<2) {
          for(j=0;j<10;j++) sub[j]=0;
          strncpy(sub,&line[22],5);
          sscanf(sub,"%s",sub_res);
          for(j=0;j<10;j++) sub[j]=0;
          ch=line[21];
          if(ch!=' ') {
            sub[0]='_';
            sub[1]=ch;
            sub[2]='\0';
          }
          strcat(sub_res,sub);
          sscanf(sub_res,"%s",res_i);
          if(!strcmp(res_i,part->fit_r)) {
            second=1;
          }
          else if(second==1) {
            second=2;
            for(k=0;k<3;k++) {
              r[0][k]=part->opt_r2[0][0][k];
              r[1][k]=part->opt_r2[0][1][k];
              r[2][k]=part->opt_r2[0][2][k];
              v[k]=part->opt_v2[0][k];
            }
          }
        }

        for(j=0;j<120;j++) {
          l_line[j]=0;
          r_line[j]=0;
        }
        i=i-54;
        strncpy(l_line,&line[0],30);
        strncpy(r_line,&line[54],i);

        for(j=0;j<10;j++) sub[j]=0;
        strncpy(sub,&line[30],8);
        sscanf(sub,"%f",&Rx);

        for(j=0;j<10;j++) sub[j]=0;
        strncpy(sub,&line[38],8);
        sscanf(sub,"%f",&Ry);

        for(j=0;j<10;j++) sub[j]=0;
        strncpy(sub,&line[46],8);
        sscanf(sub,"%f",&Rz);

        x=r[0][0]*Rx + r[1][0]*Ry + r[2][0]*Rz + v[0];
        y=r[0][1]*Rx + r[1][1]*Ry + r[2][1]*Rz + v[1];
        z=r[0][2]*Rx + r[1][2]*Ry + r[2][2]*Rz + v[2];

        fprintf(fp,"%-30s%8.3f%8.3f%8.3f%-s\n",
                l_line, x, y, z, r_line);
      }
      else {
        fprintf(fp,"%s",line);
      }
    }
    else {
      if(ok>0) {
        fprintf(fp,"%s",line);
      }
    }
  }
  fclose(fp);
  fclose(fp_in);

  return;
}
