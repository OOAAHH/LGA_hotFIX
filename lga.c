/*------------------------------------------------------------------------------
/
/                                 LGA
/                      (Local-Global Alignment)
/      A Method for Finding 3-D Similarities in Protein Structures
/                    ----------------------------
/                            ver. 09/2019
/                     (LGA first ver. 10/15/1999)
/
# Author: Adam Zemla 
# US Patent: 8024127
# Copyright: CP01155
# All rights reserved. 
# Email: adamz@llnl.gov
#
# The GDT, LCS, and the LGA program and documentation may not be distributed,
# sold or incorporated longo a commercial product, in whole or in part,
# without written consent of Adam T. Zemla and the University
# of California, Lawrence Livermore National Laboratory (UC LLNL).
#
/-------------------------------------------------------------------------------
#
# This work was produced at the University of California, Lawrence
# Livermore National Laboratory (UC LLNL) under contract
# no. W-7405-ENG-48 (Contract 48) between the U.S. Department of Energy
# (DOE) and The Regents of the University of California (University) for
# the operation of UC LLNL. Copyright is reserved to the University for
# purposes of controlled dissemination, commercialization through formal
# licensing, or other disposition under terms of Contract 48; DOE
# policies, regulations and orders; and U.S. statutes. The rights of the
# Federal Government are reserved under Contract 48 subject to the
# restrictions agreed upon by the DOE and University as allowed under
# DOE Acquisition Letter 97-1.
#
# DISCLAIMER
#
# This work was prepared as an account of work sponsored by an agency of
# the United States Government. Neither the United States Government nor
# the University of California nor any of their employees, makes any
# warranty, express or implied, or assumes any liability or
# responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents
# that its use would not infringe privately-owned rights.  Reference
# herein to any specific commercial products, process, or service by
# trade name, trademark, manufacturer or otherwise does not necessarily
# constitute or imply its endorsement, recommendation, or favoring by
# the United States Government or the University of California. The
# views and opinions of authors expressed herein do not necessarily
# state or reflect those of the United States Government or the
# University of California, and shall not be used for advertising or
# product endorsement purposes.
#
/-------------------------------------------------------------------------------
/
/ Usage: lga [<options>] file_name
/
/   <options>: [ -h | -aa | -al | -batch ]
/              [ -1 | -2 | -3 | -4 | -5 ]
/              [ -mol1:name1 | -mol2:name2 ]
/              [ -sda | -sia | -fit:b:gap:res | -stral | -stral:f | -lN:n ]
/              [ -atom:CA | -bmo:b:m:o | -cb:f | -ah:i | -ch1:c | -ch2:c ]
/              [ -aa1:n1:n2 | -aa2:n1:n2 | -gap1:n1:n2 | -gap2:n1:n2 ]
/              [ -er1:s1:s2 | -er2:s1:s2 ]
/              [ -gdc_set:s1:s2 | -gdc_sup:s1:s2 | -gdc_at:a1,a2 | -gdc_eat:e1:e2 ]
/              [ -gdc_sc | -gdc | -gdc:n | -gdc_ref | -gdc_ref:n ]
/              [ -o0 | -o1 | -o2 | -r | -rmsd | -swap | -fp | -ie ]
/              [ -d:f | -gdt | -lw:n | -lga:f | -lga_m ]
/
/   where: -h              help information
/          -1              standard RMSD
/          -2              RMSD using ISP (Iterative Superposition Procedure)
/          -3              GDT and LCS analysis
/          -4              structure alignment analysis
/          -5              structure best fit analysis: S => S-gap-S , S-gap-S => S
/          -mol1:name1     name of the molecule1 that will be used in output file
/                            (name1.name2). The alphanumeric characters and '_' 
/                            are allowed.
/          -mol2:name2     name of the molecule2 that will be used in output file
/                            (name1.name2). The alphanumeric characters and '_'
/                            are allowed.
/          -atom:CA        CA (Calpha) atoms will be used for calculations. 
/                            NOTE: to specify special character "'" use ",". 
/                            For example: use "-atom:CB" to select CB atom,
/                            use "-atom:H5,1" to select H5'1 atom.
/          -bmo:b:m:o      CB and M (M = Mid: C,CA,CB,N) atoms will be calculated.
/                            The coordinates of the polong representing amino-acid    
/                            position (BMO; backbone model) for LGA processing are 
/                            defined by the following vectors:
/                            vector CA-CB: -5.0 <= b <= 5.0
/                            vector CA-M:  -5.0 <= m <= 5.0
/                            vector CA-O:  -5.0 <= o <= 5.0
/                            For example: CA = -bmo:0:0:0 (default)
/          -cb:f           CB (Cbeta) atom position will be calculated for each amino-acid,
/                            and the coordinates of the polong representing amino-acid    
/                            position (BMO; backbone model) for LGA processing will be 
/                            defined by the vector CA-CB: -5.0 <= f <= 5.0 , (e.g. f=0 
/                            corresponds to CA position, and f=1 represents CB position)
/                            NOTE: if "-cb:f" is combined with "-atom:CB" then all 
/                            existing CB atoms are leveraged and only missing CB atoms 
/                            are calculated. This option is equvalent to -bmo:f:0:0
/          -ch1:c          chain c selected from molecule1
/          -ch2:c          chain c selected from molecule2
/          -ah:i           ATOM or HETATM records are used for calculations:
/                            i=0 both
/                            i=1 ATOM
/                            i=2 HETATM
/          -lga:f          weight LCS and GDT measures: LCS 0.0 <= f <= 1.0 GDT,
/                            default: f = 0.75
/          -lga_m          maximum value of LGA_S (LGA_M) reported in SUMMARY line
/          -d:f            DIST distance cutoff (f Angstroms; default f=5.0)
/          -opt:n          optimization parameter: 0, 1, 2. Default: 1 
/          -gdt            can be combined with "-3" option. If used then the 
/                            superposition that fits maximum number of residues under 
/                            a given distance cutoff is reported. Otherwise standard 
/                            superposition calculated using the set of identified N 
/                            residues is reported (rotated molecule1)
/          -lw:n           "Lesk window", rms calculated on residue window
/                            (length of the window = 2*n+1)
/          -lN:n           limit on N superimposed residues (if calculated N<n,
/                            then output is not generated). Useful for database 
/                            calculations.
/          -sda            sequence dependent analysis
/                            structure conformation dependent analysis (S => S-gap-S)
/          -sia            sequence independent analysis,
/                            structure conformation independent analysis (S-gap-S => S)
/          -fit:b:g:r      search for the best fit,
/               b          number of residues below 0.5 A (b - longeger 0 <= b <= 9)
/               g          length of the gap (g - longeger 0 <= g <= 99)
/               r          residue number after which the gap appears (r - string)
/          -er1:s1:s2      exact range of residues from the molecule1 used for 
/                            calculations (s1 , s2 - strings e.g.: s1 = 13L_A <= s2 = 45_B)
/                            the si pairs (ranges beg:end) can be separated by ',':
/                              -er1:s1:s2,s3:s4,s5:s6,s7:s8,s9:s10 
/                            NOTE: single residues or chains can be separated by ','(no beg:end required):
/                              -er1:s1,s2,s3,
/                            Up to 50 er1 parameters are allowed (WARNING: no overlaps)
/          -er2:s1:s2      exact range of residues from the molecule2 used for 
/                            calculations (s1 , s2 - strings e.g.: s1 = 13L_A <= s2 = 45_B)
/                            the si pairs (ranges beg:end) can be separated by ',':
/                              -er2:s1:s2,s3:s4,s5:s6,s7:s8,s9:s10 
/                            NOTE: single residues or chains can be separated by ','(no beg:end required):
/                              -er2:s1,s2,s3,
/                            Up to 50 er2 parameters are allowed (WARNING: no overlaps)
/          -gdc:n          n - number of bins used for GDC evaluation of atom pairs from the
/                            corresponding residues (1 <= n <= 20; bins: <0.5, <1.0, ... <10.0).
/                            NOTE: this option changes the default number of "bins" (n=20) for
/                            GDC calculations (GDC_all - all atoms, GDC_mc - main chain atoms, 
/                            and GDT_at - selected atoms). The default number n=20 defines bins 
/                            from 0.5 to 10.0 Angstroms.
/          -gdc            GDC score is calculated using all identical atoms from the target as 
/                            a frame of reference (equivalent to: -gdc_ref:2 -swap)
/          -gdc_ref:n      GDC score is calculated: 0 - requesting a complete set of atoms within 
/                            each residue, 1 - using atoms from the target as a frame of reference,
/                            2 - using all identical atoms from the target as a frame of reference.
/                            The default set is -gdc_ref:0
/          -gdc_sup:s1:s2  exact range of residues from the molecule2 used for 
/                            GDC superposition calculations. This additional standard (-1) 
/                            superposition is calculated on CA atoms from the set of 
/                            amino-acid ranges (s1,s2) defined by s1 and s2 strings. 
/                            e.g. -gdc_sup:s1:s2,s3:s4,s5:s6,s7:s8,s9:s10 
/                            Format is the same as for er2 parameters.
/                            NOTE: this option is applied to the molecule2 only. Corresponding
/                            residues from molecule1 are automatically determined using main 
/                            superposition. 
/          -gdc_sup        expands an option "-rmsd". If used then the superposition which is 
/                            used for GDC calculations is reported and used to rotate molecule1.
/                            Otherwise the standard LGA superposition is reported.
/          -gdc_set:s1:s2  exact range of residues from the molecule2 for which the
/                            "Global Distance Calculations" (GDC) will be performed.
/                            e.g. -gdc_set:s1:s2,s3:s4,s5:s6,s7:s8,s9:s10 
/                            Format is the same as for er2 parameters.
/                            NOTE: this option is applied to the molecule2 only. Amino-acids
/                            from the molecule2 serve as a frame of reference for GDC evaluation
/                            (corresponding amino-acids or atoms that are missing in molecule1
/                            are counted as 0 scores in GDC calculations).
/          -gdc_at:a1,a2   amino-acid atom names (one atom per one name of amino-acid) from 
/                            the molecule2 for which the GDC calculations (distances and GDC 
/                            summary) will be calculated. 
/                            Format example (aaname.atom): -gdc_at:a1,a2,a3,a4
/                            where: a1 = V.CG1, a2 = C.SG, a3 = T.OG1, a4 = H.NE2
/                            NOTE: this option is applied to the molecule2 only. The 
/                            corresponding atoms from the molecule1 will be detected based 
/                            on the calculated alignment. Up to 20 representative atoms 
/                            (one atom per each of 20 amino-acid) can be selected for 
/                            GDC evaluation. Number of identified identical "amino-acid.atom"
/                            pairs serve as a frame of reference for GDC evaluation.
/                            Results from the GDC-at calculations are reported in Dist_at and 
/                            GDC_at columns. 
/          -gdc_at:*.at    allows a selection of one mainchain or CB atom (at: N,CA,C,O,CB)
/                            the same for all amino acids (e.g. -gdc_at:*.N).
/                            NOTE: amino-acids from the molecule2 serve as a frame of reference
/                            for GDC evaluation (corresponding amino-acids or atoms that are
/                            missing in molecule1 are counted as 0 scores in GDC calculations).
/          -gdc_eat:e1:e2  exact atom "e1" from the molecule1 and "e2" from the molecule2 for 
/                            which the GDC calculations (distances and GDC summary) will be 
/                            calculated. Format example (aanumber_chain.atom):
/                            -gdc_eat:e1:e2,e3:e4,e5:e6
/                            where: for each pair (em:en) em is a selected atom from the 
/                            molecule1, and en is an atom from the molecule2.
/                            For example: e1 = 10_A.OD2, e2 = 21_B.ND2
/          -gdc_sc         automated selection of all flags required for GDC_sc calculations: 
/                            -swap -gdc:10 -gdc_at:V.CG1,L.CD1,I.CD1,P.CG,M.CE,F.CZ,W.CH2,S.OG
/                            -gdc_at:T.OG1,C.SG,Y.OH,N.OD1,Q.OE1,D.OD2,E.OE2,K.NZ,R.NH2,H.NE2
/                            NOTE: this option changes the default number of "bins" (see the
/                            selection "-gdc:n"; n=10). All GDC calculations (GDC_all - all atoms, 
/                            GDC_mc - main chain atoms, and GDT_at - selected atoms) will be 
/                            performed using n=10 as a number of bins from 0.5 to 5.0 Angstroms.
/                            Results from the GDC_sc calculations are reported in GDC_at column. 
/          -aa1:n1:n2      range of residues from the molecule1 used for calculations
/                            -9999 < n1 < n2 < 9999 (n1, n2 - longeger) 
/                            NOTE: only one aa1 parameter is allowed.
/          -aa2:n1:n2      range of residues from the molecule2 used for calculations
/                            -9999 < n1 < n2 < 9999 (n1, n2 - longeger) 
/                            NOTE: only one aa2 parameter is allowed.
/          -gap1:n1:n2     range of residues from the molecule1 removed from calculations
/                            -9999 < n1 < n2 < 9999 (n1, n2 - longeger) 
/                            NOTE: only one gap1 parameter is allowed.
/          -gap2:n1:n2     range of residues from the molecule2 removed from calculations
/                            -9999 < n1 < n2 < 9999 (n1, n2 - longeger) 
/                            NOTE: only one gap2 parameter is allowed.
/          -aa             generates a list of all residues from the molecule1 and
/                            (molecule2 AAMOL* records)
/          -al             calculations will be made only on the set of residues from the
/                            attached AAMOL* or LGA records
/          -o0             no coordinates are printed out
/          -o1             only molecule 1 (rotated) is printed out longo the
/                            subdirectory TMP
/          -o2             molecule 1 (rotated) and molecule 2 (target) both are
/                            printed out longo the subdirectory TMP
/          -r              the residue ranges of compared structures are reported in the
/                            SUMMARY line: e.g. (1_A:214_A:7_A:196_A)
/          -rmsd           additional RMSD and GDC calculations will be performed on all
/                            aligned CA, MC and ALL atoms. 
/                            RMSD is "rmsd-based" measures: see MC and ALL colums
/                            GDC is "distance-based" measures: see Dist_max, GDC_mc, and GDC_all
/          -swap           expands an option "-rmsd". RMSD and GDC calculations will be 
/                            performed with checking for swapping atoms in amino acids: 
/                            ASP, GLU, PHE, and TYR
/          -fp             full print output
/          -check          reports amino acids with missing pre-selected atoms
/          -ie             ignores errors in PDB data (force calculations). If "-ie" not 
/                            present then in case of ERROR detected in input data the 
/                            calculations are terminated
/          -stral          additional information about identified structural SPANS (regions with tight 
/                            superpositions) is reported: S_nb - number of SPANS, S_N - combined number 
/                            of residues within SPANS, S_Id - average sequence identity within SPANS 
/                            (standalone version: two output files in TMP directory are created: *.stral and *.pdb) 
/          -stral:f        cutoff for local RMSD for stral calculations (0.01 <= f <= 10.0)
/                            default: f = 0.5
/          -batch:frun     it allows to run several different lga calculations on the 
/                            same mol1.mol2 pair of structures. File frun contains a list 
/                            of parameters. Maximum number of RUN lines is limited to 400
/                            (see below).
/                            
/ Example of the frun file (NOTE: in the example below (fit search) the mol2 
/ has to be chain ' '):
/   STRUCTURES: mol1.mol2
/   RUN: mol1.mol2_L00 -5 -sda -o2 -d:2 -lw:2 -fit:0:2:21 -aa2:19:26 -gap2:22:23 
/   RUN: mol1.mol2_L10 -5 -sda -o2 -d:2 -lw:2 -fit:0:3:20 -aa2:18:26 -gap2:21:23 
/   RUN: mol1.mol2_L20 -5 -sda -o2 -d:2 -lw:2 -fit:0:4:19 -aa2:17:26 -gap2:20:23 
/   RUN: mol1.mol2_L30 -5 -sda -o2 -d:2 -lw:2 -fit:0:5:18 -aa2:16:26 -gap2:19:23 
/   RUN: mol1.mol2_L40 -5 -sda -o2 -d:2 -lw:2 -fit:0:6:17 -aa2:15:26 -gap2:18:23 
/   RUN: mol1.mol2_L50 -5 -sda -o2 -d:2 -lw:2 -fit:0:7:16 -aa2:14:26 -gap2:17:23 
/   RUN: mol1.mol2_L01 -5 -sda -o2 -d:2 -lw:2 -fit:0:3:21 -aa2:19:27 -gap2:22:24 
/   RUN: mol1.mol2_L11 -5 -sda -o2 -d:2 -lw:2 -fit:0:4:20 -aa2:18:27 -gap2:21:24 
/ 
/ Another example of the frun file:  
/   STRUCTURES: 1wdn_A.1ggg_A
/   RUN: segment_1_50 -4 -o0 -d:3 -swap -stral -aa2:1:50 
/   RUN: segment_26_75 -4 -o0 -d:3 -swap -stral -aa2:26:75 
/   RUN: segment_51_100 -4 -o0 -d:3 -swap -stral -aa2:51:100 
/   RUN: segment_76_125 -4 -o0 -d:3 -swap -stral -aa2:76:125 
/   RUN: segment_101_150 -4 -o0 -d:3 -swap -stral -aa2:101:150 
/   RUN: segment_126_175 -4 -o0 -d:3 -swap -stral -aa2:126:175 
/   RUN: segment_151_200 -4 -o0 -d:3 -swap -stral -aa2:151:200 
/   RUN: segment_176_225 -4 -o0 -d:3 -swap -stral -aa2:176:225 
/                              
/ NOTE: within the batch file "frun" the following parameters are ignored: -bmo, -cb, -atom. 
/   It is because the set of coordinates for the calculations is loaded only once from 
/   the file defined by the STRUCTURES record.
/
/ Default set of parameters: -1 , -d:5.0, -sia, -o1
/
/ Input: <file_name>
/   <file_name> - the file is located inside the subdirectory MOL2, and
/                 contains two structures in PDB format. 
/                 Each structure for LGA analysis should begin with 
/                 MOLECULE and end with END record:
/
/              MOLECULE  name1
/              ATOM ........
/                   ........
/              ATOM ........
/              END
/              MOLECULE  name2
/              ATOM ........
/                   ........
/              ATOM ........
/              END
/
/   Input files are located inside the subdirectory MOL2 (check dir_in below).
/
/ Output:  <file_name.pdb>, <file_name.lga>
/   <file_name.pdb> - contains two superimposed PDB structures: 1 => 2
/   <file_name.lga> - contains calculated residue equivalences
/   NOTE: if options: -mol1:name1 and -mol2:name2 are used 
/         then output file_name = name1.name2
/
/   Output files are written longo the subdirectory TMP (check dir_out below).
/
/-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "tools.h"

int main(int argc, char *argv[]) 
{
  typedef struct {
    char        v[1000];
    char        line[1000];
  } input_line;
  input_line   Arg_input[400];   // MAX number of RUN lines (see batch): 400
  long          i, j, k, n_atoms, ok, ignore_aa, Arg_input_c;
  long          n, n_ok, n_end, batch_number, Total_batch_number, batch_mode, sda_selector;
  char         batch_name[200], fname[200], dir_in[200], dir_out[200];
  char         fit[200], range[1000], mname[200], line[1000], keyword[200], *word;
  long        ca_atoms[2][MAXRES+1];
  float        r[3][3], v[3];
  /* NOTE: These structs are very large. Declaring them static moves
     them to static storage (BSS) instead of the stack, avoiding
     stack overflows on systems with smaller default stack limits. */
  static pdata_struct pdata;
  static part_struct  part;
  mol2_aa      aa_mol2[2];
  FILE         *fp;

  printf("\n#######################################################");
  printf("\n#                                                     #");
  printf("\n#                        LGA                          #");
  printf("\n#                  ---------------                    #");
  printf("\n#                                                     #");
  printf("\n#               Local-Global Alignment                #");
  printf("\n#        A Method for Finding 3-D Similarities        #");
  printf("\n#               in Protein Structures                 #");
  printf("\n#                                                     #");
  printf("\n#                  ------------ 09/2019               #");
  printf("\n#                                                     #");
  printf("\n#      Adam Zemla (zemla1@llnl.gov)                   #");
  printf("\n#      Lawrence Livermore National Laboratory, CA     #");
  printf("\n#                                                     #");
  printf("\n#######################################################\n");

  i=0; // i=1 for printing an expiration date
  ok=lchk(i);
  if(ok!=1) {
    printf("\n!!!   License for LGA program has expired  !!!");
    printf("\n#   For updated version of the LGA program   #");
    printf("\n# please contact Adam Zemla: zemla1@llnl.gov #\n");
//    exit(0);
  }

  clean_pdata(&pdata);

  strcpy(dir_in,"MOL2/\0");
  strcpy(dir_out,"TMP/\0");

  n_atoms=0;
  ignore_aa=1;
  sda_selector=0; // sda used as residue selector: 0 - calculations , 1 - calculations + SUMMARY (applied to: "gapn, aan, ern")
  batch_mode=0;  // 0 - no batch mode, 1 - batch (reading structures), 2 - batch (multiple runs)
  batch_number=1;
  Total_batch_number=batch_number+1;
  Arg_input_c=argc;
  for(i=1;i<Arg_input_c;i++) {
    strcpy(Arg_input[i].v,argv[i]);
  }

  /* Initialize heap workspace once */
  if (part.ws == NULL) {
    if (lga_ws_init(&part.ws) != 0) {
      fprintf(stderr, "\n# ERROR! Failed to allocate LGA workspace\n");
      exit(1);
    }
  }

  while(batch_number>0) {
    clean_part(&part);
    ok=0;
    batch_number--;
//    printf("\n#batch test %ld\n",batch_number);
    
    if(Arg_input_c<2){
      printf("\nUsage: lga [<options>] file_name\n\n");
      exit(0);
    }
    for(i=1;i<Arg_input_c;i++) {
      if(!strncmp(Arg_input[i].v,"-h\0",3) ||
         !strncmp(Arg_input[i].v,"-help\0",6)) {
        if((fp = fopen("README.lga","r"))==NULL) {
          printf("\n# ERROR! opening README.lga file\n");
          exit(0);
        }
        while(fgets(line,1000,fp)!=NULL) {
          printf("%s",line);
        }
        fclose(fp);
        printf("\n");
        exit(0);
      }
      else if(!strncmp(Arg_input[i].v,"-r\0",3)) {
        part.resranges=1;
      }
      else if(!strncmp(Arg_input[i].v,"-o0\0",4)) {
        pdata.one_output=0;
      }
      else if(!strncmp(Arg_input[i].v,"-o1\0",4)) {
        pdata.one_output=1;
      }
      else if(!strncmp(Arg_input[i].v,"-o2\0",4)) {
        pdata.one_output=2;
      }
      else if(!strncmp(Arg_input[i].v,"-fp\0",4)) {
        part.full_print=1;
      }
      else if(!strncmp(Arg_input[i].v,"-ie\0",4)) {
        pdata.ignore_errors=1;
      }
      else if(!strncmp(Arg_input[i].v,"-al\0",4)) {
        ignore_aa=0;
      }
      else if(!strncmp(Arg_input[i].v,"-aa\0",4)) {
        pdata.s_analysis=0;
      }
      else if(!strncmp(Arg_input[i].v,"-opt:",5)) {
        sscanf(Arg_input[i].v+5,"%ld",&j);
        if(j==0) {
          part.accuracy_opt=0;
          part.accuracy_gdt_n=14;
          part.accuracy_gdt_step=0.4;
          part.accuracy_lga_n=4;
          part.accuracy_lga_step=1.4;
        }
        else if(j==2) {
          part.accuracy_opt=2;
          part.accuracy_gdt_n=80;
          part.accuracy_gdt_step=0.07;
          part.accuracy_lga_n=14;
          part.accuracy_lga_step=0.4;
        }
        else { // default values: j=1
          part.accuracy_opt=1;
          part.accuracy_gdt_n=56;
          part.accuracy_gdt_step=0.1;
          part.accuracy_lga_n=8;
          part.accuracy_lga_step=0.7;
        }
      }
      else if(!strncmp(Arg_input[i].v,"-bmo:",5) || 
              !strncmp(Arg_input[i].v,"-BMO:",5) ||
              !strncmp(Arg_input[i].v,"-cb:",4)) {
        pdata.cb_calc=1;
        if(!strncmp(Arg_input[i].v,"-cb:",4)) {
          sscanf(Arg_input[i].v+4,"%f",&pdata.cb_v);
          pdata.m_v=0.0;
          pdata.o_v=0.0;
        }
        else {
          sscanf(Arg_input[i].v+5,"%s",range);
          k=strlen(range);
          for(j=0;j<k;j++) {
            if(range[j]==':') range[j]=' ';
          }
          sscanf(range,"%f %f %f",&pdata.cb_v,&pdata.m_v,&pdata.o_v);
        }
        strcpy(pdata.atoms,"BMO\0");
        pdata.len_atom=strlen(pdata.atoms);
        pdata.len_atom++;
        pdata.ca_atoms=0;
      }
      else if(!strncmp(Arg_input[i].v,"-sia\0",5)) {
        pdata.sda=0;
      }
      else if(!strncmp(Arg_input[i].v,"-sda\0",5)) {
        pdata.sda=1;
      }
      else if(!strncmp(Arg_input[i].v,"-gdt\0",5)) {
        part.gdt=1;
      }
      else if(!strncmp(Arg_input[i].v,"-gdc_sup\0",9)) {
        part.gdc=1;
        part.all_rmsd=1;
      }
      else if(!strncmp(Arg_input[i].v,"-gdc:",5)) {
        part.all_rmsd=1;
        sscanf(Arg_input[i].v+5,"%ld",&part.gdc_bin);
        if(part.gdc_bin<1 || part.gdc_bin>20) {
          part.gdc_bin=20;
        }
      }
      else if(!strncmp(Arg_input[i].v,"-gdc\0",5) ||
              !strncmp(Arg_input[i].v,"-gdc_ref\0",9)) {
        part.swap=1;
        part.all_rmsd=1;
        part.gdc_ref=2;
      }
      else if(!strncmp(Arg_input[i].v,"-gdc_ref:",9)) {
        part.all_rmsd=1;
        sscanf(Arg_input[i].v+9,"%ld",&part.gdc_ref);
        if(part.gdc_ref<0 || part.gdc_ref>2) {
          part.gdc_ref=0;
        }
      }
      else if(!strncmp(Arg_input[i].v,"-gdc_sc\0",8) || 
              !strncmp(Arg_input[i].v,"-GDC_sc\0",8) ||
              !strncmp(Arg_input[i].v,"-gdc_SC\0",8) ||
              !strncmp(Arg_input[i].v,"-GDC_SC\0",8)) {
        part.gdc_bin=10;
        part.swap=1;
        part.all_rmsd=1;
        part.gdc_at.ok=0;
        part.gdc_at.n=0;
        part.eval=1;
        strcpy(part.gdc_at.b[part.gdc_at.n].aa,"V\0");
        strcpy(part.gdc_at.e[part.gdc_at.n].aa,"CG1\0"); part.gdc_at.n++;
        strcpy(part.gdc_at.b[part.gdc_at.n].aa,"L\0");
        strcpy(part.gdc_at.e[part.gdc_at.n].aa,"CD1\0"); part.gdc_at.n++;
        strcpy(part.gdc_at.b[part.gdc_at.n].aa,"I\0");
        strcpy(part.gdc_at.e[part.gdc_at.n].aa,"CD1\0"); part.gdc_at.n++;
        strcpy(part.gdc_at.b[part.gdc_at.n].aa,"P\0");
        strcpy(part.gdc_at.e[part.gdc_at.n].aa,"CG\0"); part.gdc_at.n++;
        strcpy(part.gdc_at.b[part.gdc_at.n].aa,"M\0");
        strcpy(part.gdc_at.e[part.gdc_at.n].aa,"CE\0"); part.gdc_at.n++;
        strcpy(part.gdc_at.b[part.gdc_at.n].aa,"F\0");
        strcpy(part.gdc_at.e[part.gdc_at.n].aa,"CZ\0"); part.gdc_at.n++;
        strcpy(part.gdc_at.b[part.gdc_at.n].aa,"W\0");
        strcpy(part.gdc_at.e[part.gdc_at.n].aa,"CH2\0"); part.gdc_at.n++;
        strcpy(part.gdc_at.b[part.gdc_at.n].aa,"S\0");
        strcpy(part.gdc_at.e[part.gdc_at.n].aa,"OG\0"); part.gdc_at.n++;
        strcpy(part.gdc_at.b[part.gdc_at.n].aa,"T\0");
        strcpy(part.gdc_at.e[part.gdc_at.n].aa,"OG1\0"); part.gdc_at.n++;
        strcpy(part.gdc_at.b[part.gdc_at.n].aa,"C\0");
        strcpy(part.gdc_at.e[part.gdc_at.n].aa,"SG\0"); part.gdc_at.n++;
        strcpy(part.gdc_at.b[part.gdc_at.n].aa,"Y\0");
        strcpy(part.gdc_at.e[part.gdc_at.n].aa,"OH\0"); part.gdc_at.n++;
        strcpy(part.gdc_at.b[part.gdc_at.n].aa,"N\0");
        strcpy(part.gdc_at.e[part.gdc_at.n].aa,"OD1\0"); part.gdc_at.n++;
        strcpy(part.gdc_at.b[part.gdc_at.n].aa,"Q\0");
        strcpy(part.gdc_at.e[part.gdc_at.n].aa,"OE1\0"); part.gdc_at.n++;
        strcpy(part.gdc_at.b[part.gdc_at.n].aa,"D\0");
        strcpy(part.gdc_at.e[part.gdc_at.n].aa,"OD2\0"); part.gdc_at.n++;
        strcpy(part.gdc_at.b[part.gdc_at.n].aa,"E\0");
        strcpy(part.gdc_at.e[part.gdc_at.n].aa,"OE2\0"); part.gdc_at.n++;
        strcpy(part.gdc_at.b[part.gdc_at.n].aa,"K\0");
        strcpy(part.gdc_at.e[part.gdc_at.n].aa,"NZ\0"); part.gdc_at.n++;
        strcpy(part.gdc_at.b[part.gdc_at.n].aa,"R\0");
        strcpy(part.gdc_at.e[part.gdc_at.n].aa,"NH2\0"); part.gdc_at.n++;
        strcpy(part.gdc_at.b[part.gdc_at.n].aa,"H\0");
        strcpy(part.gdc_at.e[part.gdc_at.n].aa,"NE2\0"); part.gdc_at.n++;
      }
      else if(!strncmp(Arg_input[i].v,"-stral\0",7)) {
        part.stral=1;
        part.stral_r=0.5;
      }
      else if(!strncmp(Arg_input[i].v,"-stral:",7)) {
        part.stral=1;
        sscanf(Arg_input[i].v+7,"%f",&part.stral_r);
        if(part.stral_r<0.01 || part.stral_r>10.0) {
          part.stral_r=0.5;
        }
      }
      else if(!strncmp(Arg_input[i].v,"-check\0",7)) {
        part.check=1;
      }
      else if(!strncmp(Arg_input[i].v,"-rmsd\0",6)) {
        part.all_rmsd=1;
        part.swap=0;
      }
      else if(!strncmp(Arg_input[i].v,"-swap\0",6)) {
        part.swap=1;
        part.all_rmsd=1;
      }
      else if(!strncmp(Arg_input[i].v,"-lga_m\0",7) ||
              !strncmp(Arg_input[i].v,"-lga_M\0",7) ||
              !strncmp(Arg_input[i].v,"-LGA_m\0",7) ||
              !strncmp(Arg_input[i].v,"-LGA_M\0",7)) {
        part.lga_m=1;
      }
      else if(!strncmp(Arg_input[i].v,"-gap1:",6)) {
        sscanf(Arg_input[i].v+6,"%s",range);
        k=strlen(range);
        for(j=0;j<k;j++) {
          if(range[j]==':') range[j]=' ';
        }
        sscanf(range,"%ld %ld",&part.aar_g1[0],&part.aar_g2[0]);
        sda_selector=1;
      }
      else if(!strncmp(Arg_input[i].v,"-gap2:",6)) {
        sscanf(Arg_input[i].v+6,"%s",range);
        k=strlen(range);
        for(j=0;j<k;j++) {
          if(range[j]==':') range[j]=' ';
        }
        sscanf(range,"%ld %ld",&part.aar_g1[1],&part.aar_g2[1]);
        sda_selector=1;
      }
      else if(!strncmp(Arg_input[i].v,"-aa1:",5)) {
        sscanf(Arg_input[i].v+5,"%s",range);
        k=strlen(range);
        for(j=0;j<k;j++) {
          if(range[j]==':') range[j]=' ';
        }
        sscanf(range,"%ld %ld",&part.aar_n1[0],&part.aar_n2[0]);
        sda_selector=1;
      }
      else if(!strncmp(Arg_input[i].v,"-aa2:",5)) {
        sscanf(Arg_input[i].v+5,"%s",range);
        k=strlen(range);
        for(j=0;j<k;j++) {
          if(range[j]==':') range[j]=' ';
        }
        sscanf(range,"%ld %ld",&part.aar_n1[1],&part.aar_n2[1]);
        sda_selector=1;
      }
      else if(!strncmp(Arg_input[i].v,"-er1:",5)) {
        sscanf(Arg_input[i].v+5,"%s",range);
        k=strlen(range);
        n=0;
        n_ok=0;
        n_end=0;
        for(j=0;j<k;j++) {
          if(range[j]==',') {
            range[j]=' ';
            n_ok++;
            n_end=1;
          }
          else if(range[j]==':') {
            range[j]=' ';
            n_ok++;
          }
          if(n_ok > 1 || j == k-1 || n_end == 1) {
            if(n_ok > 1 || (j == k-1 && n_end == 0)) {
              sscanf(range+n,"%s %s",part.er[0].b[part.er[0].n].aa,part.er[0].e[part.er[0].n].aa);
            }
            else {
              sscanf(range+n,"%s",part.er[0].b[part.er[0].n].aa);
              sscanf(range+n,"%s",part.er[0].e[part.er[0].n].aa);
            }
            n_end=0;
            n_ok=0;
            n=j;
            part.er[0].b[part.er[0].n].found=0;
            part.er[0].e[part.er[0].n].found=0;
            part.er[0].n++;
            if(part.er[0].n>100) {
              printf(" ERROR! CHECK the number of -er1 ranges. MAX = 100\n");
              exit(0);
            }
          }
        }
        part.er[0].ok=0;
        sda_selector=1;
      }
      else if(!strncmp(Arg_input[i].v,"-er2:",5)) {
        sscanf(Arg_input[i].v+5,"%s",range);
        k=strlen(range);
        n=0;
        n_ok=0;
        n_end=0;
        for(j=0;j<k;j++) {
          if(range[j]==',') {
            range[j]=' ';
            n_ok++;
            n_end=1;
          }
          else if(range[j]==':') {
            range[j]=' ';
            n_ok++;
          }
          if(n_ok > 1 || j == k-1 || n_end == 1) {
            if(n_ok > 1 || (j == k-1 && n_end == 0)) {
              sscanf(range+n,"%s %s",part.er[1].b[part.er[1].n].aa,part.er[1].e[part.er[1].n].aa);
            }
            else {
              sscanf(range+n,"%s",part.er[1].b[part.er[1].n].aa);
              sscanf(range+n,"%s",part.er[1].e[part.er[1].n].aa);
            }
            n_end=0;
            n_ok=0;
            n=j;
            part.er[1].b[part.er[1].n].found=0;
            part.er[1].e[part.er[1].n].found=0;
            part.er[1].n++;
            if(part.er[1].n>100) {
              printf(" ERROR! CHECK the number of -er2 ranges. MAX = 100\n");
              exit(0);
            }
          }
        }
        part.er[1].ok=0;
        sda_selector=1;
      }
      else if(!strncmp(Arg_input[i].v,"-gdc_sup:",9)) {
        sscanf(Arg_input[i].v+9,"%s",range);
        k=strlen(range);
        n=0;
        n_ok=0;
        for(j=0;j<k;j++) {
          if(range[j]==':' || range[j]==',') {
            range[j]=' ';
            n_ok++;
          }
          if(n_ok > 1 || j == k-1) {
            sscanf(range+n,"%s %s",part.gdc_sup.b[part.gdc_sup.n].aa,part.gdc_sup.e[part.gdc_sup.n].aa);
            n_ok=0;
            n=j;
            part.gdc_sup.b[part.gdc_sup.n].found=0;
            part.gdc_sup.e[part.gdc_sup.n].found=0;
            part.gdc_sup.n++;
            if(part.gdc_sup.n>100) {
              printf(" ERROR! CHECK the number of -gdc_sup ranges. MAX = 100\n");
              exit(0);
            }
          }
        }
        part.gdc_sup.ok=0;
        part.eval=1;
      }
      else if(!strncmp(Arg_input[i].v,"-gdc_set:",9)) {
        sscanf(Arg_input[i].v+9,"%s",range);
        k=strlen(range);
        n=0;
        n_ok=0;
        for(j=0;j<k;j++) {
          if(range[j]==':' || range[j]==',') {
            range[j]=' ';
            n_ok++;
          }
          if(n_ok > 1 || j == k-1) {
            sscanf(range+n,"%s %s",part.gdc_set.b[part.gdc_set.n].aa,part.gdc_set.e[part.gdc_set.n].aa);
            n_ok=0;
            n=j;
            part.gdc_set.b[part.gdc_set.n].found=0;
            part.gdc_set.e[part.gdc_set.n].found=0;
            part.gdc_set.n++;
            if(part.gdc_set.n>100) {
              printf(" ERROR! CHECK the number of -gdc_set ranges. MAX = 100\n");
              exit(0);
            }
          }
        }
        part.gdc_set.ok=0;
        part.eval=1;
      }
      else if(!strncmp(Arg_input[i].v,"-gdc_at:",8) ||
              !strncmp(Arg_input[i].v,"-gdc_aat:",9)) {
        if(part.gdc_at_mcb<0) {
          if(!strncmp(Arg_input[i].v,"-gdc_at:",8)) sscanf(Arg_input[i].v+8,"%s",range);
          else sscanf(Arg_input[i].v+9,"%s",range);
          k=strlen(range);
          n=0;
          n_ok=0;
          for(j=0;j<k;j++) {
            if(range[j]=='.' || range[j]==',') {
              range[j]=' ';
              n_ok++;
            }
            if(n_ok > 1 || j == k-1) {
              sscanf(range+n,"%s %s",part.gdc_at.b[part.gdc_at.n].aa,part.gdc_at.e[part.gdc_at.n].aa);
              n_ok=0;
              n=j;
              if(part.gdc_at.b[part.gdc_at.n].aa[0]=='*') {
                part.gdc_at.b[0].aa[0]='*'; part.gdc_at.b[0].aa[1]='\0';
                if(!strncmp(part.gdc_at.e[part.gdc_at.n].aa,"N\0",2)) { part.gdc_at_mcb=0; }
                if(!strncmp(part.gdc_at.e[part.gdc_at.n].aa,"CA\0",3)) { part.gdc_at_mcb=1; }
                if(!strncmp(part.gdc_at.e[part.gdc_at.n].aa,"C\0",2)) { part.gdc_at_mcb=2; }
                if(!strncmp(part.gdc_at.e[part.gdc_at.n].aa,"O\0",2)) { part.gdc_at_mcb=3; }
                if(!strncmp(part.gdc_at.e[part.gdc_at.n].aa,"CB\0",3)) { part.gdc_at_mcb=4; }
                if(part.gdc_at_mcb<0) {
                  printf(" ERROR! CHECK selected atom: -gdc_at:%s.%s\n",
                         part.gdc_at.b[part.gdc_at.n].aa,part.gdc_at.e[part.gdc_at.n].aa);
                  printf("        Selecting amino-acid '*' only mainchain+CB atoms are allowed: N , CA , C , O , CB\n");
                  exit(0);
                }
                j=k;
                part.gdc_at.n=1;
              }
              part.gdc_at.b[part.gdc_at.n].found=0;
              part.gdc_at.e[part.gdc_at.n].found=0;
              part.gdc_at.n++;
              if(part.gdc_at.n>22) {
                printf(" ERROR! CHECK the number of selected amino-acid atoms (-gdc_at). MAX = 21\n");
                exit(0);
              }
            }
          }
          part.gdc_at.ok=0;
          part.eval=1;
        }
      }
      else if(!strncmp(Arg_input[i].v,"-gdc_eat:",9)) {
        sscanf(Arg_input[i].v+9,"%s",range);
        k=strlen(range);
        n=0;
        n_ok=0;
        for(j=0;j<k;j++) {
          if(range[j]=='.' || range[j]==':' || range[j]==',') {
            range[j]=' ';
            n_ok++;
          }
          if(n_ok > 3 || j == k-1) {
            sscanf(range+n,"%s %s %s %s",part.gdc_eval.b[part.gdc_eval.n].aa,part.gdc_eval.e[part.gdc_eval.n].aa,part.gdc_eval.b[part.gdc_eval.n+1].aa,part.gdc_eval.e[part.gdc_eval.n+1].aa);
            n_ok=0;
            n=j;
            part.gdc_eval.b[part.gdc_eval.n].found=0;
            part.gdc_eval.b[part.gdc_eval.n+1].found=0;
            strcpy(part.gdc_eval.b[part.gdc_eval.n].code3,"XXX\0");
            strcpy(part.gdc_eval.b[part.gdc_eval.n+1].code3,"XXX\0");
            part.gdc_eval.n=part.gdc_eval.n+2;
            if(part.gdc_eval.n>50) {
              printf(" ERROR! CHECK the number of -gdc_eat atom pairs. MAX = 50\n");
              exit(0);
            }
          }
        }
        part.gdc_eval.ok=0;
        part.eval=1;
      }
      else if(!strncmp(Arg_input[i].v,"-atom:",6)) {
        sscanf(Arg_input[i].v+6,"%s",fit);
        sscanf(fit,"%s",pdata.atoms);
        strcat(pdata.atoms,"\0");
        pdata.len_atom=strlen(pdata.atoms);
        for(j=0;j<pdata.len_atom;j++) {
          if(pdata.atoms[j]==',') pdata.atoms[j]='\'';
        }
        pdata.len_atom++;
        pdata.ca_atoms=0;
        if(!strncmp(pdata.atoms,"CB\0",3)) {
          pdata.use_CB=1;
        }
        if(pdata.cb_calc==1) {
          strcpy(pdata.atoms,"BMO\0");
        }
      }
      else if(!strncmp(Arg_input[i].v,"-ch1:",5)) {
        fit[0]='*';
        sscanf(Arg_input[i].v+5,"%s",fit);
        pdata.aa1_ch=fit[0];
        if(isalnum(pdata.aa1_ch)==0) pdata.aa1_ch=' ';
      }
      else if(!strncmp(Arg_input[i].v,"-ch2:",5)) {
        fit[0]='*';
        sscanf(Arg_input[i].v+5,"%s",fit);
        pdata.aa2_ch=fit[0];
        if(isalnum(pdata.aa2_ch)==0) pdata.aa2_ch=' ';
      }
      else if(!strncmp(Arg_input[i].v,"-fit:",5)) {
        sscanf(Arg_input[i].v+5,"%s",fit);
        sscanf(fit,"%1ld",&part.fit_b);
        if(fit[3]==':') {
          sscanf(fit+2,"%1ld",&part.fit_g);
          sscanf(fit+4,"%s",part.fit_r);
        }
        else if(fit[4]==':') {
          sscanf(fit+2,"%2ld",&part.fit_g);
          sscanf(fit+5,"%s",part.fit_r);
        }
        else ok=2;
        if(part.fit_b<0 || part.fit_b>9) part.fit_b=0;
        if(part.fit_g<0 || part.fit_g>99) part.fit_g=0;
        if(part.fit_r[0]=='#') ok=2;
      }
      else if(!strncmp(Arg_input[i].v,"-mol1:",6)) {
        sscanf(Arg_input[i].v+6,"%s",fit);
        sscanf(fit,"%s",pdata.mol1);
        k=strlen(pdata.mol1);
        for(j=0;j<k;j++) {
          if(isalnum(pdata.mol1[j])==0 && pdata.mol1[j]!='_' && pdata.mol1[j]!='\0') {
            j=k; strcpy(pdata.mol1,"*");
          }
        }
        if(k<1) strcpy(pdata.mol1,"*");
      }
      else if(!strncmp(Arg_input[i].v,"-mol2:",6)) {
        sscanf(Arg_input[i].v+6,"%s",fit);
        sscanf(fit,"%s",pdata.mol2);
        k=strlen(pdata.mol2);
        for(j=0;j<k;j++) {
          if(isalnum(pdata.mol2[j])==0 && pdata.mol2[j]!='_' && pdata.mol2[j]!='\0') {
            j=k; strcpy(pdata.mol2,"*");
          }
        }
        if(k<1) strcpy(pdata.mol2,"*");
      }
      else if(!strncmp(Arg_input[i].v,"-1\0",3)) {
        pdata.s_analysis=1;
      }
      else if(!strncmp(Arg_input[i].v,"-2\0",3)) {
        pdata.s_analysis=2;
      }
      else if(!strncmp(Arg_input[i].v,"-3\0",3)) {
        pdata.s_analysis=3;
      }
      else if(!strncmp(Arg_input[i].v,"-4\0",3)) {
        pdata.s_analysis=4;
      }
      else if(!strncmp(Arg_input[i].v,"-5\0",3)) {
        pdata.s_analysis=5;
      }
      else if(!strncmp(Arg_input[i].v,"-d:",3) || !strncmp(Arg_input[i].v,"-d_",3)) {
        sscanf(Arg_input[i].v+3,"%f",&pdata.s_distance);
      }
      else if(!strncmp(Arg_input[i].v,"-lga:",5)) {
        sscanf(Arg_input[i].v+5,"%f",&part.lga_w);
      }
      else if(!strncmp(Arg_input[i].v,"-lw:",4) || !strncmp(Arg_input[i].v,"-lw_",4)) {
        sscanf(Arg_input[i].v+4,"%ld",&part.rw_l);
      }
      else if(!strncmp(Arg_input[i].v,"-lN:",4)) {
        sscanf(Arg_input[i].v+4,"%ld",&part.lN_n);
      }
      else if(!strncmp(Arg_input[i].v,"-ah:",4)) {
        sscanf(Arg_input[i].v+4,"%ld",&pdata.ah_i);
        pdata.ca_atoms=0;
      }
      else if(!strncmp(Arg_input[i].v,"-batch:",7)) {
        batch_mode=1;
        sscanf(Arg_input[i].v+7,"%s",batch_name);
        ok++;
        if((fp = fopen(batch_name,"r"))==NULL) {
          printf("\n# ERROR! opening file %s for read\n\n",batch_name);
          exit(0);
        }
        while(fgets(line,1000,fp)!=NULL) {
          for(n=0;n<200;n++) keyword[n]=0;
          sscanf(line,"%s",keyword);
          if(!strncmp(keyword,"STRUCTURES:\0",12)) {
            sscanf(line,"%s %s",keyword,fname);
          }
          if(!strncmp(keyword,"RUN:\0",5)) {
            batch_number++;
            if(batch_number>=400) {
              printf(" ERROR! CHECK the number of RUN lines in batch file. MAX = 400\n");
              exit(0);
            }
            j=0;
            for(n=4; n<1000; n++) {
              if(line[n]!='\0' && line[n]!='\n' && line[n]!='\t' && line[n]!='\r') {
                Arg_input[batch_number].line[j]=line[n];
                j++;
              }
              if(line[n]=='\n') {
                n=1000;
              }
            }
          }
        }
        fclose(fp);
        Total_batch_number=batch_number+1;
      }
      else if(batch_mode!=1) {
        ok++;
        strcpy(fname,Arg_input[i].v);
      }
      else {
        printf(" ERROR! CHECK the parameters. Unknown parameter: %s\n",Arg_input[i].v);
        exit(0);
      }
    }
    for(j=0;j<2;j++) {
      if(part.aar_g1[j]>part.aar_g2[j]) {
        part.aar_g1[j]=9999;
        part.aar_g2[j]=-9999;
      }
      if(part.aar_n1[j]>=part.aar_n2[j]) {
        part.aar_n1[j]=-9999;
        part.aar_n2[j]=9999;
      }
    }
    if(part.rw_l<1 || part.rw_l>1000) part.rw_l=0;
    if(part.lN_n<1 || part.lN_n>1000) part.lN_n=1;
    if(pdata.ah_i<0 || pdata.ah_i>2) pdata.ah_i=0;
    if(part.lga_w<0.0 || part.lga_w>1.0) part.lga_w=0.75;
    if(pdata.s_distance<0.10 || pdata.s_distance>10.00) pdata.s_distance=5.0;
    if(pdata.cb_v<-5.0 || pdata.cb_v>5.0) pdata.cb_v=0.0;
    if(pdata.m_v<-5.0 || pdata.m_v>5.0) pdata.m_v=0.0;
    if(pdata.o_v<-5.0 || pdata.o_v>5.0) pdata.o_v=0.0;
    if(pdata.sda==1 && pdata.aa1_ch=='*' && pdata.aa2_ch!='*') pdata.aa1_ch=' ';
    if(pdata.sda==1 && pdata.aa1_ch!='*' && pdata.aa2_ch=='*') pdata.aa2_ch=' ';
    if(part.stral==1) part.rw_l=1;
    if(part.eval>0) part.all_rmsd=1;
    if(part.all_rmsd>0) part.eval=1;
/*
      printf(" TEST: aa1_n1 = %ld , aa1_n2 = %ld \n",pdata.aa1_n1,pdata.aa1_n2);
      printf(" TEST: aa2_n1 = %ld , aa2_n2 = %ld \n",pdata.aa2_n1,pdata.aa2_n2);
      exit(0);
*/
    if(ok!=1) {
      printf(" ERROR! CHECK the parameters!\n Parameters line: ");
      for(i=1;i<Arg_input_c;i++) { printf("%s ",Arg_input[i].v); }
      printf("\n\n Usage: lga [<options>] file_name");
      printf("\n HELP:  lga -h\n");
      exit(0);
    }
    pdata.parameters_n=Arg_input_c;
    sprintf(pdata.parameters[0].line,"LGA parameters: ");
    for(i=1;i<Arg_input_c;i++) {
      sprintf(pdata.parameters[i].line,"%s ",Arg_input[i].v); 
    }
    part.dist_cutoff=pdata.s_distance;
    part.gdt_cutoff[0]=part.dist_cutoff;    
    
/*  Read longo pdata
-------------------------------------------*/
    strcpy(pdata.mname,fname);
    if(batch_mode!=2) {
      strcpy(pdata.fname_in,dir_in);
      strcat(pdata.fname_in,fname);
      strcpy(part.fname_in,pdata.fname_in);
    }
      
    strcpy(mname,fname);
    if(strncmp(pdata.mol1,"*",1) || strncmp(pdata.mol2,"*",1)){
      strcpy(mname,pdata.mol1);
      strcat(mname,".");
      strcat(mname,pdata.mol2);
      if(!strncmp(pdata.mol1,"*",1)) {
        strcpy(mname,pdata.mol2);
      }
      if(!strncmp(pdata.mol2,"*",1)) {
        strcpy(mname,pdata.mol1);
      }
    }
    strcpy(pdata.fname_out,dir_out);
    strcat(pdata.fname_out,mname);
    strcpy(part.mname_lga,mname);
    strcpy(part.fname_lga,dir_out);
    strcat(part.fname_lga,mname);
    
    if(pdata.s_analysis==5) {
      part.sia=1-pdata.sda;
      pdata.sda=0;
    }
    if(batch_mode!=2) {       // batch (0 , 1) reading structures
      if(part.full_print==1) printf("# Reading molecules from the file: %s \n",fname);
      if((fp = fopen(pdata.fname_in,"r"))==NULL) {
        printf("\n# ERROR! opening file %s for read\n\n",fname);
        exit(0);
      }
      if(ignore_aa==0) {
        read_aamol2(aa_mol2, fp);
        rewind(fp);
      }
      if(pdata.s_analysis==0 || ignore_aa==1) {
        aa_mol2[0].n_aa=0;
        aa_mol2[1].n_aa=0;
      }
      read_pdb_2sets(aa_mol2, &pdata, &part, fp);
      pdata.molinput[0].n_aa_all=pdata.molinput[0].n_aa;
      pdata.molinput[1].n_aa_all=pdata.molinput[1].n_aa;
      if(pdata.sda==1) {
        pdata.ignore_errors=pdata.ignore_errors+2;
        n_atoms=select_atoms(ca_atoms, &pdata, &part);
        pdata.sda=2;
        sda_list(aa_mol2, ca_atoms, &pdata);
        rewind(fp);
        read_pdb_2sets(aa_mol2, &pdata, &part, fp);
      }
      fclose(fp);
      if(part.full_print==1) printf("\n# Reading molecules                         (DONE)\n");
    }

    strcpy(part.atoms,pdata.atoms);
    if(batch_mode!=1) {       // batch (0 , 2) processing structures
      if(part.full_print==1) printf("\n# Selecting atoms for processing .......... \n");
      n_atoms=select_atoms(ca_atoms, &pdata, &part);
      if(part.full_print==1) printf("\n# Selecting atoms                           (DONE)\n");

      printf("\n# Molecule1: number of %s atoms %4ld (%5ld),  selected %4ld , name %s",
              pdata.atoms, pdata.molinput[0].n_aa_all, pdata.molinput[0].n_atoms_all,ca_atoms[0][0],
              pdata.molecule[0].molecule_name);
      printf("\n# Molecule2: number of %s atoms %4ld (%5ld),  selected %4ld , name %s",
              pdata.atoms, pdata.molinput[1].n_aa_all, pdata.molinput[1].n_atoms_all,ca_atoms[1][0],
              pdata.molecule[1].molecule_name);

      if(pdata.sda>0 && pdata.s_analysis<5 && sda_selector==0) {
        part.m_n_aa=pdata.molinput[0].n_aa_all;
        part.t_n_aa=pdata.molinput[1].n_aa_all;
      }
      else {
        part.m_n_aa=ca_atoms[0][0];
        part.t_n_aa=ca_atoms[1][0];
      }
      
      printf("\n# PARAMETERS: ");
      for(i=1;i<pdata.parameters_n;i++) {
        printf("%s ",pdata.parameters[i].line);
      }
//      printf("\n");

      if(pdata.s_analysis==4) {
        printf("\n# Search for Atom-Atom correspondence");
      }
      else {
        printf("\n# FIXED Atom-Atom correspondence");
      }
      if(pdata.s_analysis==0) {
        print_aamol2(ca_atoms, &pdata);
      }
      else if(pdata.s_analysis==1) {
        printf("\n# Standard RMSD calculation \n");
        ca_atoms[0][0]=n_atoms;
        ca_atoms[1][0]=n_atoms;
        part.rmsd=rmsd(n_atoms, ca_atoms, pdata.molecule, r, v);
      }
      else if(pdata.s_analysis==2) {
        printf("\n# RMSD calculation using Iterative Superposition Procedure \n");
        part.isp=1;
        ca_atoms[0][0]=n_atoms;
        ca_atoms[1][0]=n_atoms;
        part.rmsd_isp=rmsd_isp(n_atoms, ca_atoms, pdata.molecule, r, v, &part);
      }
      else if(pdata.s_analysis==3) {
        printf("\n# GDT and LCS analysis \n");
        ca_atoms[0][0]=n_atoms;
        ca_atoms[1][0]=n_atoms;
        lcs_gdt_analysis(n_atoms, ca_atoms, pdata.molecule, r, v, &part);
      }
      else if(pdata.s_analysis==4) {
        printf("\n# Structure alignment analysis \n");
        alignment_best_match(ca_atoms, pdata.molecule, r, v, &part);
      }
      else if(pdata.s_analysis==5) {
        if(part.sia==0) {
          printf("\n# Structure best fit analysis (S => S-gap-S) \n");
        }
        else {
          printf("\n# Structure best fit analysis (S-gap-S => S) \n");
        }
        if(part.fit_g>0 && part.fit_n==0) {
          if(part.sia==0) {
            printf(" ERROR! CHECK the numbering of residues in Molecule2 and the residue\n");
            printf("        number %s after which GAP %2ld residues long follows.\n",part.fit_r,part.fit_g);
          }
          else {
            printf(" ERROR! CHECK the numbering of residues in Molecule1 and the residue\n");
            printf("        number %s after which GAP %2ld residues long follows.\n",part.fit_r,part.fit_g);
          }
          exit(0);
        }
        if(ca_atoms[0][0]<ca_atoms[1][0]+part.fit_g && part.sia==0) {
          printf(" ERROR! CHECK the number of residues:\n");
          printf("        Molecule1 >= Molecule2 + GAP \n");
          printf("          %4ld    <    %4ld    + %2ld   ERROR! \n",ca_atoms[0][0],ca_atoms[1][0],part.fit_g);
          exit(0);
        }
        if(ca_atoms[1][0]<ca_atoms[0][0]+part.fit_g && part.sia==1) {
          printf(" ERROR! CHECK the number of residues:\n");
          printf("        Molecule1 + GAP <= Molecule2 \n");
          printf("            %4ld  + %2ld  >   %4ld   ERROR! \n",ca_atoms[0][0],part.fit_g,ca_atoms[1][0]);
          exit(0);
        }
        best_fit(ca_atoms, pdata.molecule, r, v, &part, &pdata);
      }
      if(part.full_print==1) printf("\n# Processing data                           (DONE)\n");

      if(pdata.s_analysis>0 && pdata.s_analysis<5)
        write_output(&part, &pdata, ca_atoms, r, v, pdata.fname_out);
    }
    if(batch_mode==1) batch_mode=2; // batch: switch to runs after reading structures
    if(batch_mode==2 && batch_number>0) {
      n=0;
      word=strtok(Arg_input[Total_batch_number - batch_number].line," ");
      while(word!=NULL) {
        n++;
        strcpy(Arg_input[n].v,word);
        word=strtok(NULL," ");
      }
      Arg_input_c=n+1;
//      printf("\n# test: %ld # %s\n",Arg_input_c,Arg_input[Total_batch_number - batch_number].line);
    }

  }
  
  printf("\n# END of job\n");
  /* Free workspace */
  lga_ws_free(&part.ws);
  exit(0);
}
