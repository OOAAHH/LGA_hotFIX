#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "tools.h"

/*-----------------------------------------------------------
/
/   check_aa_atoms - checks the atoms for each amino acid
/
/------------------------------------------------------------*/
long check_aa_atoms(char code1, char* atom, long *n_def)
{
  long id_atom;
  
  id_atom=-1;
  *n_def=0;

  if(code1=='A') {                               /*   ALA   */
    if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"CB\0",3)) {
      id_atom=4;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=5;
  } else
  if(code1=='V') {                               /*   VAL   */
    if(!strncmp(atom,"CG1\0",4)) {
      id_atom=5;
    }
    else if(!strncmp(atom,"CG2\0",4)) {
      id_atom=6;
    }
    else if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"CB\0",3)) {
      id_atom=4;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=7;
  } else 
  if(code1=='L') {                               /*   LEU   */
    if(!strncmp(atom,"CD1\0",4)) {
      id_atom=6;
    }
    else if(!strncmp(atom,"CD2\0",4)) {
      id_atom=7;
    }
    else if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"CB\0",3)) {
      id_atom=4;
    }
    else if(!strncmp(atom,"CG\0",3)) {
      id_atom=5;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=8;
  } else 
  if(code1=='I') {                               /*   ILE   */
    if(!strncmp(atom,"CD1\0",4) || !strncmp(atom,"CD\0",3)) {  /* PDB standard is CD1 while CHARMM uses CD */
      id_atom=7;
    }
    else if(!strncmp(atom,"CG1\0",4)) {
      id_atom=5;
    }
    else if(!strncmp(atom,"CG2\0",4)) {
      id_atom=6;
    }
    else if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"CB\0",3)) {
      id_atom=4;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=8;
  } else 
  if(code1=='P') {                               /*   PRO   */
    if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"CB\0",3)) {
      id_atom=4;
    }
    else if(!strncmp(atom,"CG\0",3)) {
      id_atom=5;
    }
    else if(!strncmp(atom,"CD\0",3)) {
      id_atom=6;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=7;
  } else 
  if(code1=='M') {                               /*   MET   */
    if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"CB\0",3)) {
      id_atom=4;
    }
    else if(!strncmp(atom,"CG\0",3)) {
      id_atom=5;
    }
    else if(!strncmp(atom,"CE\0",3)) {
      id_atom=7;
    }
    else if(!strncmp(atom,"SD\0",3)) {
      id_atom=6;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=8;
  } else 
  if(code1=='F') {                               /*   PHE   */
    if(!strncmp(atom,"CD1\0",4)) {
      id_atom=6;
    }
    else if(!strncmp(atom,"CD2\0",4)) {
      id_atom=7;
    }
    else if(!strncmp(atom,"CE1\0",4)) {
      id_atom=8;
    }
    else if(!strncmp(atom,"CE2\0",4)) {
      id_atom=9;
    }
    else if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"CB\0",3)) {
      id_atom=4;
    }
    else if(!strncmp(atom,"CG\0",3)) {
      id_atom=5;
    }
    else if(!strncmp(atom,"CZ\0",3)) {
      id_atom=10;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=11;
  } else 
  if(code1=='W') {                               /*   TRP   */
    if(!strncmp(atom,"CD1\0",4)) {
      id_atom=6;
    }
    else if(!strncmp(atom,"CD2\0",4)) {
      id_atom=7;
    }
    else if(!strncmp(atom,"CE2\0",4)) {
      id_atom=9;
    }
    else if(!strncmp(atom,"CE3\0",4)) {
      id_atom=10;
    }
    else if(!strncmp(atom,"CH2\0",4)) {
      id_atom=13;
    }
    else if(!strncmp(atom,"CZ2\0",4)) {
      id_atom=11;
    }
    else if(!strncmp(atom,"CZ3\0",4)) {
      id_atom=12;
    }
    else if(!strncmp(atom,"NE1\0",4)) {
      id_atom=8;
    }
    else if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"CB\0",3)) {
      id_atom=4;
    }
    else if(!strncmp(atom,"CG\0",3)) {
      id_atom=5;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=14;
  } else 
  if(code1=='G') {                               /*   GLY   */
    if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=4;
  } else 
  if(code1=='S') {                               /*   SER   */
    if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"CB\0",3)) {
      id_atom=4;
    }
    else if(!strncmp(atom,"OG\0",3)) {
      id_atom=5;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=6;
  } else 
  if(code1=='T') {                               /*   THR   */
    if(!strncmp(atom,"CG2\0",4)) {
      id_atom=6;
    }
    else if(!strncmp(atom,"OG1\0",4)) {
      id_atom=5;
    }
    else if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"CB\0",3)) {
      id_atom=4;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=7;
  } else 
  if(code1=='C') {                               /*   CYS   */
    if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"CB\0",3)) {
      id_atom=4;
    }
    else if(!strncmp(atom,"SG\0",3)) {
      id_atom=5;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=6;
  } else 
  if(code1=='Y') {                               /*   TYR   */
    if(!strncmp(atom,"CD1\0",4)) {
      id_atom=6;
    }
    else if(!strncmp(atom,"CD2\0",4)) {
      id_atom=7;
    }
    else if(!strncmp(atom,"CE1\0",4)) {
      id_atom=8;
    }
    else if(!strncmp(atom,"CE2\0",4)) {
      id_atom=9;
    }
    else if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"CB\0",3)) {
      id_atom=4;
    }
    else if(!strncmp(atom,"CG\0",3)) {
      id_atom=5;
    }
    else if(!strncmp(atom,"CZ\0",3)) {
      id_atom=10;
    }
    else if(!strncmp(atom,"OH\0",3)) {
      id_atom=11;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=12;
  } else 
  if(code1=='N') {                               /*   ASN   */
    if(!strncmp(atom,"ND2\0",4)) {
      id_atom=7;
    }
    else if(!strncmp(atom,"OD1\0",4)) {
      id_atom=6;
    }
    else if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"CB\0",3)) {
      id_atom=4;
    }
    else if(!strncmp(atom,"CG\0",3)) {
      id_atom=5;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=8;
  } else 
  if(code1=='Q') {                               /*   GLN   */
    if(!strncmp(atom,"OE1\0",4)) {
      id_atom=7;
    }
    else if(!strncmp(atom,"NE2\0",4)) {
      id_atom=8;
    }
    else if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"CB\0",3)) {
      id_atom=4;
    }
    else if(!strncmp(atom,"CD\0",3)) {
      id_atom=6;
    }
    else if(!strncmp(atom,"CG\0",3)) {
      id_atom=5;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=9;
  } else 
  if(code1=='D') {                               /*   ASP   */
    if(!strncmp(atom,"OD1\0",4)) {
      id_atom=6;
    }
    else if(!strncmp(atom,"OD2\0",4)) {
      id_atom=7;
    }
    else if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"CB\0",3)) {
      id_atom=4;
    }
    else if(!strncmp(atom,"CG\0",3)) {
      id_atom=5;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=8;
  } else 
  if(code1=='E') {                               /*   GLU   */
    if(!strncmp(atom,"OE1\0",4)) {
      id_atom=7;
    }
    else if(!strncmp(atom,"OE2\0",4)) {
      id_atom=8;
    }
    else if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"CB\0",3)) {
      id_atom=4;
    }
    else if(!strncmp(atom,"CD\0",3)) {
      id_atom=6;
    }
    else if(!strncmp(atom,"CG\0",3)) {
      id_atom=5;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=9;
  } else 
  if(code1=='K') {                               /*   LYS   */
    if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"CB\0",3)) {
      id_atom=4;
    }
    else if(!strncmp(atom,"CD\0",3)) {
      id_atom=6;
    }
    else if(!strncmp(atom,"CE\0",3)) {
      id_atom=7;
    }
    else if(!strncmp(atom,"CG\0",3)) {
      id_atom=5;
    }
    else if(!strncmp(atom,"NZ\0",3)) {
      id_atom=8;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=9;
  } else 
  if(code1=='R') {                               /*   ARG   */
    if(!strncmp(atom,"NH1\0",4)) {
      id_atom=9;
    }
    else if(!strncmp(atom,"NH2\0",4)) {
      id_atom=10;
    }
    else if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"CB\0",3)) {
      id_atom=4;
    }
    else if(!strncmp(atom,"CD\0",3)) {
      id_atom=6;
    }
    else if(!strncmp(atom,"CG\0",3)) {
      id_atom=5;
    }
    else if(!strncmp(atom,"CZ\0",3)) {
      id_atom=8;
    }
    else if(!strncmp(atom,"NE\0",3)) {
      id_atom=7;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=11;
  } else 
  if(code1=='H') {                               /*   HIS   */
    if(!strncmp(atom,"CD2\0",4)) {
      id_atom=7;
    }
    else if(!strncmp(atom,"CE1\0",4)) {
      id_atom=8;
    }
    else if(!strncmp(atom,"ND1\0",4)) {
      id_atom=6;
    }
    else if(!strncmp(atom,"NE2\0",4)) {
      id_atom=9;
    }
    else if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"CB\0",3)) {
      id_atom=4;
    }
    else if(!strncmp(atom,"CG\0",3)) {
      id_atom=5;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=10;
  }
  if(code1=='X') {                               /*   N CA C O CB  */
    if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    else if(!strncmp(atom,"CB\0",3)) {
      id_atom=4;
    }
    *n_def=5;
  }
  if(code1=='#') {                               /* Unknown amino-acid which has at least mainchain */
    if(!strncmp(atom,"CA\0",3)) {
      id_atom=1;
    }
    else if(!strncmp(atom,"C\0",2)) {
      id_atom=2;
    }
    else if(!strncmp(atom,"N\0",2)) {
      id_atom=0;
    }
    else if(!strncmp(atom,"O\0",2)) {
      id_atom=3;
    }
    *n_def=4;
  }
  
  if(!strncmp(atom,"OXT\0",4)) {                /* N-terminal oxygen */
    id_atom=14;
    *n_def=*n_def+1;
  }

  return id_atom;
}

/*-----------------------------------------------------------
/
/   buger - checks floating polongs numbers
/
/------------------------------------------------------------*/
long buger(char* sub, char* line, long b, long e)
{
  long j, bug, bugd, bugs, gap, gapl, gapr;

  bug=2;
  bugd=0;
  bugs=0;
  for(j=0;j<10;j++) {
    if(isdigit(sub[j])!=0) {bug=0;}
    else if(sub[j]=='-' || sub[j]=='+') {if(bug==0) {bugs=2;} else {bugs++;}}
    else if(sub[j]=='.') {bug=0; bugd++;}
    else if(sub[j]==' ') {}
    else if(sub[j]=='\0') {}
    else if(sub[j]=='\r') {}
    else if(sub[j]=='\n') {}
    else if(sub[j]=='\t') {}
    else {bugd=2;}
  }
  gap=0;
  gapl=0;
  gapr=0;
  for(j=b;j<e;j++) {
    if(line[j]!=' ' && line[j+1]==' ') {
      gap++;
      gapr=1;
    }
    if(line[j]==' ' && line[j+1]!=' ') {
      gap++;
      if(gapr==0) gapl=1;
    }
  }
  if(bugd>1 || bugs>1 || gap!=gapl+gapr || (gap==0 && line[e-2]==' ')) {bug=2;}

  return bug;
}
