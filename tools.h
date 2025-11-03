/*------------------------------------------------------------
/
/ tools.h 
/
/-----------------------------------------------------------*/

#include <stdio.h>  /* for FILE* used in function prototypes */

#define MAXATOMS  95001 // 65536
#define MAXRES    9501  // 4501  // 8001 (4001 allows comparison of 2000 residue long structures)
#define MAXGAP2   65536 // 32761 = MAXGAP**2
#define MAXGAP    256   // 181
#define MAXITR    121   // 121

typedef struct { // Check copy_atom procedure
  float       x;
  float       y;
  float       z;
} vector;        // Check copy_atom procedure

typedef struct {
  vector R;
  long  id_aa;   /* consecutive number of aa where atom belongs to       */
  long  id_atom; /* counts per atom type for current amino acid          */
                  /*      0  1 2 3  4  5   6   7   8   9  10  11  12   13 */
                  /* A:   N CA C O CB                                     */
                  /* V:   N CA C O CB CG1 CG2                             */
                  /* L:   N CA C O CB CG  CD1 CD2                         */
                  /* I:   N CA C O CB CG1 CG2 CD1                         */
                  /* P:   N CA C O CB CG  CD                              */
                  /* M:   N CA C O CB CG  SD  CE                          */
                  /* F:   N CA C O CB CG  CD1 CD2 CE1 CE2 CZ              */
                  /* W:   N CA C O CB CG  CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2 */
                  /* G:   N CA C O                                        */
                  /* S:   N CA C O CB OG                                  */
                  /* T:   N CA C O CB OG1 CG2                             */
                  /* C:   N CA C O CB SG                                  */
                  /* Y:   N CA C O CB CG  CD1 CD2 CE1 CE2 CZ  OH          */
                  /* N:   N CA C O CB CG  OD1 ND2                         */
                  /* Q:   N CA C O CB CG  CD  OE1 NE2                     */
                  /* D:   N CA C O CB CG  OD1 OD2                         */
                  /* E:   N CA C O CB CG  CD  OE1 OE2                     */
                  /* K:   N CA C O CB CG  CD  CE  NZ                      */
                  /* R:   N CA C O CB CG  CD  NE  CZ  NH1 NH2             */
                  /* H:   N CA C O CB CG  ND1 CD2 CE1 NE2                 */
                  /* X:   N CA C O CB                                     */
                  /* #:   N CA C O                                        */
                  /* *: OXT - 14          terminal oxygen                 */
} atom_coords;

typedef struct { // Check copy_atom procedure
  long         serial;
  char        name[10];
  char        alt_loc[10];
  char        res_name[10];
  char        chain_id[10];
  long         res_seq;
  char        i_code[10];
  vector      R;
  float       occupancy;
  float       temp_factor;
  char        seg_id[10];
  char        element[10];
  char        charge[10];
  char        res_i[10];
  long         res_n_local;
} atom_struct;   // Check copy_atom procedure

typedef struct {
  long         n_atoms_all;
  long         n_atoms;
  long         n_aa_all;
  long         n_aa;
  long         n_chains;
  char        chain_name[200];
  char        molecule_name[200];
  long       nb_ca[MAXRES+1];
  atom_struct atom[MAXRES+1];   
} pdb_struct;

typedef struct {
  char        code1;          /* aa name 1 letter code */
  char        code3[10];      /* aa name 3 letter code */
  char        aa[11];
  long       at_beg;         /* aa beg atom number */
  long       at_end;         /* aa end atom number */
  long       at_eq_ca;       /* ca atom number     */
  long       at_eq_cb;       /* cb atom number     */
  long       at_eq_mc_beg;   /* aa mc beg atom number */
  long       at_eq_mc_end;   /* aa mc end atom number */
  long       at_eq_al_beg;   /* aa all beg atom number */
  long       at_eq_al_end;   /* aa all end atom number */
  long       n_def;          /* standard number of atoms in aa */
  long       found;
} aa_num;

typedef struct {
  long       n_aa;
  long       n_atom;
  aa_num      equiv_aa[MAXRES+1];
  long       equiv_atom[MAXATOMS+1];
  long       equiv_atom_ca[MAXATOMS+1];
  long       equiv_atom_mc[MAXATOMS+1];
  atom_coords coord[MAXATOMS+1];
} check_mol2;

typedef struct {
  char        line[MAXRES+1];
} stext;                       

typedef struct {
  long         ok;               /* exact range                    */
  long         n;                /* exact range                    */
  aa_num      b[101];           /* exact range                    */
  aa_num      e[101];           /* exact range                    */
} aa_range;

typedef struct {
  char        mname[200];
  char        mol1[200];
  char        mol2[200];
  char        fname_in[200];
  char        fname_out[200];
  long         ignore_errors;    /* ignore errors                       */
  char        aa1_ch;           /* chain id                            */
  char        aa2_ch;           /* chain id                            */
  long         one_output;       /* output: only molecule 1             */
  char        atoms[10];        /* atoms used for calculations         */
  long         ca_atoms;         /* CA atoms used for calculations      */
  long         len_atom;         /* length of the atoms used            */
  long         ah_i;             /* ATOM or HETATM records              */
  long         n_atoms;          /* common number of atoms              */
  long         s_analysis;       /* selection of the type of analysis   */
  long         sda;              /* sequence dependent analysis         */
  long         use_CB;           /* use existing CB                     */
  long         cb_calc;          /* cb calculations                     */
  float       cb_v;             /* cb calculations, value of ca-cb     */
  float       m_v;              /* mean calculations, value of m-cbv   */
  float       o_v;              /* oxygen-mean-cbv (o-m-cbv) value     */
  float       s_distance;       /* selection of the distance for ISP   */
  long         parameters_n;
  stext       parameters[50];
  pdb_struct  molinput[2];            
  pdb_struct  molecule[2];            
} pdata_struct;                       
                                      
typedef struct {
  long       equiv[MAXRES+1];     
  long       align[MAXRES+1];     
  float       rmsd[MAXRES+1];     
} part_ear;                       

/* Workspace for heap-allocated large buffers to avoid stack usage */
typedef struct lga_workspace {
  /* Common 2 x (MAXRES+1) buffers */
  long  *calphas2;   /* used in tools_3.c, tools_4.c */
  long  *l_atoms2;   /* used in tools_3.c */
  long  *out2;       /* used in tools_3.c */

  /* Index/helper arrays */
  long  *ind3;       /* 3 x (MAXRES+1), used in tools_3.c */
  long  *tmp;        /* (MAXRES+1), used in tools_3.c */
  long  *rwa1;       /* (MAXRES+1), used in tools_4.c */
  long  *rwa2;       /* (MAXRES+1), used in tools_4.c */

  /* Geometry buffers: (MAXRES+1) x 3 */
  float *xa3;        /* used in tools_3.c */
  float *xb3;        /* used in tools_3.c */

  /* best_fit support buffers: each 2 x (MAXRES+1) */
  long  *sup;        /* used in tools_4.c */
  long  *sup1;       /* used in tools_4.c */
  long  *sup2;       /* used in tools_4.c */

  /* Large check_mol2 buffer used by print_lga/check_all_atoms */
  check_mol2 *mol;   /* array of 2 elements */

  /* STRAL lines buffer: 200 lines of stext */
  stext *straline;   /* array of 200 elements */
} lga_workspace;

typedef struct {
  long         n_aa;
  aa_num      num[MAXRES+1];
} mol2_aa;                       

typedef struct {
  long         n;
  float       seq_id;
  float       rms_local;
  float       lcs_gdt_ts;
  float       q_norm;
  stext       s[7];   // [5-6] reserved for REMARK in PDB output
} sline;                       

typedef struct {                      
  part_ear    ear[11];
  sline       summary;
  stext       printlga[MAXRES+1];
  long         n_printlga;
  char        atoms[10];          /* atoms used for calculations (see pdata) */
  char        fit_r[200];
  long         error[2];
  long         fit_b;
  long         fit_n;
  long         fit_g;
  long         fit_si;
  long         accuracy_opt;
  long         accuracy_gdt_n;
  float       accuracy_gdt_step;
  long         accuracy_lga_n;
  float       accuracy_lga_step; 
  aa_num      resrangaa[8];       /* residue ranges are reported in SUMMARY line */
  long         resranges;          /* residue ranges are reported in SUMMARY line */
  long         all_rmsd;           /* rmsd are calculated on CA, MC and ALL atoms */
  long         eval;
  long         gdc_bin;            /* number of bins: <0.5, <1.0, ... <10.0       */
  long         gdc_ref;            /* define frame of reference (0, 1, 2) for GDC */
  aa_range    gdc_eval;           /* exact aanumber.atoms for evaluation         */
  aa_range    gdc_at;             /* exact aansme.atoms for evaluation           */
  long         gdc_at_mcb;         /* exact any.mcb_atoms for evaluation          */
  aa_range    gdc_sup;            /* exact range for superposition               */
  aa_range    gdc_set;            /* exact range for evaluation                  */
  aa_range    er[2];              /* exact range                         */
  long         aar_g1[2];          /* range          aar_g1 - aar_g2      */
  long         aar_g2[2];          /* range          aar_g1 - aar_g2      */
  long         aar_n1[2];          /* range          aar_n1 - aar_n2      */
  long         aar_n2[2];          /* range          aar_n1 - aar_n2      */
  long         sia;
  long         best_ind[11];
  long         rw_l;
  long         lN_n;
  float       rw_rms[MAXRES+1];
  float       rmsd_align[11];     /* RMSD_ALIGN calculated on alignment  */
  char        fname_lga[200];
  char        fname_in[200];
  char        mname_lga[200];
  aa_num      model_index[101];
  long         lcs_gdt_print;
  long         stral;              /* STRAL on - 1 , off - 0              */
  float       stral_r;            /* STRAL cutoff for local RMSD         */
  long         isp;                /* ISP on - 1 , off - 0                */
  float       rmsd;               /* RMSD value calculated on all atoms  */
  float       rmsd_isp;           /* RMSD_ISP calculated on all atoms    */
  long         m_n_aa;
  long         t_n_aa;
  long         full_print;
  long         check;
  long         isp_iter;
  long         isp_iter_start;
  long         gdt;                /* reports: 0 - rmsd (std), 1 - gdt (fit) superposition */
  long         gdc;                /* reports: 0 - LGA rmsd (std), 1 - GDC superposition   */
  long         lga_m;              /* report LGA_M (maximum LGA_S value): no - 0 , yes - 1 */
  long         swap;               /* consider swapping: no - 0 , yes - 1                  */
  float       lga_w;              /* weight between GDT and LCS: default 0.75             */
  float       isp_dist;
  float       dist_cutoff;
  float       rms_lcs_cutoff;         
  float       isp_cutoff[11];     /* cutoffs for ISP procedure           */
  float       gdt_cutoff[21];     /* distance cutoffs for GDT procedure  */
  long       opt_sup[21][2][MAXRES+1];
  long       opt_sup1[21][MAXRES+1];
  long       opt_sup2[21][MAXRES+1];
  long       opt_set_nb[21][MAXRES+1];
  long       opt_lcs[MAXRES+1];
  float       opt_sup_rms[21];
  float       opt_all_rms[21];
  float       opt_r[21][3][3];
  float       opt_v[21][3];
  float       opt_r1[21][3][3];
  float       opt_v1[21][3];
  float       opt_r2[21][3][3];
  float       opt_v2[21][3];
  float       gdt_ts;
  float       lcs_ts;
  long       pairs[2][MAXGAP2+1];
  
  /* Heap workspace pointer (allocated once per process) */
  lga_workspace *ws;
} part_struct;

void clean_pdb(pdb_struct *);
void clean_pdata(pdata_struct *);
void clean_part(part_struct *);
int  lga_ws_init(lga_workspace **);
void lga_ws_free(lga_workspace **);
void read_aamol2(mol2_aa[2], FILE*);
void read_pdb_2sets(mol2_aa[2], pdata_struct *, part_struct *, FILE*);
void C_321(char*, char*);
long select_atoms(long[2][MAXRES+1], pdata_struct *, part_struct *);
void copy_atom(atom_struct *, atom_struct *);
void print_aamol2(long[2][MAXRES+1], pdata_struct *);
float rmsd(long, long[2][MAXRES+1], pdb_struct[2],
           float[3][3], float[3]);
float rmsd_any(long, long[MAXATOMS+1], long[MAXATOMS+1], 
           atom_coords[MAXATOMS+1], atom_coords[MAXATOMS+1],
           float[3][3], float[3]);
double qkfit(double[4][4], double[4][4]);
void eigen(double[7], double[4][4], long[4]);
void esort(double[7], double[4][4], long[4]);
float rmsd_isp(long, long[2][MAXRES+1], pdb_struct[2],
               float[3][3], float[3], part_struct *);
float rmsd_dist(long, long[2][MAXRES+1], pdb_struct[2],
                float[3][3], float[3], part_struct *);
void lcs_gdt_analysis(long, long[2][MAXRES+1], pdb_struct[2],
                      float[3][3], float[3], part_struct *);
float lcs_run(long, long[2][MAXRES+1], pdb_struct[2],
              float[3][3], float[3], part_struct *);
void gdt_run(long, long[2][MAXRES+1], pdb_struct[2],
             float[3][3], float[3], part_struct *);
void write_output(part_struct *, pdata_struct *,
                  long[2][MAXRES+1], float[3][3], float[3], char*);
void euler_angles(float[3][3], float[3][4]);
long align_search(long[2], long[2][MAXRES+1], pdb_struct[2],
                 float[3][3], float[3], part_struct *, part_ear *);
long alignment_run(long, long[2][MAXRES+1], long[2][MAXRES+1], pdb_struct[2],
                  float[3][3], float[3], part_struct *);
void alignment_best_match(long[2][MAXRES+1], pdb_struct[2],
                          float[3][3], float[3], part_struct *);
void best_fit(long[2][MAXRES+1], pdb_struct[2],
              float[3][3], float[3], part_struct *, pdata_struct *);
long check_align_lga(long[2][MAXRES+1], long[2][MAXRES+1],
                    pdb_struct[2], float[3][3], float[3], part_struct *, long);
long print_lga(long[2][MAXRES+1], long[2][MAXRES+1], long[2][MAXRES+1],
              pdb_struct[2], float[3][3], float[3], part_struct *, long, long);
void rms_window(long[2][MAXRES+1], long[2][MAXRES+1], long[2][MAXRES+1],
                pdb_struct[2], part_struct *, long, long);
float gdc_dist(long, long, 
           atom_coords[MAXATOMS+1], atom_coords[MAXATOMS+1],
           float[3][3], float[3]);
float distance_calc(long[2][MAXRES+1],
                    pdb_struct[2], float[3][3], float[3], long, long);
long align_select(part_ear *, part_ear *, part_ear *);
void align_collect(part_ear *, part_ear *, part_ear *,
                   long[2][MAXRES+1], pdb_struct[2], 
                   float[3][3], float[3], float);
void align_final(part_ear *, long[2][MAXRES+1], pdb_struct[2],
                 float[3][3], float[3], float);
void sda_list(mol2_aa[2], long[2][MAXRES+1], pdata_struct *);
long lchk(long);
void calc_cb(float, float, float, pdb_struct[2], long);
float check_all_atoms(pdb_struct[2], part_struct *, check_mol2[2],
                      float[3][3], float[3]);
long check_aa_atoms(char, char*, long*);
float rot_sqr(vector, vector, float[3][3], float[3]);
long buger(char*, char*, long, long);
