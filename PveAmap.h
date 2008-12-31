#include  <bicpl.h>

struct dartel_prm {
  int rform;         // regularization form: 0 - linear elastic energy; 1 - membrane energy; 2 - bending energy
  double rparam[6];  // regularization parameters
  double lmreg;      // LM regularization
  int cycles;        // number of cycles for full multi grid (FMG)
  int its;           // Number of relaxation iterations in each multigrid cycle
  int k;             // time steps for solving the PDE
  int code;          // objective function: 0 - sum of squares; 1 - symmetric sum of squares; 2 - multinomial
};

int n_pure_classes = 3;
int n_classes = 5;
int Nflips = 50;
int Niters = 200;
int subsample = 8;
int iters_nu = 40;
int pve = 2;
int use_watershed = FALSE;
int correct_nu = TRUE;
int write_fuzzy = FALSE;
int write_nu = FALSE;
int write_label = TRUE;
int warp_priors = TRUE;
double thresh_brainmask = 0.05;
double thresh_kmeans = 0.5;