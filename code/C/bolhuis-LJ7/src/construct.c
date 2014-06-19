#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "path.h"

Replica trial_replica, *replica[MAXREPLICA];
Sys sys;
Slice *slice, *trial, *fwtrial, *bwtrial;
Path path;
Langevin langevin;
Stats nulstat = { 0 };
vector nulvec = { 0 };
Scr scr, nulscr = { 0 };
Bivar bivar;
State state[MAXSTATES];

int input(char *);
void print_header();
void new_conf();
void replica_init();
void langevin_init();
void moment_init();
void setup_lj();
void read_lambda();
void targetinput();
void printconf(Slice *);
void dosinput();

void set_up() {
  long itime;
  int i, j;
  double x, ener0;

  input("path.inp");
  sys.filep = stdout;

  slice = (Slice *)calloc(MAXPATH, sizeof(Slice));
  fwtrial = (Slice *)calloc(MAXPATH, sizeof(Slice));
  bwtrial = (Slice *)calloc(MAXPATH, sizeof(Slice));
  trial = (Slice *)calloc(MAXPATH, sizeof(Slice));

  path.ntarget = 1;
  sys.nstates = 1;
  sys.initial_state = 1;
  targetinput();
  moment_init();
  choose_new_momenta(&slice[0]);
  if (sys.nstates == 0) {
    confinput();
    add_newstate(&slice[0]);
  }

  langevin_init();

  sys.sqrttemp = 1. / sqrt(sys.beta);

  ener0 = total_energy(&slice[0]);
  gprint(ener0);
  if (fabs(ener0 - sys.fixed_energy) > 0.0001)
    fix_energy(&slice[0]);
  ener0 = total_energy(&slice[0]);
  gprint(ener0);
  printconf(&slice[0]);

  sys.nreplica = 20;
  sys.current_replica = 1;
  sys.scalefactor = 0.01;
  sys.maxlength = MAXPATH / 2;
  replica_init();

  printf("hier\n");
  if (sys.start_type == 1) {
    printf("reading previous path configuration\n");
    confinput();
    dosinput();
  }

  replica[sys.current_replica]->pathlen = path.nslices;
  dprint(sys.current_replica);
  dprint(replica[sys.current_replica]->pathlen);
  dprint(sys.initial_state);

  for (i = 0; i < replica[sys.current_replica]->pathlen; i++) {
    printf("%d ", in_state(&slice[i]));
    print_rc(&slice[i], sys.initial_state);
  }

  if (sys.current_replica == 0) {
    replica[sys.current_replica]->type = analyse_state_i(
        slice, replica[sys.current_replica],
        replica[sys.current_replica]->pathlen, sys.initial_state);
  } else {
    replica[sys.current_replica]->type =
        analyse(slice, replica[sys.current_replica],
                replica[sys.current_replica]->pathlen, sys.initial_state);
  }
  printf("replica %d  type %d\n", sys.current_replica,
         replica[sys.current_replica]->type);
  dprint(sys.initial_state);
  reset_center(&slice[0]);

  print_header();
  scr = nulscr;
  return;
}

void replica_init() {
  Replica *prep;
  int i, j, pathlen, maxlength, len, type;

  printf("in replica initialisation\n");

  read_lambda();

  for (j = 0; j < MAXSTATES; j++) {
    for (i = 0; i < sys.nreplica; i++) {
      state[j].srep[i].index = i;
      state[j].srep[i].swapindex = i;
      state[j].srep[i].lambda = sys.lambda[i];
    }
  }

  for (i = 0; i < sys.nreplica; i++) {
    replica[i] = &state[sys.initial_state - 1].srep[i];

    replica[i]->index = i;
    replica[i]->swapindex = i;
    replica[i]->lambda = sys.lambda[i];
  }

  prep = replica[0];
  maxlength = MAXPATH / 2;
  //  fwtrial[0]=slice[i];
  len = trajectory_state_i(&slice[0], prep, maxlength, sys.initial_state);
  //  len = maxlength;
  dprint(len);
  for (i = 0; i <= len; i++)
    slice[i + len + 1] = slice[i];
  for (i = 0; i <= len; i++)
    slice[len - i] = slice[i + len + 1];

  pathlen = 2 * len + 2;
  dprint(pathlen);
  for (i = 0; i < pathlen; i++)
    print_rc(&slice[i], sys.initial_state);

  type = analyse_state_i(slice, prep, pathlen, sys.initial_state);
  dprint(type);
  replica[0]->type = type;
  replica[0]->pathlen = pathlen;
  path.nslices = pathlen;
  sys.current_replica = 0;

  return;
}

void langevin_init() {
  double fdt, ex, temp, dt, s11, s22, s12, S;

  dt = sys.delta_t;
  fdt = sys.delta_t * langevin.gamma;
  ex = exp(-fdt / 2.0);
  temp = 1. / sys.beta;

  langevin.mu = 1;
  langevin.c1 = ex;
  langevin.c2 = (1 - langevin.c1) / langevin.gamma;
  langevin.c3 =
      sqrt(2 * (1 - langevin.c1 * langevin.c1) / (sys.beta * langevin.mu));

  gprint(langevin.mu);
  gprint(langevin.c1);

  gprint(langevin.c2);
  gprint(langevin.c3);
  gprint(langevin.gamma);
}

void moment_init() {
  double norm;
  int dim, i;

  dim = 3 * sys.npart;
  sys.norm = (double **)calloc(6, sizeof(double *));
  for (i = 0; i < 6; i++)
    sys.norm[i] = (double *)calloc(dim, sizeof(double));
  norm = 1. / sqrt(sys.npart);
  for (i = 0; i < sys.npart; i++) {
    sys.norm[0][i * 3] = norm;
    sys.norm[1][i * 3 + 1] = norm;
    sys.norm[2][i * 3 + 2] = norm;
  }
}

int input(char name[]) {

  char dum[40];
  FILE *fp;
  int i, idum;

  printf("in input\n");
  if ((fp = fopen(name, "r")) == NULL) {
    printf("input: can't open %s\n", name);
    return 1;
  } else {
    fscanf(fp, "%s", dum);
    fscanf(fp, "%d%s", &sys.sim_type, dum);
    fscanf(fp, "%d%s", &sys.start_type, dum);
    fscanf(fp, "%d%s", &sys.ensemble_type, dum);
    fscanf(fp, "%d%s", &sys.ncycle1, dum);
    fscanf(fp, "%d%s", &sys.ncycle2, dum);
    fscanf(fp, "%d%s", &sys.nshoot, dum);
    fscanf(fp, "%d%s", &sys.nswap, dum);
    fscanf(fp, "%d%s", &sys.nswapstates, dum);
    fscanf(fp, "%d%s", &sys.nreverse, dum);
    fscanf(fp, "%lf%s", &sys.beta, dum);
    fscanf(fp, "%lf%s", &sys.delta_t, dum);
    fscanf(fp, "%lf%s", &langevin.gamma, dum);
    fscanf(fp, "");
    fscanf(fp, "%d%s", &sys.freq_check, dum);
    fscanf(fp, "%d%s", &sys.freq_graphics, dum);
    fscanf(fp, "");
    fscanf(fp, "%s", dum);
    fscanf(fp, "%d%s", &path.nslices, dum);
    fscanf(fp, "%d%s", &path.ninter, dum);
    fscanf(fp, "%lf%s", &sys.fixed_energy, dum);
    fscanf(fp, "%lf%s", &sys.dvmax, dum);
    fscanf(fp, "%lf%s", &sys.min_stable, dum);

    fscanf(fp, "");
    fscanf(fp, "%s", dum);
    fscanf(fp, "%d%s", &sys.npart, dum);
    fscanf(fp, "%lf%s", &sys.boxl.x, dum);
    fscanf(fp, "%lf%s", &sys.boxl.y, dum);
    fscanf(fp, "%lf%s", &sys.boxl.z, dum);

    fclose(fp);
    printf("read input\n");

    sprintf(sys.block_stats.aver[0].name, "potential energy  ");
    sprintf(sys.block_stats.aver[1].name, "kinetic   energy   ");
    sprintf(sys.block_stats.aver[2].name, "total     energy   ");
    sprintf(sys.block_stats.aver[3].name, "momentum x         ");
    sprintf(sys.block_stats.aver[4].name, "momentum y         ");
    sprintf(sys.block_stats.aver[5].name, "momentum z         ");
    sprintf(sys.block_stats.aver[6].name, "angular momentum x ");
    sprintf(sys.block_stats.aver[7].name, "angular momentum y ");
    sprintf(sys.block_stats.aver[8].name, "angular momentum z ");
    sprintf(sys.block_stats.aver[9].name, "temp first slice   ");
    sprintf(sys.block_stats.aver[10].name, "temp last  slice   ");
    sprintf(sys.block_stats.aver[12].name, "hahb               ");
    sprintf(sys.block_stats.aver[13].name, "dhb               ");
    sprintf(sys.block_stats.aver[14].name, "number of blips   ");

    for (i = 0; i < STATINTERVAL; i++)
      sprintf(sys.block_stats.shoot[i].name, "slices in replica %3d ", i);
    for (i = 0; i < STATINTERVAL; i++)
      sprintf(sys.block_stats.mcacc[i].name, "shots replica     %3d ", i);
    for (i = 0; i < STATINTERVAL; i++)
      sprintf(sys.block_stats.mcacc[i + STATINTERVAL].name,
              "swaps replica     %3d ", i);
    for (i = 0; i < STATINTERVAL; i++)
      sprintf(sys.block_stats.mcacc[i + 2 * STATINTERVAL].name,
              "revs  replica     %3d ", i);
    for (i = 0; i < STATINTERVAL; i++)
      sprintf(sys.block_stats.mcacc[i + 3 * STATINTERVAL].name,
              "swap  states      %3d ", i);

    printf("set average names\n");

    nulstat = sys.final_stats = sys.block_stats;
  }
  return 1;
}

void print_header() {

  fprintf(sys.filep, "DYNAMICAL SAMPLING OF LANGEVIN ACTION\n\n");
  fprintf(sys.filep, "Peter Bolhuis 1997\n\n");
  fprintf(sys.filep, "SYSTEM PARAMETERS\n\n");
  fprintf(sys.filep, "Number of particles                %12d\n", sys.npart);
  fprintf(sys.filep, "Box length x                       %12g\n", sys.boxl.x);
  fprintf(sys.filep, "           y                       %12g\n", sys.boxl.y);
  fprintf(sys.filep, "           z                       %12g\n", sys.boxl.z);
  fprintf(sys.filep, "Fixed energy per time slice        %12g\n",
          sys.fixed_energy);
  fprintf(sys.filep, "Friction                           %12g\n", sys.gamma);
  fprintf(sys.filep, "Temperature (beta)                 %12g\n", sys.beta);
  fprintf(sys.filep, "Timestep                           %12g\n", sys.delta_t);
  fprintf(sys.filep, "Integration step bDdt              %12g\n", sys.bDdt);
  fprintf(sys.filep, "Sigma_r                            %12g\n", bivar.sig_r);
  fprintf(sys.filep, "Sigma_v                            %12g\n", bivar.sig_v);
  fprintf(sys.filep, "c_rv                               %12g\n", bivar.c_rv);

  fprintf(sys.filep, "\nPATH PARAMETERS\n\n");
  fprintf(sys.filep, "Number of timesteps                %12d\n", path.nslices);
  fprintf(sys.filep, "Number of iterations between slices%12d\n", path.ninter);
  fprintf(sys.filep, "Number of states                   %12d\n", sys.nstates);
  fprintf(sys.filep, "Initial total energy               %12g\n", path.energy);
  fprintf(sys.filep, "Quench                             %12d\n", path.quench);
  fprintf(sys.filep, "Constraint particle in A           %12d\n", path.tagged1);
  fprintf(sys.filep, "Number of required neighbors       %12d\n", path.nneigh1);
  fprintf(sys.filep, "Region A: minimum diff from target %12g\n",
          path.lower_moi_1);
  fprintf(sys.filep, "Region A: maximum diff from target %12g\n",
          path.upper_moi_1);
  fprintf(sys.filep, "Constraint particle in B           %12d\n", path.tagged2);
  fprintf(sys.filep, "Number of required neighbors       %12d\n", path.nneigh2);
  fprintf(sys.filep, "Region B: minimum diff from target %12g\n",
          path.lower_moi_2);
  fprintf(sys.filep, "Region B: maximum diff from target %12g\n",
          path.upper_moi_2);

  fprintf(sys.filep, "\nSIMULATION PARAMETERS\n\n");

  fprintf(sys.filep,
          ((sys.start_type == 0) ? "starting from random configuration"
                                 : "Read configuration from file fluid.inp"));
  fprintf(sys.filep, "  -  istart = %d\n", sys.start_type);

  fprintf(sys.filep, "Number of simulation blocks        %12d\n", sys.ncycle1);
  fprintf(sys.filep, "Number of cycles per block         %12d\n", sys.ncycle2);
  fprintf(sys.filep, "Total number of cycles             %12d\n",
          sys.ncycle1 * sys.ncycle2);
  fprintf(sys.filep, "Frequency of check                 %12d\n",
          sys.freq_check);
  fprintf(sys.filep, "Frequency of graphics              %12d\n",
          sys.freq_graphics);
}

void read_lambda() {

  char dum[240];
  FILE *fp;
  int i;

  if ((fp = fopen("lambda.inp", "r")) == NULL) {
    printf("inppt: can't open lambda.inp\n");
    return;
  } else {
    fscanf(fp, "%d %s", &sys.nreplica, dum);
    dprint(sys.nreplica);

    for (i = 0; i < sys.nreplica; i++) {
      // dprint(i);
      fscanf(fp, "%lf ", &sys.lambda[i]);
      printf("interface %d has lambda = %g\n", i, sys.lambda[i]);
    }
    fclose(fp);
  }
}

void targetinput() {
  Pts ptsdum, *psi, *psj;
  vector r_cm, dr, rij;
  double r2, r6, dr2, en = 0;
  char filename[200], dum[240];
  FILE *fp;
  int i, j, k, begin, end, nworm, not_done;

  printf("in target input\n");
  sprintf(filename, "target.inp");
  if ((fp = fopen(filename, "r")) == NULL) {
    fprintf(sys.filep, "input: can't open target.inp\n");
    return;
  } else {
    fscanf(fp, "%d %s", &path.ntarget, dum);
    fscanf(fp, "");
    sys.nstates = path.ntarget;
    for (i = 0; i < path.ntarget; i++) {
      for (j = 0; j < sys.npart; j++) {
        fscanf(fp, "%lf %lf %lf", &(state[i].target.pts[j].r.x),
               &(state[i].target.pts[j].r.y), &(state[i].target.pts[j].r.z));
      }
    }
    fclose(fp);
  }

  for (i = 0; i < path.ntarget; i++) {
    r_cm = nulvec;
    for (j = 0; j < sys.npart; j++) {
      vector_add(r_cm, state[i].target.pts[j].r, r_cm);
    }
    scalar_divide(r_cm, sys.npart, r_cm);
    for (j = 0; j < sys.npart; j++) {
      vector_minus(state[i].target.pts[j].r, r_cm, state[i].target.pts[j].r);
    }

    for (j = 0; j < sys.npart; j++) {
      state[i].target.pts[j].dr2 = 0;
      for (k = 0; k < sys.npart; k++)
        if (j != k) {
          vector_minus(state[i].target.pts[j].r, state[i].target.pts[k].r, dr);
          dr2 = vector_inp(dr, dr);
          state[i].target.pts[j].dr2 += dr2;
        }
      printf("target %d, particle %d, distance_sq = %g\n", i, j,
             state[i].target.pts[j].dr2);
    }

    printf("target %d, particle %d, total distance_sq = %g\n", i, j,
           state[i].target.pts[0].dr2);

    for (j = 0; j < sys.npart; j++)
      state[i].target.ordered_dr2[j] = state[i].target.pts[j].dr2;

    printf("Unsorted elements: \n");
    for (j = 0; j < sys.npart; j++) {
      printf("target %d, particle %d, distance_sq = %g\n", i, j,
             state[i].target.ordered_dr2[j]);
    }
    quicksort(state[i].target.ordered_dr2, 0, sys.npart - 1);
    printf("Sorted elements: \n");

    for (j = 0; j < sys.npart; j++) {
      printf("target %d, particle %d, distance_sq = %g\n", i, j,
             state[i].target.ordered_dr2[j]);
    }
  }

  for (i = 0; i < sys.nstates; i++) {
    state[i].min = sys.min_stable;
    state[i].min2 = sys.min_stable;
  }

  slice[0] = state[sys.initial_state - 1].target;

  slice[0].xi = 0;
  for (j = 1; j < path.nslices; j++)
    slice[j] = slice[0];
  printconf(&slice[0]);
  check_energy(&slice[0]);
  reset_center(&slice[0]);
  print_rc(&slice[0], sys.initial_state);
  printf("exiting target input\n");

  return;
}

void confinput() {
  char dum[40];
  char filename[240];
  double fdum;
  FILE *fp;
  int i, j, m, begin, end, nworm;

  sprintf(filename, "path.con");
  if ((fp = fopen(filename, "r")) == NULL) {
    printf("input: can't open %s\n", filename);
    return;
  } else {
    fscanf(fp, "%d %s", &path.nslices, dum);
    fscanf(fp, "%d %s", &sys.current_replica, dum);
    fscanf(fp, "%d %s", &sys.initial_state, dum);

    fscanf(fp, "%lf %lf %lf  %s", &sys.boxl.x, &sys.boxl.y, &sys.boxl.z, dum);
    for (j = 0; j < path.nslices; j++) {
      for (i = 0; i < sys.npart; i++) {
        fscanf(fp, "%lf %lf %lf %lf %lf %lf", &(slice[j].pts[i].r.x),
               &(slice[j].pts[i].r.y), &(slice[j].pts[i].r.z),
               &(slice[j].pts[i].v.x), &(slice[j].pts[i].v.y),
               &(slice[j].pts[i].v.z));
      }
    }
    fclose(fp);
  }
  dprint(path.nslices);
  replica[sys.current_replica]->pathlen = path.nslices;
  vprint(slice[path.nslices - 1].pts[sys.npart - 1].r);
  vprint(slice[path.nslices - 1].pts[sys.npart - 1].v);

  for (i = 0; i < path.nslices; i++)
    create_all_rc(&slice[i]);
}

void confoutput() {
  FILE *fp;
  char filename[240];
  int i, j;

  sprintf(filename, "path.con");
  if ((fp = fopen(filename, "w")) == NULL) {
    printf("output: can't open %s\n", filename);
    return;
  } else {
    fprintf(fp, "%d  slices \n", path.nslices);
    fprintf(fp, "%d  current_replica  \n", sys.current_replica);
    fprintf(fp, "%d  initial_state \n", sys.initial_state);

    fprintf(fp, "%lf %lf %lf   boxl,dr\n", sys.boxl.x, sys.boxl.y, sys.boxl.z);
    for (j = 0; j < path.nslices; j++) {
      for (i = 0; i < sys.npart; i++) {
        fprintf(fp, "%12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf\n",
                slice[j].pts[i].r.x, slice[j].pts[i].r.y, slice[j].pts[i].r.z,
                slice[j].pts[i].v.x, slice[j].pts[i].v.y, slice[j].pts[i].v.z);
      }
    }
    fclose(fp);
  }
}

void dosinput() {
  FILE *fp;
  char filename[240];
  int i, j, idum;

  sprintf(filename, "dos_all.inp");
  if ((fp = fopen(filename, "r")) == NULL) {
    printf("output: can't open %s\n", filename);
    return;
  } else {
    printf("reading dos data\n");
    fscanf(fp, "%lf", &sys.scalefactor);
    fscanf(fp, "");
    for (i = 0; i < sys.nstates; i++) {
      for (j = 0; j < sys.nreplica; j++) {
        fscanf(fp, "%d %lf\n", &idum, &state[i].srep[j].dos);
        printf("%d %lf\n", j, state[i].srep[j].dos);
      }
      fscanf(fp, "");
    }
    fclose(fp);
  }
}

void targetoutput() {
  FILE *fp;
  char filename[240];
  int i, j;

  sprintf(filename, "target.out");
  if ((fp = fopen(filename, "w")) == NULL) {
    printf("output: can't open %s\n", filename);
    return;
  } else {
    fprintf(fp, "%d nstates\n\n", sys.nstates);

    for (i = 0; i < sys.nstates; i++) {
      for (j = 0; j < sys.npart; j++) {
        fprintf(fp, "%lf %lf %lf\n", state[i].target.pts[j].r.x,
                state[i].target.pts[j].r.y, state[i].target.pts[j].r.z);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
}

void distoutput() {
  FILE *fp;
  char filename[100];
  double x, y[4];
  int i, j, k;

  for (k = 0; k < sys.nstates; k++) {
    sprintf(filename, "crosshistAB_%d.dat", k);
    if ((fp = fopen(filename, "w")) == NULL) {
      printf("output: can't open file.dat\n");
      return;
    } else {
      for (j = 0; j < sys.nreplica; j++) {
        for (i = 0; i < MBIN; i++) {

          x = (double)(i) / 50.;
          fprintf(fp, "%lf %g\n", x, scr.crosshistAB[k][j][i]);
        }
        fprintf(fp, "\n");
      }
      fclose(fp);
    }
  }

  for (k = 0; k < sys.nstates; k++) {
    sprintf(filename, "crosshistBA_%d.dat", k);
    if ((fp = fopen(filename, "w")) == NULL) {
      printf("output: can't open file.dat\n");
      return;
    } else {
      for (j = 0; j < sys.nreplica; j++) {
        for (i = 0; i < MBIN; i++) {
          x = (double)(i) / 50.;
          fprintf(fp, "%lf %g\n", x, scr.crosshistBA[k][j][i]);
        }
        fprintf(fp, "\n");
      }
      fclose(fp);
    }
  }

  for (k = 0; k < sys.nstates; k++) {
    sprintf(filename, "fehistAB_%d.dat", k);
    if ((fp = fopen(filename, "w")) == NULL) {
      printf("output: can't open file.dat\n");
      return;
    } else {
      for (j = 0; j < sys.nreplica; j++) {
        for (i = 0; i < MBIN; i++) {

          x = 0.2 * (i - MBIN / 2);
          fprintf(fp, "%lf %d\n", x, scr.fehistAB[k][j][i]);
        }
        fprintf(fp, "\n");
      }
      fclose(fp);
    }
  }

  for (k = 0; k < sys.nstates; k++) {
    sprintf(filename, "fehistBA_%d.dat", k);
    if ((fp = fopen(filename, "w")) == NULL) {
      printf("output: can't open file.dat\n");
      return;
    } else {
      for (j = 0; j < sys.nreplica; j++) {
        for (i = 0; i < MBIN; i++) {

          x = 0.2 * (i - MBIN / 2);
          fprintf(fp, "%lf %d\n", x, scr.fehistBA[k][j][i]);
        }
        fprintf(fp, "\n");
      }
      fclose(fp);
    }
  }
  if ((fp = fopen("boundhistAB.dat", "w")) == NULL) {
    printf("output: can't open file.dat\n");
    return;
  } else {
    for (j = 0; j < sys.nreplica; j++) {
      for (i = 0; i < MBIN; i++) {

        x = 0.2 * (i - MBIN / 2);
        fprintf(fp, "%lf %d\n", x, scr.boundhistAB[j][i]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
  if ((fp = fopen("boundhistBA.dat", "w")) == NULL) {
    printf("output: can't open file.dat\n");
    return;
  } else {
    for (j = 0; j < sys.nreplica; j++) {
      for (i = 0; i < MBIN; i++) {

        x = 0.2 * (i - MBIN / 2);
        fprintf(fp, "%lf %d\n", x, scr.boundhistBA[j][i]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  if ((fp = fopen("densmat.dat", "w")) == NULL) {
    printf("output: can't open file.dat\n");
    return;
  } else {
    for (i = 0; i < LATBIN; i++) {
      if (i == 0)
        fprintf(fp, "{");
      for (j = 0; j < LATBIN; j++) {
        if (j == 0)
          fprintf(fp, "{");
        fprintf(fp, "%d", scr.densmat[j][i]);
        if (j < LATBIN - 1)
          fprintf(fp, ", ");
      }
      if (i < LATBIN - 1)
        fprintf(fp, "},\n");
    }
    fprintf(fp, "}}\n");

    fclose(fp);
  }

  if ((fp = fopen("currentmat.dat", "w")) == NULL) {
    printf("output: can't open file.dat\n");
    return;
  } else {

    for (i = 0; i < LATBIN; i++) {
      for (j = 0; j < LATBIN; j++) {
        if (scr.ncurrentmat[j][i] > 10) {
          scr.currentmat[j][i][0] /= scr.ncurrentmat[j][i];
          scr.currentmat[j][i][1] /= scr.ncurrentmat[j][i];
        } else {
          scr.currentmat[j][i][0] = 0;
          scr.currentmat[j][i][1] = 0;
        }
      }
    }
    for (i = 0; i < LATBIN; i++) {
      if (i == 0)
        fprintf(fp, "{");
      for (j = 0; j < LATBIN; j++) {
        if (j == 0)
          fprintf(fp, "{");
        fprintf(fp, "{%g,%g}", scr.currentmat[j][i][0],
                scr.currentmat[j][i][1]);
        if (j < LATBIN - 1)
          fprintf(fp, ", ");
      }
      if (i < LATBIN - 1)
        fprintf(fp, "},\n");
    }
    fprintf(fp, "}}\n");

    fclose(fp);
  }

  if ((fp = fopen("dos.dat", "w")) == NULL) {
    printf("output: can't open file.dat\n");
    return;
  } else {
    for (j = 0; j < sys.nreplica; j++) {
      fprintf(fp, "%d %g\n", j, replica[j]->dos);
    }
    fclose(fp);
  }

  if ((fp = fopen("dos_all.dat", "w")) == NULL) {
    printf("output: can't open file.dat\n");
    return;
  } else {
    fprintf(fp, "%lf\n\n", sys.scalefactor);

    for (i = 0; i < sys.nstates; i++) {
      for (j = 0; j < sys.nreplica; j++) {
        fprintf(fp, "%d %g\n", j, state[i].srep[j].dos);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  if ((fp = fopen("dos_lambda.dat", "w")) == NULL) {
    printf("output: can't open file.dat\n");
    return;
  } else {
    for (i = 0; i < sys.nstates; i++) {
      for (j = 1; j < sys.nreplica; j++) {
        fprintf(fp, "%g %g\n", sys.lambda[j],
                state[i].srep[j].dos - state[i].srep[1].dos);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  if ((fp = fopen("replicacount.dat", "w")) == NULL) {
    printf("output: can't open file.dat\n");
    return;
  } else {
    for (j = 0; j < sys.nreplica; j++) {
      fprintf(fp, "%d %d\n", j, replica[j]->ntotal);
    }
    fclose(fp);
  }

  if ((fp = fopen("replicacount_all.dat", "w")) == NULL) {
    printf("output: can't open file.dat\n");
    return;
  } else {
    for (i = 0; i < sys.nstates; i++) {
      for (j = 0; j < sys.nreplica; j++) {
        fprintf(fp, "%d %d\n", j, state[i].srep[j].ntotal);
      }
      fprintf(fp, "\n");
    }

    fclose(fp);
  }
}

void printconf(Slice *psl) {
  int j;

  for (j = 0; j < sys.npart; j++) {
    vprint(psl->pts[j].r);
    vprint(psl->pts[j].v);
    vprint(psl->pts[j].f);
  }
  return;
}
