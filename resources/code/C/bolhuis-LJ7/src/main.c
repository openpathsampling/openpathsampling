#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "path.h"

void do_pathsampling();
void do_cycle();
void accumulate_stats();
void terminate_block(Stats *, Stats *);
void finalstat(Stats *);
void printstat(Stats *);
void printstatusline();

void crossinghistogram(Replica *);
void fe_histogram(Replica *);
void get_flux(Replica *);
void followroute(int, int);
void write_hists(int);
void update_densityplot();
void openpathfiles();
void dump_paths(FILE *);
void closepathfiles();

int initial_path();
void reset_center(Slice *);
void rescale_velocities();
void fixed_energy();
void spectrum(Slice *psl);

int main(int argc, char *argv[]) {
  FILE *fp, *fp2, *fp3, *fp4;
  int i, j, k, l, min, max, count = 0;
  ;
  double x, y;
  vector rtmp;

  set_up();
  if (sys.freq_graphics) {
    Init_Graphics(argc, argv, &slice[0]);
  }
  //  check_energy();

  if ((fp = fopen("swap.dat", "w")) == NULL) {
    printf("output: can't open swap.out\n");
    return 0;
  }

  if ((fp2 = fopen("route.dat", "w")) == NULL) {
    printf("output: can't open route.out\n");
    return 0;
  }

  if ((fp3 = fopen("rc.dat", "w")) == NULL) {
    printf("output: can't open rc.out\n");
    return 0;
  }
  fprintf(fp3, "{");
  if ((fp4 = fopen("dos_track.dat", "w")) == NULL) {
    printf("output: can't open dos_track.out\n");
    return 0;
  }

  for (j = 1; j <= sys.ncycle1; j++) {
    for (i = 1; i <= sys.ncycle2; i++) {
      if (sys.sim_type == 1) {
        do_pathsampling();
        fprintf(fp, "%d ", count++);
        for (k = 0; k < sys.nreplica; k++)
          fprintf(fp, "%d ", replica[k]->swapindex);
        fprintf(fp, "\n");
        for (k = 0; k < sys.nreplica; k++)
          fprintf(fp2, "%2d ", scr.route[k]);
        fprintf(fp2, "\n");

        if (i % sys.freq_check == 0)
          dump_paths(fp3);
      } else {
        do_cycle();
        fprintf(fp3, "%g \n", get_rc(&slice[0], sys.initial_state));
        check_energy(&slice[0]);
      }
      accumulate_stats();
      if (i % sys.freq_check == 0) {
        check_energy(&slice[0]);
        //	reset_center(&slice[0]);
        // check_energy(&slice[0]);
      }
    }

    printf("\n\nBLOCK NR  %d\n\n", j);
    terminate_block(&sys.block_stats, &sys.final_stats);
    //  write_hists(j);
    printf("hier \n");
    dprint(sys.nreplica);
    gprint(replica[1]->dos);
    for (k = 0; k < sys.nreplica; k++)
      printf("%d %g\n", k, replica[k]->dos);

    for (l = 0; l < sys.nstates; l++) {
      for (k = 0; k < sys.nreplica; k++)
        fprintf(fp4, "%g\n", state[l].srep[k].dos);
    }
    fprintf(fp4, "\n");

    max = 0;
    min = 10000000;
    for (l = 0; l < sys.nstates; l++) {
      for (k = 0; k < sys.nreplica; k++) {
        if (state[l].srep[k].ntotal > max)
          max = state[l].srep[k].ntotal;
        if (state[l].srep[k].ntotal < min)
          min = state[l].srep[k].ntotal;
      }
    }
    if ((double)max / min < 1.1) {
      printf("histogram flat within 10 percent now reducing scalefactor\n");
      sys.scalefactor /= 2.0;
      gprint(sys.scalefactor);
      printf("resetting histogram\n");
      for (l = 0; l < sys.nstates; l++) {
        for (k = 0; k < sys.nreplica; k++)
          state[l].srep[k].ntotal = 0;
      }
    } else {
      printf("histogram not flat yet:  min = %d, max = %d \n", min, max);
    }
  }
  printf("\n\nFINAL STATISTICS  %d\n\n", j);
  finalstat(&sys.final_stats);
  confoutput();
  distoutput();
  targetoutput();

  fclose(fp);
  fclose(fp2);
  fprintf(fp3, "}");
  fclose(fp3);
  fclose(fp4);
  printf("turnovers:");
  for (k = 0; k < sys.nreplica; k++)
    printf("%3d ", scr.turnovers[k]);
  printf("\n");
  return 1;
}

void mainloop_for_graphics() {
  static int icycle1, icycle2, initial;
  int i;

  if (initial == 0) {
    initial = 1;
    icycle1 = icycle2 = 1;
  }

  if (sys.sim_type == 1) {
    do_pathsampling();
  } else {
    do_cycle();
  }
  accumulate_stats();

  if (icycle2 % sys.freq_check == 0) {
	  check_energy(&slice[0]);
    dprint(path.nslices);
  }

  if (icycle2 % sys.ncycle2 == 0) {

    fprintf(sys.filep, "\n\nBLOCK NR  %d\n\n", icycle1);
    terminate_block(&sys.block_stats, &sys.final_stats);
    icycle2 = 0;
    if (icycle1 % sys.ncycle1 == 0) {
      fprintf(sys.filep, "\n\nFINAL STATISTICS \n\n");
      finalstat(&sys.final_stats);
      confoutput();
      distoutput();
      targetoutput();
      exit(1);
    }
    icycle1++;
  }

  icycle2++;
  return;
}

void do_pathsampling() {

  int i, j, iwhich, irep, jrep, drep;
  double en;

  for (i = 0; i < sys.nshoot + sys.nswap + sys.nreverse + sys.nswapstates; i++) {
    iwhich = (int)(ran3() * (sys.nshoot + sys.nswap + sys.nreverse + sys.nswapstates));

    if (iwhich < sys.nshoot) {
      j = sys.current_replica;

      if (shoot(replica[j]))
        sys.block_stats.mcacc[j].acc++;
      sys.block_stats.mcacc[j].try ++;

    } else {
      if (iwhich < sys.nshoot + sys.nswap) {
        irep = sys.current_replica;

        drep = 2 * ((int)(2 * ran3())) - 1;
        if (abs(drep) != 1)
          printf("error:drep = %d\n", drep);
        jrep = irep + drep;
        if (swap_replica(irep, jrep))
          sys.block_stats.mcacc[irep + STATINTERVAL].acc++;
        sys.block_stats.mcacc[irep + STATINTERVAL].try ++;

      } else {
        if (iwhich < sys.nshoot + sys.nswap + sys.nswapstates) {
			irep = sys.current_replica;
          jrep = (int)(ran3() * (sys.nreplica - 1)) + 1;
          if (irep == sys.current_replica) {
            if (swap_states(irep, jrep))
              sys.block_stats.mcacc[irep + 3 * STATINTERVAL].acc++;
            sys.block_stats.mcacc[irep + 3 * STATINTERVAL].try ++;
          }

        } else {
          if (iwhich <
              sys.nshoot + sys.nswap + +sys.nswapstates + sys.nreverse) {
            irep = sys.current_replica;
            if (reverse_replica(irep))
              sys.block_stats.mcacc[irep + 2 * STATINTERVAL].acc++;
            sys.block_stats.mcacc[irep + 2 * STATINTERVAL].try ++;
          }
        }
      }
    }

    crossinghistogram(replica[sys.current_replica]);
    fe_histogram(replica[sys.current_replica]);
    get_flux(replica[sys.current_replica]);
  }
  printstatusline();

}

void do_cycle() {

  int i, j, iwhich, fix;
  double en;

  if (sys.sim_type == 3) {
    for (i = 0; i < 1; i++) {
      propagate_NHL(&slice[0]);
    }
  }

  if (sys.sim_type == 4) {
    path.nslices = 100;
    for (i = 0; i < path.nslices; i++) {
      propagate_NHL(&slice[i]);
      if (i < path.nslices - 1)
        slice[i + 1] = slice[i];
    }
    slice[0] = slice[path.nslices - 1];
  }
}

void followroute(int irep, int jrep) {
  Replica *prep;
  int i, j, k, swapindex, ilam;
  return;
}

void crossinghistogram(Replica *prep) {
  int irep, type, i, imax, imin, ilam, k;
  double lambdamax, lambdamin, z, weight;

  ilam = -1;
  prep->ntotal++;
  prep->navlen++;
  prep->avlen += prep->pathlen;
  irep = prep->index;

  type = prep->type;

  scr.pathtype[irep][type - 1]++;
  scr.type_mat[sys.initial_state - 1][irep][type - 1]++;

  update_average(sys.block_stats.shoot[irep], prep->pathlen);

  if (type == 0) {
    printf("error: incorrect path in crossinghist replica %d\n", irep);
  }
  if ((type == 1) || (type == 2)) {

    lambdamax = -100000;
    for (i = 0; i < prep->pathlen; i++) {
      z = get_rc(&slice[i], sys.initial_state);

      if (z > lambdamax)
        lambdamax = z;
    }

    if (type == 2)
      lambdamax = 15;

    k = irep;
    imax = (int)(lambdamax * 50);
    imin = (int)(replica[k]->lambda * 50);
    weight = 1;

    for (i = imin; i <= imax; i++) {
      if ((i > 0) && (i < MBIN))
        scr.crosshistAB[sys.initial_state - 1][k][i] += weight;
    }

    for (i = imin; i <= imax; i++) {
      if ((i > 0) && (i < MBIN))
        scr.crosshistAB_block[irep][i]++;
    }
  }

  if ((type == 3) || (type == 4)) {

    lambdamin = 100000;
    for (i = 0; i < prep->pathlen; i++) {
      z = get_rc(&slice[i], sys.initial_state);
      if (z < lambdamin)
        lambdamin = z;
    }

    imin = (int)(lambdamin * 1000);
    imax = (int)(prep->lambda * 1000);
    for (i = imin; i <= imax; i++) {
      if ((i > 0) && (i < MBIN))
        scr.crosshistBA[sys.initial_state - 1][irep][i]++;
    }
    for (i = imin; i <= imax; i++) {
      if ((i > 0) && (i < MBIN))
        scr.crosshistBA_block[irep][i]++;
    }
  }
}

void fe_histogram(Replica *prep) {
  int irep, type, i, ix, reachednext, imax;
  double z;

  irep = prep->index;

  if (irep == 0) {
    for (i = 0; i < prep->pathlen; i++) {
      z = get_rc(&slice[i], sys.initial_state);
      if (z < prep->lambda) {
        ix = (int)(z * 100);
        if ((ix > 0) && (ix < MBIN))
          scr.fehistAB[sys.initial_state - 1][irep][ix]++;
      }
    }
    return;
  }

  if (irep == sys.nreplica - 1) {
    for (i = 0; i < prep->pathlen; i++) {
      z = get_rc(&slice[i], sys.initial_state);
      if (z > prep->lambda) {
        ix = (int)(z * 100);
        if ((ix > 0) && (ix < MBIN))
          scr.fehistBA[sys.initial_state - 1][irep][ix]++;
      }
    }

    analyse_full(slice, prep, prep->pathlen, sys.initial_state);
    return;
  }

  type = prep->type;
  if (type == 0) {
    printf("error: incorrect path in crossinghist replica %d\n", irep);
  }
  if (type == 1) {
    for (i = 0; i < prep->pathlen; i++) {
      z = get_rc(&slice[i], sys.initial_state);
      if (z > prep->lambda) {
        ix = (int)(z * 5);
        if ((ix > 0) && (ix < MBIN))
          scr.fehistAB[sys.initial_state - 1][irep][ix]++;
      }
    }
  }
  if (type == 2) {
    for (i = 0; i < prep->pathlen; i++) {
      z = get_rc(&slice[i], sys.initial_state);
      if (z > prep->lambda) {
        ix = (int)(z * 5);
        if ((ix > 0) && (ix < MBIN))
          scr.fehistAB[sys.initial_state - 1][irep][ix]++;

        if (z < replica[irep + 1]->lambda) {
          if ((ix > 0) && (ix < MBIN))
            scr.boundhistAB[irep][ix]++;
        }
      }
    }
  }

  if (type == 4) {
    for (i = 0; i < prep->pathlen; i++) {
      z = get_rc(&slice[i], sys.initial_state);
      if (z < prep->lambda) {
        ix = (int)(z * 5);
        if ((ix > 0) && (ix < MBIN))
          scr.fehistBA[sys.initial_state - 1][irep][ix]++;
      }
    }
  }

  if (type == 3) {
    for (i = 0; i < prep->pathlen; i++) {
      z = get_rc(&slice[i], sys.initial_state);
      if (z < prep->lambda) {
        ix = (int)(z * 5);
        if ((ix > 0) && (ix < MBIN))
          scr.fehistBA[sys.initial_state - 1][irep][ix]++;
        if ((z > replica[irep - 1]->lambda)) {
          if ((ix > 0) && (ix < MBIN))
            scr.boundhistBA[irep][ix]++;
        }
      }
    }
  }
}

void get_flux(Replica *prep) {
  int i, len;

  if (prep->index == 0) {
    for (i = prep->pathlen - 1; i >= 0; i--) {
      if (in_state(&slice[i])) {
        len = i;
        break;
      }
    }

    scr.flux0[sys.initial_state - 1] += len;
    scr.nflux0[sys.initial_state - 1]++;
  }

  if (prep->index == 1) {
    for (i = prep->pathlen - 1; i >= 0; i--) {
      if (in_upper_window(&slice[i], replica[1]->lambda, sys.initial_state)) {
        len = i;
        break;
      }
    }
    scr.flux1[sys.initial_state - 1] += len;
    scr.nflux1[sys.initial_state - 1]++;
  }
}

void accumulate_stats() {
  vector r_cm, dr, v, r_tag, r_far, r_av, r, momentum, angular, out;
  int i, j, nn, na, nb, still_in_A, not_yet_in_B, t_A, t_B, t_r;
  double en, kin_ener, dr2, l, temp0, tempn;

  path.tot_ener = path.kin_ener = path.energy = 0;
  if (sys.sim_type == 3) {
    en = energy(&slice[0]);
    path.energy += en;
    kin_ener = 0;
    for (j = 0; j < sys.npart; j++) {
      v = slice[0].pts[j].v;
      kin_ener += 0.5 * vector_inp(v, v);
    }
    path.kin_ener += kin_ener;
    en += kin_ener;
    path.tot_ener += en;
  } else {

    for (i = 0; i <= 0 * path.nslices; i++) {
      en = energy(&slice[i]);
      scr.energy[i] += en;
      scr.energy2[i] += en * en;
      path.energy += en;
      kin_ener = 0;
      for (j = 0; j < sys.npart; j++) {
        v = slice[i].pts[j].v;
        kin_ener += 0.5 * vector_inp(v, v);
      }
      scr.kin_ener[i] += kin_ener;
      scr.kin_ener2[i] += kin_ener * kin_ener;
      path.kin_ener += kin_ener;
      en += kin_ener;
      scr.tot_ener[i] += en;
      scr.tot_ener2[i] += en * en;
      path.tot_ener += en;
    }
  }
  scr.fcount++;

  r_cm = nulvec;
  for (j = 0; j < sys.npart; j++) {
    vector_add(r_cm, slice[0].pts[j].r, r_cm);
  }
  scalar_divide(r_cm, sys.npart, r_cm);

  temp0 = 0;
  angular = momentum = nulvec;
  for (j = 0; j < sys.npart; j++) {
    v = slice[0].pts[j].v;
    r = slice[0].pts[j].r;
    vector_minus(r, r_cm, r);
    vector_add(momentum, v, momentum);
    vector_cross(r, v, out);
    vector_add(angular, out, angular);
    temp0 += vector_inp(v, v);
  }
  temp0 /= 3 * sys.npart - 3;
  scr.n++;

  POTENERGY.now = path.energy;
  KIN_ENER.now = path.kin_ener;
  TOT_ENER.now = path.tot_ener;
  TEMP.now = temp0;
  TEMP_L.now = tempn;

  PX_MOM.now = momentum.x;
  PY_MOM.now = momentum.y;
  PZ_MOM.now = momentum.z;
  LX_MOM.now = angular.x;
  LY_MOM.now = angular.y;
  LZ_MOM.now = angular.z;

  for (i = 0; i < NSTAT; i++) {
    sys.block_stats.aver[i].sum += sys.block_stats.aver[i].now;
    sys.block_stats.aver[i].sumsq +=
        sys.block_stats.aver[i].now * sys.block_stats.aver[i].now;
    sys.block_stats.aver[i].n++;
  }
}

void terminate_block(Stats *psb, Stats *psav) {
  int i;

  for (i = 0; i < NSTAT; i++) {
    psb->aver[i].sum /= sys.ncycle2;
    psb->aver[i].sumsq = sqrt(psb->aver[i].sumsq / sys.ncycle2 -
                              psb->aver[i].sum * psb->aver[i].sum);
  }
  for (i = 0; i < NSHOOTAV; i++) {
    psb->shoot[i].sum /= psb->shoot[i].n;
    psb->shoot[i].sumsq = sqrt(psb->shoot[i].sumsq / psb->shoot[i].n -
                               psb->shoot[i].sum * psb->shoot[i].sum);
  }

  for (i = 0; i < NACC; i++)
    psb->mcacc[i].ratio = (double)psb->mcacc[i].acc / psb->mcacc[i].try;

  printstat(psb);

  for (i = 0; i < NSTAT; i++) {
    psav->aver[i].sum += psb->aver[i].sum;
    psav->aver[i].sumsq += psb->aver[i].sum * psb->aver[i].sum;
    psb->aver[i].sum = 0;
    psb->aver[i].sumsq = 0;
    psb->aver[i].now = 0;
  }
  for (i = 0; i < NSHOOTAV; i++) {
    psav->shoot[i].sum += psb->shoot[i].sum;
    psav->shoot[i].sumsq += psb->shoot[i].sum * psb->shoot[i].sum;
    psb->shoot[i].sum = 0;
    psb->shoot[i].sumsq = 0;
    psb->shoot[i].now = 0;
    psb->shoot[i].n = 0;
  }

  for (i = 0; i < NACC; i++) {
    psav->mcacc[i].acc += psb->mcacc[i].acc;
    psav->mcacc[i].try += psb->mcacc[i].try;
    psb->mcacc[i].acc = 0;
    psb->mcacc[i].try = 0;
  }
}

void finalstat(Stats *psav) {
  int i;

  for (i = 0; i < NSTAT; i++) {
    psav->aver[i].sum /= sys.ncycle1;
    psav->aver[i].sumsq = sqrt(psav->aver[i].sumsq / sys.ncycle1 -
                               psav->aver[i].sum * psav->aver[i].sum);
  }
  for (i = 0; i < NSHOOTAV; i++) {
    psav->shoot[i].sum /= sys.ncycle1;
    psav->shoot[i].sumsq = sqrt(psav->shoot[i].sumsq / sys.ncycle1 -
                                psav->shoot[i].sum * psav->shoot[i].sum);
  }

  for (i = 0; i < NACC; i++)
    psav->mcacc[i].ratio = (double)psav->mcacc[i].acc / psav->mcacc[i].try;
  printstat(psav);
}

void printstat(Stats *psav) {
  int i, j, k, ntot[MAXSTATES];
  double f0, f1, flux[MAXSTATES];

  printf("\nAcceptances\n");
  for (i = 0; i < NACC; i++) {
    printf("%7d accepted %sout of %7d trial moves, ratio = %g\n",
           psav->mcacc[i].acc, psav->mcacc[i].name, psav->mcacc[i].try,
           psav->mcacc[i].ratio);
  }
  printf("Slice number averages\n");
  for (i = 0; i < NSHOOTAV; i++) {
    printf("     %s   = %10g +- %7g\n", psav->shoot[i].name, psav->shoot[i].sum,
           psav->shoot[i].sumsq);
  }

  printf("\nPathlength averages        ");
  for (i = 1; i <= sys.nstates; i++) {
    printf("%6d ", i);
  }
  printf("\n");

  for (k = 0; k < MAXREPLICA; k++) {
    printf("     Path length replica %2d ", k);
    for (i = 0; i < sys.nstates; i++) {
      printf("%6.1lf ",
             (double)state[i].srep[k].avlen / state[i].srep[k].navlen);
    }
    printf("\n");
  }

  printf("\nPathtype averages        ");
  for (i = 1; i <= sys.nstates; i++) {
    printf("  %d->%d  %d->j", i, i, i);
  }
  printf("\n");

  for (k = 0; k < MAXREPLICA; k++) {
    printf("     Path type replica %2d ", k);
    for (i = 0; i < sys.nstates; i++) {
      for (j = 0; j < 2; j++) {
        printf("%5d ", scr.type_mat[i][k][j]);
      }
    }
    printf("\n");
  }

  printf("\nMSTIS matrix  ");
  for (i = 0; i < sys.nstates; i++) {
    printf("    %2d", i + 1);
  }
  printf("\n");
  for (i = 0; i < sys.nstates; i++)
    ntot[i] = 0;
  for (k = 0; k < sys.nstates; k++) {
    printf("     state  %2d ", k + 1);
    for (i = 0; i < sys.nstates; i++) {
      printf("%5d ", scr.mstis_mat[i][k]);
      ntot[i] += scr.mstis_mat[i][k];
    }
    printf("\n");
  }
  printf("     total     ");
  for (i = 0; i < sys.nstates; i++)
    printf("%5d ", ntot[i]);
  printf("\n");

  printf("\nNormalized matrix   ");
  for (i = 0; i < sys.nstates; i++) {
    printf("%2d      ", i);
  }
  printf("\n");
  for (k = 0; k < sys.nstates; k++) {
    printf("     state  %2d ", k);
    for (i = 0; i < sys.nstates; i++) {
      printf("%6.5lf ", (double)scr.mstis_mat[i][k] / ntot[i]);
    }
    printf("\n");
  }

  printf("\nFluxes     0        1     total     flux\n");
  for (i = 0; i < sys.nstates; i++) {
    f0 = (double)scr.flux0[i] / scr.nflux0[i];
    f1 = (double)scr.flux1[i] / scr.nflux1[i];
    flux[i] = 1. / (f0 + f1);
    printf("%2d      ", i);
    printf("%7.2lf %7.2lf %7.2lf %7.6lf \n", f0, f1, (f0 + f1), 1. / (f0 + f1));
  }
  printf("\n");

  printf("\nCrossing prob from dos at replica %d\n", sys.nreplica - 1);
  for (i = 0; i < sys.nstates; i++) {
    printf("%2d      ", i);
    printf("%7.2lf %7.2lf %12.5g \n", state[i].srep[sys.nreplica - 1].dos,
           state[i].srep[sys.nreplica - 1].dos - state[i].srep[1].dos,
           exp(state[i].srep[sys.nreplica - 1].dos - state[i].srep[1].dos));
  }
  printf("\n");

  printf("\nRate matrix   ");
  for (i = 0; i < sys.nstates; i++) {
    printf("%2d      ", i);
  }
  printf("\n");
  for (k = 0; k < sys.nstates; k++) {
    printf("     state  %2d ", k);
    for (i = 0; i < sys.nstates; i++) {
      scr.rate[i][k] =
          exp(state[i].srep[sys.nreplica - 1].dos - state[i].srep[1].dos) *
          flux[i] * scr.mstis_mat[i][k] / ntot[i];
      printf("%12.5g ", (double)scr.rate[i][k]);
    }
    printf("\n");
  }

  printf("\nAverages\n");
  for (i = 0; i < NSTAT; i++) {
    printf("     %s   = %10g +- %7g\n", psav->aver[i].name, psav->aver[i].sum,
           psav->aver[i].sumsq);
  }

  fflush(NULL);
}

void printstatusline() {
  int hb_first, ha_last, nbar, i, j;

  i = in_state(&slice[0]);
  j = in_state(&slice[path.nslices - 1]);

  printf("status: state= %d replica=%d path-type=%d path-length %d initial %d, "
         "final %d\n",
         sys.initial_state, sys.current_replica,
         replica[sys.current_replica]->type, path.nslices, i, j);
  return;
}

void write_hists(int blocknum) {
  FILE *fp;
  double x, y[4];
  int i, j;
  char filename[100], num[10];

  printf("begin\n");
  sprintf(num, "%d", 100 + blocknum);
  sprintf(filename, "chistAB_%s.dat", num + 1);
  printf("now printing histogram %d to %s\n", blocknum, filename);

  if ((fp = fopen(filename, "w")) == NULL) {
    printf("output: can't open file.dat\n");
    return;
  } else {
    for (j = 0; j < sys.nreplica; j++) {
      for (i = 0; i < MBIN; i++) {

        x = 0.01 * (i - MBIN / 2);
        fprintf(fp, "%lf %d\n", x, scr.crosshistAB_block[j][i]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  sprintf(filename, "chistBA_%s.dat", num + 1);
  if ((fp = fopen(filename, "w")) == NULL) {
    printf("output: can't open file.dat\n");
    return;
  } else {
    for (j = 0; j < sys.nreplica; j++) {
      for (i = 0; i < MBIN; i++) {

        x = 0.01 * (i - MBIN / 2);
        fprintf(fp, "%lf %d\n", x, scr.crosshistBA_block[j][i]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  sprintf(filename, "chistAB_L_%s.dat", num + 1);
  if ((fp = fopen(filename, "w")) == NULL) {

    printf("output: can't open file.dat\n");
    return;
  } else {
    for (j = 0; j < sys.nreplica; j++) {
      for (i = 0; i < MBIN; i++) {

        x = 0.01 * (i - MBIN / 2);
        fprintf(fp, "%lf %d\n", x, scr.crosshistAB_L_block[j][i]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
  sprintf(filename, "chistBA_L_%s.dat", num + 1);
  if ((fp = fopen(filename, "w")) == NULL) {

    printf("output: can't open file.dat\n");
    return;
  } else {
    for (j = 0; j < sys.nreplica; j++) {
      for (i = 0; i < MBIN; i++) {

        x = 0.01 * (i - MBIN / 2);
        fprintf(fp, "%lf %d\n", x, scr.crosshistBA_L_block[j][i]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
  sprintf(filename, "chistAB_H_%s.dat", num + 1);
  if ((fp = fopen(filename, "w")) == NULL) {

    printf("output: can't open file.dat\n");
    return;
  } else {
    for (j = 0; j < sys.nreplica; j++) {
      for (i = 0; i < MBIN; i++) {

        x = 0.01 * (i - MBIN / 2);
        fprintf(fp, "%lf %d\n", x, scr.crosshistAB_H_block[j][i]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
  sprintf(filename, "chistBA_H_%s.dat", num + 1);
  if ((fp = fopen(filename, "w")) == NULL) {

    printf("output: can't open file.dat\n");
    return;
  } else {
    for (j = 0; j < sys.nreplica; j++) {
      for (i = 0; i < MBIN; i++) {

        x = 0.01 * (i - MBIN / 2);
        fprintf(fp, "%lf %d\n", x, scr.crosshistBA_H_block[j][i]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  for (j = 0; j < sys.nreplica; j++) {
    for (i = 0; i < MBIN; i++) {
      scr.crosshistAB_block[j][i] = 0;
      scr.crosshistAB_L_block[j][i] = 0;
      scr.crosshistAB_H_block[j][i] = 0;
      scr.crosshistBA_block[j][i] = 0;
      scr.crosshistBA_L_block[j][i] = 0;
      scr.crosshistBA_H_block[j][i] = 0;
    }
  }

  printf("end\n");
}

void reset_center(Slice *psl) {
  int i, j;
  vector r_cm, momentum, v, r;
  double angular;

  r_cm = nulvec;
  for (j = 0; j < sys.npart; j++) {
    vector_add(r_cm, psl->pts[j].r, r_cm);
  }
  scalar_divide(r_cm, sys.npart, r_cm);
  for (i = 0; i < path.nslices; i++) {
    for (j = 0; j < sys.npart; j++) {
      vector_minus(slice[i].pts[j].r, r_cm, slice[i].pts[j].r);
    }
  }
  return;
}

void fixed_energy() {
  double fac, kin_ener, kin_fixed, fixed_energy, pot_ener, pot_tot, kin_tot,
      maxener;

  int i, j;

  fixed_energy = sys.fixed_energy;
  maxener = -10000;
  for (i = 0; i < path.nslices; i++) {
    if ((pot_ener = energy(&slice[i])) > maxener) {
      maxener = pot_ener;
    }
  }
  if (maxener > sys.fixed_energy) {
    fprintf(sys.filep, "Warning: maximum energy larger than fixed energy %g\n",
            sys.fixed_energy - maxener);
    fprintf(sys.filep, "rescaling to %g\n", maxener);
    fixed_energy = maxener;
  }

  for (i = 0; i < path.nslices; i++) {
    kin_ener = 0;
    for (j = 0; j < sys.npart; j++) {
      kin_ener += 0.5 * vector_inp(slice[i].pts[j].v, slice[i].pts[j].v);
    }
    pot_ener = energy(&slice[i]);

    kin_fixed = fixed_energy - pot_ener;
    if (kin_fixed < 0) {
      fprintf(sys.filep, "Warning: kinetic energy smaller than zero...\n");
      kin_fixed = 0;
    }
    fac = sqrt(kin_fixed / kin_ener);
    kin_ener = 0;
    for (j = 0; j < sys.npart; j++) {
      slice[i].pts[j].v.x *= fac;
      slice[i].pts[j].v.y *= fac;
      slice[i].pts[j].v.z *= fac;
      kin_ener += 0.5 * vector_inp(slice[i].pts[j].v, slice[i].pts[j].v);
    }
  }
  sys.nfix++;
}

void fix_energy(Slice *psl) {
  double fac, kin_ener, kin_fixed, fixed_energy, pot_ener, pot_tot, kin_tot,
      maxener;

  int i, j;

  fixed_energy = sys.fixed_energy;
  maxener = energy(psl);
  if (maxener > sys.fixed_energy) {
    fprintf(sys.filep, "Warning: maximum energy larger than fixed energy %g\n",
            sys.fixed_energy - maxener);
    fprintf(sys.filep, "rescaling to %g\n", maxener);
    fixed_energy = maxener;
  }

  kin_ener = 0;
  for (j = 0; j < sys.npart; j++) {
    kin_ener += 0.5 * vector_inp(psl->pts[j].v, psl->pts[j].v);
  }
  pot_ener = energy(psl);

  kin_fixed = fixed_energy - pot_ener;
  if (kin_fixed < 0) {
    fprintf(sys.filep, "Warning: kinetic energy smaller than zero...\n");
    kin_fixed = 0;
  }
  fac = sqrt(kin_fixed / kin_ener);
  kin_ener = 0;
  for (j = 0; j < sys.npart; j++) {
    psl->pts[j].v.x *= fac;
    psl->pts[j].v.y *= fac;
    psl->pts[j].v.z *= fac;
    kin_ener += 0.5 * vector_inp(psl->pts[j].v, psl->pts[j].v);
  }
}

void dump_paths(FILE *fp) {
  int irep, j, k, type;
  double sum, x, y, z, a[MAXSTATES];
  Replica *prep;

  printf("in dump\n");
  irep = sys.current_replica;
  prep = replica[irep];
  if (irep > 0) {
    type = analyse(slice, prep, prep->pathlen, sys.initial_state);
    dprint(type);
    if (type == 2) {
      fprintf(fp, "{");
      for (j = 0; j < prep->pathlen; j++) {
        if (j > 0)
          fprintf(fp, ",\n");
        sum = 0;
        for (k = 0; k < sys.nstates; k++)
          sum += 1. / (pow(slice[j].state_rc[k], 2));
        for (k = 0; k < sys.nstates; k++)
          a[k] = 1. / (pow(slice[j].state_rc[k], 2)) / sum;

        x = ((a[0] - 0.25) + 0.5 * (a[1] - 0.25) + 0.5 * (a[3] - 0.25)) /
            0.8165;
        y = ((a[1] - 0.25) + (a[3] - 0.25) / 3.) / 0.9428;
        z = a[3] - 0.25;

        fprintf(fp, "{%7.5lf, %7.5lf, %7.5lf}", x, y, z);
      }
      fprintf(fp, "},\n");
    }
  }
  printf("out dump\n");
  return;
}