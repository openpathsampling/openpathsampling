#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include "path.h"

int check_stability(int);

double get_rc(Slice *psl, int state_index) {

  return psl->state_rc[state_index - 1];
}

double create_all_rc(Slice *psl) {
  vector dr;
  double diff, dr2;
  int i, j;

  for (i = 0; i < sys.npart; i++) {
    dr2 = 0;
    for (j = 0; j < sys.npart; j++)
      if (i != j) {
        vector_minus(psl->pts[i].r, psl->pts[j].r, dr);
        dr2 += vector_inp(dr, dr);
      }
    psl->ordered_dr2[i] = dr2;
  }

  quicksort(psl->ordered_dr2, 0, sys.npart - 1);

  for (j = 0; j < sys.nstates; j++) {
    diff = 0;
    for (i = 0; i < sys.npart; i++) {
      dr2 = psl->ordered_dr2[i] - state[j].target.ordered_dr2[i];
      diff += fabs(dr2);
    }
    psl->state_rc[j] = diff;
  }
  return 0;
}

double print_rc(Slice *psl, int state_index) {

  double dr2;

  dr2 = get_rc(psl, state_index);

  printf("rc = %g\n", dr2);
  return dr2;
}

int in_state(Slice *psl) {
  int in, i;
  double dx, dy, r2;

  for (i = 0; i < sys.nstates; i++) {
    r2 = get_rc(psl, i + 1);
    if (r2 < state[i].min2) {
      return i + 1;
    }
  }
  return 0;
}

int trajectory(Slice *psl, int maxlength) {
  int i, istate;

  for (i = 1; i < maxlength; i++) {
    psl[i] = psl[i - 1];
    propagate_NHL(&psl[i]);
    create_all_rc(&psl[i]);
    if (in_state(&psl[i])) {
      return i;
    }
  }
  printf("maxlength %d reached ", maxlength);
  return 0;
}

int in_upper_window(Slice *psl, double lambda, int state_index) {
  int in;
  double r2new, dx, dy, r2, lambda2;
  r2 = get_rc(psl, state_index);
  in = (r2 > lambda);
  return in;
}

int trajectory_state_i(Slice *psl, Replica *prep, int maxlength,
                       int state_index) {
  int i, istate;

  for (i = 1; i < maxlength; i++) {
    psl[i] = psl[i - 1];
    propagate_NHL(&psl[i]);
    create_all_rc(&psl[i]);
    if (in_upper_window(&psl[i], prep->lambda, state_index))
      return i;
  }
  printf("maxlength %d reached\n", maxlength);
  return 0;
}

int analyse(Slice *psl, Replica *prep, int len, int state_index) {
  int i, initial_state, final_state;

  initial_state = in_state(&psl[0]);
  final_state = in_state(&psl[len - 1]);

  if (initial_state == state_index) {
    if (final_state == state_index) {
      for (i = 0; i < len; i++) {
        if (in_upper_window(&psl[i], prep->lambda, state_index)) {
          return 1;
        }
      }
      return 0;
    } else {
      if (final_state == 0) {
        printf("final state unknown\n");
        return 0;
      }
      for (i = 0; i < len; i++) {
        if (in_upper_window(&psl[i], prep->lambda, state_index)) {
          return 2;
        }
      }
      return 2;
    }
  }

  printf("analyse: path corrupted in replica %d, initial %d  final %d \n",
         prep->index, initial_state, final_state);

  for (i = 0; i < len; i++)
    print_rc(&psl[i], state_index);

  return 0;
}

void analyse_full(Slice *psl, Replica *prep, int len, int state_index) {
  int i, initial_state, final_state;

  initial_state = in_state(&psl[0]);
  final_state = in_state(&psl[len - 1]);

  if (initial_state != state_index) {
    printf("error in path intial state not equal to stateindex\n");
  }

  scr.mstis_mat[initial_state - 1][final_state - 1]++;
  return;
}

int analyse_state_i(Slice *psl, Replica *prep, int len, int state_index) {
  int i, visitedA, visited2;
  double **w;

  visitedA = visited2 = 0;
  if (in_upper_window(&psl[0], prep->lambda, state_index)) {
    for (i = 1; i < len; i++) {
      if (in_state(&psl[i]) == state_index) {
        visitedA = 1;
      }
      if ((visitedA) && (in_upper_window(&psl[i], prep->lambda, state_index))) {
        visited2 = 1;
      }
    }

    if ((visitedA) &&
        (in_upper_window(&psl[len - 1], prep->lambda, state_index))) {
      return 1;
    }
    if ((visitedA) && (visited2)) {
      printf("path starts in 1, visits A and then 1 but doesn't end in 1\n");
      return 1;
    }
  }

  for (i = 0; i < prep->pathlen; i++) {
    print_rc(&slice[i], state_index);
    printf(" lambda = %g ", prep->lambda);
    dprint(in_state(&slice[i]));
  }
  printf("analyse i: path corrupted in replica %d   \n", prep->index);
  return 0;
}

void choose_new_momenta(Slice *psl) {
  Pts *psi;
  vector r_cm, momentum, angular, out, r, v;
  double dv[NPART * 3], vn[NPART * 3], e[NPART * 3], inp, norm;
  int i, j, k, dim;

  for (i = 0; i < sys.npart; i++) {
    psi = &psl->pts[i];
    sys.norm[3][i * 3 + 0] = 0;
    sys.norm[3][i * 3 + 1] = -psi->r.z;
    sys.norm[3][i * 3 + 2] = psi->r.y;
    sys.norm[4][i * 3 + 0] = psi->r.z;
    sys.norm[4][i * 3 + 1] = 0;
    sys.norm[4][i * 3 + 2] = -psi->r.x;
    sys.norm[5][i * 3 + 0] = -psi->r.y;
    sys.norm[5][i * 3 + 1] = psi->r.x;
    sys.norm[5][i * 3 + 2] = 0;
  }

  /*gram-schmidt algorithm*/

  dim = sys.npart * 3;
  for (i = 3; i < 6; i++) {
    for (k = 0; k < dim; k++)
      e[k] = sys.norm[i][k];
    for (j = 0; j < i; j++) {
      inp = 0;
      for (k = 0; k < dim; k++)
        inp += sys.norm[i][k] * sys.norm[j][k];
      for (k = 0; k < dim; k++)
        e[k] -= inp * sys.norm[j][k];
    }
    norm = 0;
    for (k = 0; k < dim; k++)
      norm += e[k] * e[k];
    norm = sqrt(norm);
    for (k = 0; k < dim; k++)
      sys.norm[i][k] = e[k] / norm;
  }

  for (i = 0; i < 6; i++) {
    norm = 0;
    for (k = 0; k < dim; k++)
      norm += sys.norm[i][k] * sys.norm[i][k];
    for (j = 0; j < i; j++) {
      inp = 0;
      for (k = 0; k < dim; k++)
        inp += sys.norm[i][k] * sys.norm[j][k];
    }
  }

  for (k = 0; k < dim; k++)
    vn[k] = dv[k] = gssran() * sys.dvmax;
  for (i = 0; i < 6; i++) {
    inp = 0;
    for (k = 0; k < dim; k++)
      inp += sys.norm[i][k] * dv[k];
    for (k = 0; k < dim; k++)
      vn[k] -= inp * sys.norm[i][k];
  }

  k = 0;
  for (i = 0; i < sys.npart; i++) {
    psl->pts[i].v.x += vn[k++];
    psl->pts[i].v.y += vn[k++];
    psl->pts[i].v.z += vn[k++];
  }

  return;
}

int shoot(Replica *prep) {
  vector v;
  int i, j, index, fwlen, bwlen, pathlen, maxlength, type, statenotfound,
      initial, final, chosenrandomindex = 0;
  double dr2, kin_ener, pot_ener, old_kin_ener, fac, en_diff, aux;

  if (prep->index == 0)
    return 0;

  index = 0;

  for (i = 0; i < prep->pathlen; i++) {
    if (in_upper_window(&slice[i], prep->lambda, sys.initial_state)) {
      index = i;
      i = prep->pathlen;
    }
  }

  if (in_upper_window(&slice[index], prep->lambda, sys.initial_state) == 0) {
    index = 1 + (int)(ran3() * (prep->pathlen - 1));
    chosenrandomindex = 1;
  }

  index = 1 + (int)(ran3() * (prep->pathlen - 1));
  chosenrandomindex = 1;

  type = analyse(slice, prep, prep->pathlen, sys.initial_state);
  if (type == 0) {
    printf("error: type 0 is not correct\n");
    dprint(sys.initial_state);
    dprint(prep->pathlen);
    dprint(prep->index);
    gprint(prep->lambda);
    for (i = 0; i < prep->pathlen; i++) {
      print_rc(&slice[i], sys.initial_state);
    }
    exit(1);
  }

  fwtrial[0] = slice[index];

  old_kin_ener = 0;
  for (i = 0; i < sys.npart; i++) {
    v = fwtrial[0].pts[i].v;
    old_kin_ener += 0.5 * vector_inp(v, v);
  }

  choose_new_momenta(&fwtrial[0]);

  kin_ener = 0;
  for (i = 0; i < sys.npart; i++) {
    v = fwtrial[0].pts[i].v;
    kin_ener += 0.5 * vector_inp(v, v);
  }

  pot_ener = energy(&fwtrial[0]);
  if (sys.ensemble_type == NVT) {
    fac = sqrt(old_kin_ener / kin_ener);
  } else {
    fac = sqrt((sys.fixed_energy - pot_ener) / kin_ener);
  }
  for (i = 0; i < sys.npart; i++) {
    fwtrial[0].pts[i].v.x *= fac;
    fwtrial[0].pts[i].v.y *= fac;
    fwtrial[0].pts[i].v.z *= fac;
  }

  bwtrial[0] = fwtrial[0];

  for (i = 0; i < sys.npart; i++) {
    bwtrial[0].pts[i].v.x *= -1;
    bwtrial[0].pts[i].v.y *= -1;
    bwtrial[0].pts[i].v.z *= -1;
  }
  bwtrial[0].xi = -bwtrial[0].xi;

  maxlength = MAXPATH;
  bwlen = trajectory(&bwtrial[0], maxlength);
  if (bwlen == 0) {
    printf("rejected because bwlen=0, too long\n");
    return 0;
  }
  initial = in_state(&bwtrial[bwlen]);
  if (initial != sys.initial_state) {

    return 0;
  }

  for (i = 0; i <= bwlen; i++) {
    for (j = 0; j < sys.npart; j++) {
      bwtrial[i].pts[j].v.x *= -1;
      bwtrial[i].pts[j].v.y *= -1;
      bwtrial[i].pts[j].v.z *= -1;
    }
    bwtrial[i].xi = -bwtrial[i].xi;
  }

  maxlength = MAXPATH - bwlen - 1;
  fwlen = trajectory(&fwtrial[0], maxlength);
  if (fwlen == 0) {

    printf("rejected because fwlen=0, too long, maxlenght = %d\n", maxlength);

    printf("this trial trajectory is a candidate for new state\n");
    find_minimum(&fwtrial[maxlength - 1]);

    statenotfound = 1;
    for (i = 1; i <= sys.nstates; i++) {
      dr2 = get_rc(&fwtrial[maxlength - 1], i);
      printf("distance to state %d = %g\n", i, dr2);
      if (dr2 < 0.1)
        statenotfound = 0;
    }
    if (statenotfound) {
      if (add_newstate(&fwtrial[maxlength - 1])) {
        if (statenotfound == 0) {
          printf("the state was found in the database, there is something "
                 "wrong\n");
          return 0;
        } else {
          printf("new state was found and added to the database, might as well "
                 "accept path\n");
        }
        printf("but we don't for safety, check path for new state\n");

        for (i = 0; i < prep->pathlen; i++) {
          create_all_rc(&slice[i]);
          j = in_state(&slice[i]);
          if ((j != sys.initial_state) && (j != 0)) {
            printf("found new state in old trajectory at slice %d. Breaking "
                   "trajectory. \n",
                   i);
            prep->pathlen = i + 1;
            dprint(prep->pathlen);
            path.nslices = prep->pathlen;
            prep->type = analyse(slice, prep, prep->pathlen, sys.initial_state);
            printf("old trajectory now type %d\n", prep->type);
            return 0;
          }
        }
      }
    }
    return 0;
  }

  final = in_state(&fwtrial[fwlen]);
  if (final == 0) {
    printf("forward path rejected because not reached stat\n");
    return 0;
  }

  en_diff = total_energy(&fwtrial[fwlen]) - total_energy(&bwtrial[bwlen]);
  if ((fabs(en_diff) > 20) || (en_diff / en_diff != 1)) {
    printf("shoot: failure in integration detected, de =%g\n", en_diff);
    return 0;
  }

  pathlen = bwlen + fwlen + 1;
  if (chosenrandomindex) {
    aux = (double)(prep->pathlen - 2) / (pathlen - 2);
    if (aux < ran3()) {
      return 0;
    }
  }

  for (i = 0; i < bwlen; i++)
    trial[i] = bwtrial[bwlen - i];
  for (i = 0; i <= fwlen; i++)
    trial[i + bwlen] = fwtrial[i];

  type = analyse(trial, prep, pathlen, sys.initial_state);
  if (type == 0) {
    return 0;
  }

  for (i = 0; i < pathlen; i++)
    slice[i] = trial[i];

  prep->pathlen = pathlen;
  path.nslices = pathlen;
  prep->type = type;

  return 1;
}

int shoot_oneway(Replica *prep) {
  int i, j, index, len, pathlen, maxlength, trial_pathlen, type, statenotfound;
  double dr2;
  printf("in shoot state %d rep %d lambda %g\n", sys.initial_state, prep->index,
         prep->lambda);

  if (prep->index == 0)
    return 0;

  index = 0;
  for (i = 0; i < prep->pathlen; i++) {
    if (in_upper_window(&slice[i], prep->lambda, sys.initial_state)) {
      index = i;
      i = prep->pathlen;
    }
  }

  if (in_upper_window(&slice[index], prep->lambda, sys.initial_state) == 0) {
    printf("error: path not reached lambda %g\n", prep->lambda);
    dprint(index);
    dprint(sys.initial_state);
    dprint(prep->index);
    return 0;
  }

  for (i = 0; i <= index; i++)
    trial[i] = slice[i];
  maxlength = MAXPATH - index;
  len = trajectory(&trial[index], maxlength);
  if (len == 0) {
    printf("rejected because len=0, too long\n");
    printf("this trial trajectory is a candidate for new state\n");
    find_minimum(&trial[maxlength - 1]);

    statenotfound = 1;
    for (i = 1; i <= sys.nstates; i++) {
      dr2 = get_rc(&trial[maxlength - 1], i);
      printf("distance to state %d = %g\n", i, dr2);
      if (dr2 < 0.1)
        statenotfound = 0;
    }
    if (statenotfound)
      add_newstate(&trial[maxlength - 1]);
    find_minimum(&trial[0]);
    if (statenotfound == 0) {
      printf("the state was found in the database, there is something wrong\n");
      return 0;
    } else {
      printf("new state was found and added to the database, might as well "
             "accept path\n");
    }
  }
  trial_pathlen = index + len + 1;
  if (trial_pathlen > MAXPATH) {
    printf("error: pathlen = %d\n", trial_pathlen);
    return 0;
  }

  type = analyse(trial, prep, trial_pathlen, sys.initial_state);
  if (type == 0) {
    printf("trial path rejected\n");

    return 0;
  }

  for (i = 0; i < trial_pathlen; i++)
    slice[i] = trial[i];

  prep->pathlen = trial_pathlen;
  path.nslices = trial_pathlen;
  prep->type = type;

  i = in_state(&slice[0]);
  j = in_state(&slice[trial_pathlen - 1]);

  printf("path accepted state= %d, replica =%d , length = %d type = %d, "
         "initial= %d, final = %d\n",
         sys.initial_state, prep->index, trial_pathlen, type, i, j);

  return 1;
}

int swap_replica(int irep, int jrep) {
  Replica *prepi, *prepj;
  Slice *swaptrial, h;
  int i, pathlen, ABexchOK, ireptype, jreptype, type;
  double aux;

  prepi = replica[irep];
  prepj = replica[jrep];

  prepi->dos += sys.scalefactor;

  if ((jrep < 0))
    return 0;

  if ((irep == 0) || (jrep == 0)) {
    return swap_replica_0(irep, jrep);
  }
  if ((irep == sys.nreplica - 1) && (jrep == sys.nreplica))
    return 0;
  if (jrep > sys.nreplica - 1)
    return 0;

  aux = prepi->dos - prepj->dos;

  if (ran3() > exp(aux)) {
    return 0;
  }

  type = analyse(slice, prepj, prepi->pathlen, sys.initial_state);
  if (type == 0) {
    return 0;
  }

  prepj->pathlen = prepi->pathlen;
  prepj->type = prepi->type;
  sys.current_replica = jrep;

  return 1;
}

int swap_replica_0(int irep, int jrep) {
  Replica *prep;
  Slice *trial0;
  Slice *trial1;
  int i, j, pathlen0, pathlen1, inA, in1, index, maxlength, start, len, type,
      type0, type1, reversal;

  trial0 = fwtrial;
  trial1 = bwtrial;

  if (irep == 0) {

    prep = replica[0];
    type = analyse_state_i(slice, prep, prep->pathlen, sys.initial_state);
    if (prep->type != type) {
      printf("error : type %d and replica type %d do not match\n", type,
             prep->type);
      return 0;
    }
    if (type == 0) {
      printf("error: replica 0 has wrong type\n");
      return 0;
    }

    reversal = (ran3() < 0.5) ? 1 : 0;
    if (reversal) {
      pathlen0 = prep->pathlen;
      for (j = 0; j < pathlen0; j++) {
        trial1[j] = slice[pathlen0 - 1 - j]; // use trial1 for temp storage, to
                                             // not affect slice
        for (i = 0; i < sys.npart; i++) {
          trial1[j].pts[i].v.x *= -1;
          trial1[j].pts[i].v.y *= -1;
          trial1[j].pts[i].v.z *= -1;
        }
        trial1[j].xi = -trial1[j].xi;
      }
    } else {
      for (i = 0; i < prep->pathlen; i++) {
        trial1[i] =
            slice[i]; // use trial1 for temp storage, to not affect slice
      }
    }

    for (i = prep->pathlen - 1; i >= 0; i--) {
      if (in_state(&trial1[i])) {
        start = i;
        break;
      }
    }
    for (i = start; i < prep->pathlen; i++)
      trial0[i - start] = trial1[i];
    index = prep->pathlen - start - 1;
    maxlength = MAXPATH - index;
    len = trajectory(&trial0[index], maxlength);
    if (len == 0)
      return 0;
    pathlen0 = index + len + 1;
    if (pathlen0 > MAXPATH) {
      printf("error: pathlen = %d\n", pathlen0);
      return 0;
    }
    type = analyse(trial0, replica[1], pathlen0, sys.initial_state);
    if (type == 0)
      return 0;
    type0 = type;
    for (i = 0; i < pathlen0; i++)
      slice[i] = trial0[i];
    replica[1]->pathlen = pathlen0;
    replica[1]->type = type0;

    sys.current_replica = jrep;
    path.nslices = pathlen0;
    prep = replica[1];
    type = analyse(slice, prep, prep->pathlen, sys.initial_state);
    if (prep->type != type) {
      printf("swap 0->1 : type %d  and replica type %d do not match\n", type,
             prep->type);
      return 0;
    }
    if (type == 0) {
      printf("swap 0->1: replica 0 has wrong type\n");
      return 0;
    }
  }

  if (jrep == 0) {
    prep = replica[1];
    type = analyse(slice, prep, prep->pathlen, sys.initial_state);
    if (prep->type != type) {
      printf("error : type %d and replica type %d do not match\n", type,
             prep->type);
      return 0;
    };
    if (type == 0) {
      printf("reject: replica 1 has wrong type\n");
      return 0;
    }

    reversal = (ran3() < 0.5) ? 1 : 0;
    if ((type == 2) && (reversal == 0)) {
      printf("type is 2, and reversal =0, reject\n");
      return 0;
    }

    if (reversal) {
      pathlen1 = prep->pathlen;
      for (j = 0; j < pathlen1; j++) {
        trial0[j] = slice[pathlen1 - 1 - j]; // use trial0 for temp storage
        for (i = 0; i < sys.npart; i++) {
          trial0[j].pts[i].v.x *= -1;
          trial0[j].pts[i].v.y *= -1;
          trial0[j].pts[i].v.z *= -1;
        }
        trial0[j].xi = -trial0[j].xi;
      }
    } else {
      for (i = 0; i < prep->pathlen; i++) {
        trial0[i] =
            slice[i]; // use trial0 for temp storage, to not affect slice
      }
    }

    start = -1;
    for (i = prep->pathlen - 1; i >= 0; i--) {
      if (in_upper_window(&trial0[i], replica[0]->lambda, sys.initial_state)) {
        start = i;
        break;
      }
    }
    if (start == -1) {
      printf("path type 2 does not pass lambda %g\n", replica[0]->lambda);
      return 0;
    }

    for (i = start; i < prep->pathlen; i++)
      trial1[i - start] = trial0[i];
    index = prep->pathlen - start - 1;

    maxlength = MAXPATH - index;
    len =
        trajectory_state_i(&trial1[index], prep, maxlength, sys.initial_state);
    if (len == 0)
      return 0;
    pathlen1 = index + len + 1;

    if (pathlen1 > MAXPATH) {
      printf("error: pathlen = %d\n", pathlen1);
      return 0;
    }

    type = analyse_state_i(trial1, replica[0], pathlen1, sys.initial_state);
    if (type == 0) {
      return 0;
    }
    type1 = type;
    for (i = 0; i < pathlen1; i++)
      slice[i] = trial1[i];
    replica[0]->pathlen = pathlen1;
    replica[0]->type = type1;
    sys.current_replica = jrep;
    path.nslices = pathlen1;

    prep = replica[0];

    type = analyse_state_i(slice, prep, prep->pathlen, sys.initial_state);
    // dprint(type);
    if (prep->type != type) {
      printf("swap 1->0 : type %d and replica type %d do not match\n", type,
             prep->type);
      return 0;
    }
    if (type == 0) {
      printf("swap 1->0: replica 0 has wrong type\n");
      return 0;
    }
  }

  SWAP(replica[irep]->swapindex, replica[jrep]->swapindex, i);
  return 1;
}

int swap_replica_n(int irep) {
  Replica *prep;
  int i, j, initial_state, final_state, pathlen, type, start, index, maxlength,
      len;

  prep = replica[irep];
  if (prep->type == 1)
    return 0;
  pathlen = prep->pathlen;

  initial_state = in_state(&slice[0]);
  final_state = in_state(&slice[pathlen - 1]);

  for (i = 0; i < pathlen; i++) {
    trial[i] = slice[pathlen - 1 - i];
  }

  type = analyse(trial, prep, pathlen, final_state);
  if (type == 0)
    return 0;

  printf("accepting swap n, new initial state = %d\n", final_state);

  sys.initial_state = final_state;

  for (i = 0; i < sys.nreplica; i++) {
    replica[i] = &state[sys.initial_state - 1].srep[i];
  }
  prep = replica[irep];

  prep->pathlen = pathlen;
  for (j = 0; j < pathlen; j++) {
    for (i = 0; i < sys.npart; i++) {
      trial[j].pts[i].v.x *= -1;
      trial[j].pts[i].v.y *= -1;
      trial[j].pts[i].v.z *= -1;
    }
    trial[j].xi = -trial[j].xi;
    slice[j] = trial[j];
  }

  type = analyse(slice, prep, prep->pathlen, sys.initial_state);
  printf("new initital state %d replica %d  has type %d \n", sys.initial_state,
         irep, type);
  prep->type = type;

  return 1;
}

int swap_states(int irep, int jrep) {
  Replica *prepi, *prepj;
  int i, j, initial_state, final_state, pathlen, type, start, index, maxlength,
      len;
  double aux;

  prepi = replica[irep];

  if (prepi->type == 1)
    return 0;
  pathlen = prepi->pathlen;

  initial_state = in_state(&slice[0]);
  final_state = in_state(&slice[pathlen - 1]);

  if (initial_state != sys.initial_state) {
    printf("error: initial state not correct\n");
    // exit(1);
    return 0;
  }
  if ((final_state < 1) || (final_state > MAXSTATES)) {
    printf("warning: final slice not ending in known state, rejecting\n");
    return 0;
  }

  prepj = &state[final_state - 1].srep[jrep];

  aux = prepi->dos - prepj->dos;

  if (ran3() > exp(aux)) {
    return 0;
  }

  for (i = 0; i < pathlen; i++) {
    trial[i] = slice[pathlen - 1 - i];
  }

  type = analyse(trial, prepj, pathlen, final_state);
  if (type == 0)
    return 0;

  sys.initial_state = final_state;

  for (i = 0; i < sys.nreplica; i++) {
    replica[i] = &state[sys.initial_state - 1].srep[i];
  }

  prepj->pathlen = pathlen;
  for (j = 0; j < pathlen; j++) {
    for (i = 0; i < sys.npart; i++) {
      trial[j].pts[i].v.x *= -1;
      trial[j].pts[i].v.y *= -1;
      trial[j].pts[i].v.z *= -1;
    }
    trial[j].xi = -trial[j].xi;
    slice[j] = trial[j];
  }

  type = analyse(slice, prepj, prepj->pathlen, sys.initial_state);
  prepj->type = type;
  sys.current_replica = jrep;

  return 1;
}

int reverse_replica(int irep) {
  Replica *prep;
  int i, j, pathlen, type;

  if (irep == 0)
    return 0;
  prep = replica[irep];

  type = analyse(slice, prep, prep->pathlen, sys.initial_state);
  if (prep->type != type)
    printf("error: type %d and replica type %d do not match\n", type,
           prep->type);
  if (prep->type != 1)
    return 0;

  pathlen = prep->pathlen;
  for (i = 0; i < pathlen; i++) {
    trial[i] = slice[pathlen - 1 - i];
  }

  type = analyse(trial, prep, pathlen, sys.initial_state);
  if (type == 0)
    return 0;
  if (type == 2)
    return 0;

  for (j = 0; j < pathlen; j++) {
    for (i = 0; i < sys.npart; i++) {
      trial[j].pts[i].v.x *= -1;
      trial[j].pts[i].v.y *= -1;
      trial[j].pts[i].v.z *= -1;
    }
    trial[j].xi = -trial[j].xi;
    slice[j] = trial[j];
  }

  return 1;
}

void save_trialpath(int ifile, Slice *psl) {}

int add_newstate(Slice *psl) {

  vector r_cm, dr;
  double dr2;
  int i, j, k;

  i = sys.nstates;
  state[i].target = *psl;

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

  for (j = 0; j < sys.npart; j++) {
    vprint(state[i].target.pts[j].r);
  }

  state[i].min = sys.min_stable;
  state[i].min2 = sys.min_stable;

  if (sys.nstates < MAXSTATES) {
    sys.nstates++;
    if (check_stability(i) == 0) {
      printf("proposed state not stable enough, rejecting\n");
      sys.nstates--;
      return 0;
    } else {
      printf("added state %d to database\n", i);
    }
  }

  printf("number of states = %d\n", sys.nstates);

  return 1;
}

int check_stability(int istate) {
  vector v;
  double kin_ener, pot_ener, fac;
  int i, j, relaxed2other, trial_state;

  trial[0] = state[istate].target;
  trial_state = istate + 1;
  dprint(istate);
  dprint(trial_state);

  choose_new_momenta(&trial[0]);

  kin_ener = 0;
  for (i = 0; i < sys.npart; i++) {
    v = trial[0].pts[i].v;
    kin_ener += 0.5 * vector_inp(v, v);
  }

  gprint(kin_ener);
  pot_ener = energy(&trial[0]);
  fac = sqrt((sys.fixed_energy - pot_ener) / kin_ener);
  fac = 0;
  gprint(fac);
  for (i = 0; i < sys.npart; i++) {
    trial[0].pts[i].v.x *= fac;
    trial[0].pts[i].v.y *= fac;
    trial[0].pts[i].v.z *= fac;
  }

  relaxed2other = 0;
  for (i = 1; i < MAXPATH; i++) {
    trial[i] = trial[i - 1];
    propagate_NHL(&trial[i]);
    create_all_rc(&trial[i]);
    j = in_state(&trial[i]);
    if ((j != trial_state) && (j != 0)) {
      relaxed2other = 1;
    }
  }

  dprint(relaxed2other);

  if (relaxed2other) {
    printf("state not very stable\n");
    return 0;
  } else {
    printf("state stable\n");
    return 1;
  }
}