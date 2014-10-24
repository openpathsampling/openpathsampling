#include <stdio.h>
#include <math.h>
#include "path.h"

int check_endpoint(Slice *);
int check_beginpoint(Slice *);
vector force(Slice *, int ipart);
void lowe_andersen(Slice *);
double dist(vector *, vector *);
double log_prob(int, vector *, vector *, vector *);
int beadsinA(void);
int beadsinB(void);

double total_energy(Slice *psl) {
  Pts *psi, *psj;
  vector rij;
  double r2, r6, kin_ener = 0, pot_ener = 0;
  int i, j;

  for (i = 0; i < sys.npart; i++) {
    psi = &psl->pts[i];
    for (j = 0; j < i; j++) {
      psj = &psl->pts[j];
      vector_minus(psi->r, psj->r, rij);
      r2 = vector_inp(rij, rij);
      r6 = r2 * r2 * r2;
      pot_ener += 4 * (1. / (r6 * r6) - 1. / r6);
    }
    kin_ener += vector_inp(psi->v, psi->v);
  }

  return pot_ener + 0.5 * kin_ener;
}

double energy(Slice *psl) {
  Pts *psi, *psj;
  vector rij;
  double r2, r6, en = 0;
  int i, j;

  for (i = 0; i < sys.npart; i++) {
    psi = &psl->pts[i];
    for (j = 0; j < i; j++) {
      psj = &psl->pts[j];
      vector_minus(psi->r, psj->r, rij);
      r2 = vector_inp(rij, rij);
      r6 = r2 * r2 * r2;
      en += 4 * (1. / (r6 * r6) - 1. / r6);
    }
  }
  return en;
}

vector force(Slice *psl, int ipart) {
  Pts *psi, *psj;
  vector rij, f;
  double r2, r6, rt, r, fmag, fac;
  int i;

  f = nulvec;
  psi = &psl->pts[ipart];
  for (i = 0; i < sys.npart; i++) {
    if (i != ipart) {
      psj = &psl->pts[i];
      vector_minus(psi->r, psj->r, rij);

      r2 = vector_inp(rij, rij);
      r6 = r2 * r2 * r2;
      fmag = -12. / (r6 * r6 * r2) + 6. / (r6 * r2);
      fac = -4 * fmag;
      scalar_plustimes(rij, fac, f);
    }
  }

  return f;
}

void calculate_forces(Slice *psl) {
  Pts *psi, *psj;
  vector rij;
  double r2, r6, fmag;
  int i, j;

  for (i = 0; i < sys.npart; i++)
    psl->pts[i].f = nulvec;

  for (i = 0; i < sys.npart; i++) {
    psi = &psl->pts[i];
    for (j = 0; j < i; j++) {
      psj = &psl->pts[j];
      vector_minus(psi->r, psj->r, rij);
      r2 = vector_inp(rij, rij);
      r6 = r2 * r2 * r2;
      fmag = 48. / (r6 * r6 * r2) - 24. / (r6 * r2);
      scalar_plustimes(rij, fmag, psi->f);
      scalar_mintimes(rij, fmag, psj->f);
    }
  }

  return;
}

void propagate_NHL(Slice *psl) {
  Pts *psi;
  vector dr, dv;
  double dt, dt2, expxi, en_kin2;
  int i, j, k;

  dt = sys.delta_t;
  dt2 = 0.5 * dt;

  for (k = 0; k < path.ninter; k++) {

    psl->xi = langevin.c1 * psl->xi + langevin.c3 * gssran();
    en_kin2 = 0;
    for (i = 0; i < sys.npart; i++) {
      psi = &psl->pts[i];

      expxi = exp(-psl->xi * dt2);
      psi->v.x *= expxi;
      psi->v.y *= expxi;
      psi->v.z *= expxi;

      psi->v.x += dt2 * psi->f.x;
      psi->v.y += dt2 * psi->f.y;
      psi->v.z += dt2 * psi->f.z;

      en_kin2 += vector_inp(psi->v, psi->v);

      psi->r.x += dt2 * psi->v.x;
      psi->r.y += dt2 * psi->v.y;
      psi->r.z += dt2 * psi->v.z;
    }

    psl->xi += dt * (en_kin2 - (3 * sys.npart - 3) / sys.beta) / langevin.mu;
    psl->s += dt * psl->xi;

    for (i = 0; i < sys.npart; i++) {
      psi = &psl->pts[i];
      psi->r.x += dt2 * psi->v.x;
      psi->r.y += dt2 * psi->v.y;
      psi->r.z += dt2 * psi->v.z;
    }

    calculate_forces(psl);

    for (i = 0; i < sys.npart; i++) {
      psi = &psl->pts[i];
      psi->v.x += dt2 * psi->f.x;
      psi->v.y += dt2 * psi->f.y;
      psi->v.z += dt2 * psi->f.z;

      expxi = exp(-psl->xi * dt2);
      psi->v.x *= expxi;
      psi->v.y *= expxi;
      psi->v.z *= expxi;
    }

    psl->xi = langevin.c1 * psl->xi + langevin.c3 * gssran();
  }
}

void propagate(Slice *psl) {
  Pts *psi;
  vector dr, dv;
  double dt, dt2;
  int i, j, k;

  dt = sys.delta_t;
  dt2 = dt * 0.5;
  for (k = 0; k < path.ninter; k++) {
    for (i = 0; i < sys.npart; i++) {
      psi = &psl->pts[i];
      psi->r.x += dt * (psi->v.x + dt2 * psi->f.x);
      psi->r.y += dt * (psi->v.y + dt2 * psi->f.y);
      psi->r.z += dt * (psi->v.z + dt2 * psi->f.z);
      psi->v.x = psi->v.x + dt2 * psi->f.x;
      psi->v.y = psi->v.y + dt2 * psi->f.y;
      psi->v.z = psi->v.z + dt2 * psi->f.z;
    }
    calculate_forces(psl);
    for (i = 0; i < sys.npart; i++) {
      psi = &psl->pts[i];
      psi->v.x += dt2 * psi->f.x;
      psi->v.y += dt2 * psi->f.y;
      psi->v.z += dt2 * psi->f.z;
    }
    if (sys.ensemble_type == NVT)
      if (ran3() < 10 * dt) {
        lowe_andersen(psl);
      }
  }
}

void lowe_andersen(Slice *psl) {

  vector dr, dv;
  double l, dvdr, lagauss, mag;
  int i, j;

  do {
    i = (int)(ran3() * sys.npart);
    j = (int)(ran3() * sys.npart);
  } while (i == j);

  vector_minus(psl->pts[i].v, psl->pts[j].v, dv);
  vector_minus(psl->pts[i].r, psl->pts[j].r, dr);
  l = sqrt(vector_inp(dr, dr));
  scalar_divide(dr, l, dr);
  dvdr = vector_inp(dv, dr);
  lagauss = gssran() * sys.sqrttemp;
  mag = 0.5 * (lagauss - dvdr);
  scalar_plustimes(dr, mag, psl->pts[i].v);
  scalar_mintimes(dr, mag, psl->pts[j].v);
  return;
}

void check_energy(Slice *psl) {
  vector r_cm, momentum, angular, out, v, r;
  double kin_ener, pot_en, NHL_term, NHL_tot, temp;
  int i, j, inext;

  kin_ener = 0;
  pot_en = energy(psl);
  for (j = 0; j < sys.npart; j++) {
    kin_ener += 0.5 * vector_inp(psl->pts[j].v, psl->pts[j].v);
  }

  NHL_term = psl->xi * psl->xi * langevin.mu / 2. +
             (3 * sys.npart - 3) * psl->s / sys.beta;
  NHL_tot = kin_ener + pot_en + NHL_term;

  temp = (2 * kin_ener) / (3 * sys.npart - 3);

  fprintf(sys.filep, "E = %9.2lf, K = %9.2lf, E+K = %9.2lf  NHL = %9.2lf , "
                     "E+K+NHL= %9.2lf T= %9.2lf\n",
          pot_en, kin_ener, kin_ener + pot_en, NHL_term, NHL_tot, temp);

  r_cm = nulvec;
  for (j = 0; j < sys.npart; j++) {
    vector_add(r_cm, psl->pts[j].r, r_cm);
  }
  scalar_divide(r_cm, sys.npart, r_cm);
  vprint(r_cm);

  angular = momentum = nulvec;
  for (j = 0; j < sys.npart; j++) {
    v = psl->pts[j].v;
    r = psl->pts[j].r;

    vector_minus(r, r_cm, r);
    vector_add(momentum, v, momentum);
    vector_cross(r, v, out);
    vector_add(angular, out, angular);
  }
  vprint(momentum);
  vprint(angular);
  print_rc(psl, sys.initial_state);
  return;
}

void propagate_zerovelocity(Slice *psl) {
  Pts *psi;
  vector dr, dv;
  double dt, dt2;
  int i, j, k;

  dt = sys.delta_t;
  dt2 = dt * 0.5;
  for (k = 0; k < path.ninter; k++) {
    for (i = 0; i < sys.npart; i++) {
      psi = &psl->pts[i];
      psi->r.x += dt * (psi->v.x + dt2 * psi->f.x);
      psi->r.y += dt * (psi->v.y + dt2 * psi->f.y);
      psi->r.z += dt * (psi->v.z + dt2 * psi->f.z);
      psi->v.x = psi->v.x + dt2 * psi->f.x;
      psi->v.y = psi->v.y + dt2 * psi->f.y;
      psi->v.z = psi->v.z + dt2 * psi->f.z;
    }
    calculate_forces(psl);
    for (i = 0; i < sys.npart; i++) {
      psi = &psl->pts[i];
      psi->v.x += dt2 * psi->f.x;
      psi->v.y += dt2 * psi->f.y;
      psi->v.z += dt2 * psi->f.z;
    }
  }
}

void find_minimum(Slice *psl) {
  Slice rold;
  vector dr;
  double dr2;
  int i, j;

  do {
    rold = *psl;
    for (i = 0; i < 100; i++) {
      for (j = 0; j < sys.npart; j++) {
        psl->pts[j].v = nulvec;
      }
      propagate_zerovelocity(psl);
    }
    dr2 = 0;
    for (j = 0; j < sys.npart; j++) {

      vector_minus(psl->pts[j].r, rold.pts[j].r, dr);
      dr2 += vector_inp(dr, dr);
    }
  } while (dr2 > 1e-15);

  gprint(dr2);

  return;
}
