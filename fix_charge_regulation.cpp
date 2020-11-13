/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Tine Curk and Jiaxing Yuan
------------------------------------------------------------------------- */
#include "fix_charge_regulation.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "comm.h"
#include "compute.h"
#include "group.h"
#include "domain.h"
#include "region.h"
#include "random_park.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "neighbor.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

// large energy value used to signal overlap
#define MAXENERGYSIGNAL 1.0e100
#define MAXENERGYTEST 1.0e50
#define small 0.0000001
#define PI 3.1415926

enum {
    EXCHATOM, EXCHMOL
}; // exchmode
enum {
    MOVEATOM, MOVEMOL
}; // movemode

/* ---------------------------------------------------------------------- */
/* fix ID group-ID charge_regulation nevery nexchanges lb pKa pKb pH pI^+ pI^- reservoir_temperature seed acid_type cation_type base_type anion_type reaction_cutoff*/
/*If reaction_cutoff is set to zero, it means there is no geometry constraint when inserting and deleting an atom*/
Fix_charge_regulation::Fix_charge_regulation(LAMMPS *lmp, int narg, char **arg) :
        Fix(lmp, narg, arg),
        idregion(NULL), full_flag(0), ngroups(0), groupstrings(NULL), ngrouptypes(0), grouptypestrings(NULL),
        grouptypebits(NULL), grouptypes(NULL), local_gas_list(NULL), molcoords(NULL), molq(NULL), molimage(NULL),
        random_equal(NULL), random_unequal(NULL),
        idftemp(NULL), ptype_ID(NULL) {

    dynamic_group_allow = 1;
    vector_flag = 1;
    size_vector = 8;
    global_freq = 1;
    extvector = 0;
    restart_global = 1;
    time_depend = 1;
    gcmc_nmax = 0;
    regionflag = 0;
    overlap_flag = 0;
    energy_stored = 0;

    pH = force->numeric(FLERR, arg[3]);
    pKa = force->numeric(FLERR, arg[4]);
    pKb = force->numeric(FLERR, arg[5]);

    pI_plus = force->numeric(FLERR, arg[6]);
    pI_minus = force->numeric(FLERR, arg[7]);

    acid_type = force->inumeric(FLERR, arg[8]);
    base_type = force->inumeric(FLERR, arg[9]);
    cation_type = force->inumeric(FLERR, arg[10]);
    anion_type = force->inumeric(FLERR, arg[11]);

    // read optional arguments
    options(narg-12,&arg[12]);

    if (nevery <= 0) error->all(FLERR, "Illegal fix charge_regulation command");
    if (nexchanges < 0) error->all(FLERR, "Illegal fix charge_regulation command");
    if (lb < 0.0) error->all(FLERR, "Illegal fix charge_regulation command");
    if (pH < 0.0) error->all(FLERR, "Illegal fix charge_regulation command");
    if (pH < 0.0) error->all(FLERR, "Illegal fix charge_regulation command");
    if (pKa < 0.0) error->all(FLERR, "Illegal fix charge_regulation command");
    if (pKb < 0.0) error->all(FLERR, "Illegal fix charge_regulation command");
    if (pI_plus < 0.0) error->all(FLERR, "Illegal fix charge_regulation command");
    if (pI_minus < 0.0) error->all(FLERR, "Illegal fix charge_regulation command");
    if (*target_temperature_tcp < 0.0) error->all(FLERR, "Illegal fix charge_regulation command");
    if (seed <= 0) error->all(FLERR, "Illegal fix charge_regulation command");
  //  if (acid_type <= 0) error->all(FLERR, "Illegal fix charge_regulation command");
  //  if (cation_type <= 0) error->all(FLERR, "Illegal fix charge_regulation command");
  //  if (base_type <= 0) error->all(FLERR, "Illegal fix charge_regulation command");
  //  if (anion_type <= 0) error->all(FLERR, "Illegal fix charge_regulation command");
    if (reaction_distance < 0.0) error->all(FLERR, "Illegal fix charge_regulation command");

    force_reneighbor = 1;
    next_reneighbor = update->ntimestep + 1;
    random_equal = new RanPark(lmp, seed);
    random_unequal = new RanPark(lmp, seed);
    nacid_attempts = 0.0;
    nacid_successes = 0.0;
    nbase_attempts = 0.0;
    nbase_successes = 0.0;
    nsalt_attempts = 0.0;
    nsalt_successes = 0.0;
    exclusion_group_bit = 0.0;
}

int Fix_charge_regulation::setmask() {
    int mask = 0;
    mask |= PRE_EXCHANGE;
    return mask;
}

void Fix_charge_regulation::init() {

    triclinic = domain->triclinic;

    char *id_pe = (char *) "thermo_pe";
    int ipe = modify->find_compute(id_pe);
    c_pe = modify->compute[ipe];


    int *type = atom->type;

    if (atom->molecule_flag) {
        tagint *molecule = atom->molecule;
        int flag = 0;
        for (int i = 0; i < atom->nlocal; i++)
            if (type[i] == ngcmc_type)
                if (molecule[i]) flag = 1;
        int flagall;
        MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);
        if (flagall && comm->me == 0)
            error->all(FLERR,
                       "Fix charge regulation cannot exchange individual atoms belonging to a molecule");
    }

    if (domain->dimension == 2)
        error->all(FLERR, "Cannot use fix charge regulation in a 2d simulation");

    gas_mass = atom->mass[cation_type];
    if (gas_mass <= 0.0)
        error->all(FLERR, "Illegal fix charge regulation gas mass <= 0");

}

void Fix_charge_regulation::pre_exchange() {

    if (next_reneighbor != update->ntimestep) return;
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
    if (triclinic) {
        sublo = domain->sublo_lamda;
        subhi = domain->subhi_lamda;
    } else {
        sublo = domain->sublo;
        subhi = domain->subhi;
    }
    if (regionflag) volume = region_volume;
    else volume = domain->xprd * domain->yprd * domain->zprd;
    if (triclinic) domain->x2lamda(atom->nlocal);


    domain->pbc();
    comm->exchange();
    atom->nghost = 0;
    comm->borders();


    if (triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
    energy_stored = energy_full();


    if (energy_stored > MAXENERGYTEST)
        error->warning(FLERR, "Energy of old configuration in fix charge_regulation is > MAXENERGYTEST.");
    for (int i = 0; i < nexchanges; i++) {

        double rand_number = random_equal->uniform();

        if (rand_number < 1.0 / 6) {
            forward_acid();
            nacid_attempts = nacid_attempts + 1;
        }
        else if (rand_number < 2.0 / 6) {
            backward_acid();
            nacid_attempts = nacid_attempts + 1;
        }
        else if (rand_number < 3.0 / 6) {
            forward_base();
            nbase_attempts = nbase_attempts + 1;
        }
        else if ( rand_number < 4.0 / 6) {
            backward_base();
            nbase_attempts = nbase_attempts + 1;
        }
        else if (rand_number < 5.0 / 6) {
            forward_salt();
            nsalt_attempts = nsalt_attempts + 1;
        }
        else {
            backward_salt();
            nsalt_attempts = nsalt_attempts + 1;
        }

        if (add_tags_flag && atom->tag_enable) {
              assign_tags();
        }
    }



    if (triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    comm->exchange();
    atom->nghost = 0;
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
    next_reneighbor = update->ntimestep + nevery;

}

void Fix_charge_regulation::forward_acid() {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
    volume = (xhi - xlo) * (yhi - ylo) * (zhi - zlo) * pow(lb, 3) / 1.66;
    beta = 1.0 / (force->boltz * *target_temperature_tcp);
    double energy_before = energy_stored;
    double factor;
    if (reaction_distance < small) {
        factor = particle_number_A() * volume * pow(10, -pKa)
                 * (1 + pow(10, pH - pI_plus)) / ((1 + particle_number_A_minus()) * (1 + particle_number_S_plus()));
    } else {
        factor = particle_number_A() * (4.0 * PI * pow(reaction_distance, 3) / 3.0) * pow(10, -pKa)
                 * (1 + pow(10, pH - pI_plus)) / ((1 + particle_number_A_minus()) * (1 + particle_number_S_plus()));
    }
    int nlocal_temp = atom->nlocal;
    int *type_temp = atom->type;
    double *charge_temp = atom->q;
    double **x_temp = atom->x;
    double *position = new double[3];
    double *check = new double[7];
    double check_all = 0;
    double flag = 0;
    check = delete_particle(acid_type, 0);
    MPI_Allreduce(&check[0], &check_all, 1, MPI_DOUBLE, MPI_SUM, world);
    if (check_all > 0) {
        position[0] = check[2];
        position[1] = check[3];
        position[2] = check[4];
        if (check[5] > 0) {
            atom->q[int(check[1])] = -1;
        }
        flag = insert_particle_cutoff(cation_type, 1, position);
        double energy_after = energy_full();
        if (energy_after < MAXENERGYTEST &&
            random_equal->uniform() < factor * exp(beta * (energy_before - energy_after))) {
            energy_stored = energy_after;
            nacid_successes += 1;
        } else {
            energy_stored = energy_before;
            atom->natoms--;
            if (flag > 0) {
                atom->nlocal--;
            }
            if (check[5] > 0) {
                atom->q[int(check[1])] = 0;
                atom->mask[int(check[1])] = int(check[6]);
            }
            if (force->kspace) force->kspace->qsum_qsq();
            if (force->pair->tail_flag) force->pair->reinit();
        }
    }
}

void Fix_charge_regulation::backward_acid() {
    int successes_temp = 0;
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
    volume = (xhi - xlo) * (yhi - ylo) * (zhi - zlo) * pow(lb, 3) / 1.66;
    beta = 1.0 / (force->boltz * *target_temperature_tcp);
    double factor;
    if (reaction_distance < small) {
        factor = (1 + particle_number_A()) * volume * pow(10, -pKa)
                 * (1 + pow(10, pH - pI_plus)) / (particle_number_A_minus() * particle_number_S_plus());
    } else {
        factor = (1 + particle_number_A()) * (4.0 * PI * pow(reaction_distance, 3) / 3.0) * pow(10, -pKa)
                 * (1 + pow(10, pH - pI_plus)) / (particle_number_A_minus() * particle_number_S_plus());
    }
    double energy_before = energy_stored;
    int nlocal_temp = atom->nlocal;
    int *type_temp = atom->type;
    double *charge_temp = atom->q;
    double **x_temp = atom->x;
    double *position = new double[3];
    double *check1 = new double[7];
    double *check2 = new double[7];
    double check1_all = 0;
    double check2_all = 0;
    check1 = delete_particle(acid_type, -1);
    MPI_Allreduce(&check1[0], &check1_all, 1, MPI_DOUBLE, MPI_SUM, world);
    successes_temp = 0;
    if (check1_all > 0) {
        position[0] = check1[2];
        position[1] = check1[3];
        position[2] = check1[4];
        check2 = delete_particle_cutoff(cation_type, 1, position);
        MPI_Allreduce(&check2[0], &check2_all, 1, MPI_DOUBLE, MPI_SUM, world);
        if (check2_all > 0) {
            successes_temp = 1;
            position[0] = check1[2];
            position[1] = check1[3];
            position[2] = check1[4];
            if (check1[5] > 0) {
                atom->q[int(check1[1])] = 0;
            }
        } else {
            if (check1[5] > 0)
                atom->q[int(check1[1])] = -1;
        }
    }
    double energy_after = energy_full();
    if (successes_temp > 0 && energy_after < MAXENERGYTEST &&
        random_equal->uniform() < (1.0 / factor) * exp(beta * (energy_before - energy_after))) {
        nacid_successes += 1;
        atom->natoms--;
        energy_stored = energy_after;
        if (check2[5] > 0) {
            atom->avec->copy(atom->nlocal - 1, int(check2[1]), 1);
            atom->nlocal--;
        }
    } else {
        energy_stored = energy_before;
        if (successes_temp > 0) {
            if (check1[5] >= 0) {
                atom->q[int(check1[1])] = -1;
                atom->mask[int(check1[1])] = int(check1[6]);
            }
            if (check2[5] >= 0) {
                atom->q[int(check2[1])] = 1;
                atom->mask[int(check2[1])] = int(check2[6]);
            }
        }
    }
}

void Fix_charge_regulation::forward_base() {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
    volume = (xhi - xlo) * (yhi - ylo) * (zhi - zlo) * pow(lb, 3) / 1.66;
    beta = 1.0 / (force->boltz * *target_temperature_tcp);
    double factor;
    if (reaction_distance < small) {
        factor = particle_number_B() * volume * pow(10, -pKb)
                 * (1 + pow(10, 14 - pH - pI_minus)) /
                 ((1 + particle_number_B_plus()) * (1 + particle_number_S_minus()));
    } else {
        factor = particle_number_B() * (4.0 * PI * pow(reaction_distance, 3) / 3.0) * pow(10, -pKb)
                 * (1 + pow(10, 14 - pH - pI_minus)) /
                 ((1 + particle_number_B_plus()) * (1 + particle_number_S_minus()));
    }
    double energy_before = energy_stored;
    int nlocal_temp = atom->nlocal;
    int *type_temp = atom->type;
    double *charge_temp = atom->q;
    double **x_temp = atom->x;
    double *position = new double[3];
    double *check = new double[7];
    double check_all = 0;
    double flag = 0;
    check = delete_particle(base_type, 0);
    MPI_Allreduce(&check[0], &check_all, 1, MPI_DOUBLE, MPI_SUM, world);
    if (check_all > 0) {
        position[0] = check[2];
        position[1] = check[3];
        position[2] = check[4];
        if (check[5] > 0) {
            atom->q[int(check[1])] = 1;
        }
        flag = insert_particle_cutoff(anion_type, -1, position);
        double energy_after = energy_full();
        if (energy_after < MAXENERGYTEST &&
            random_equal->uniform() < factor * exp(beta * (energy_before - energy_after))) {
            energy_stored = energy_after;
            nbase_successes += 1;
        } else {
            energy_stored = energy_before;
            atom->natoms--;
            if (flag > 0) {
                atom->nlocal--;
            }
            if (check[5] > 0) {
                atom->q[int(check[1])] = 0;
                atom->mask[int(check[1])] = int(check[6]);
            }
            if (force->kspace) force->kspace->qsum_qsq();
            if (force->pair->tail_flag) force->pair->reinit();
        }
    }
}

void Fix_charge_regulation::backward_base() {
    int successes_temp = 0;
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
    volume = (xhi - xlo) * (yhi - ylo) * (zhi - zlo) * pow(lb, 3) / 1.66;
    beta = 1.0 / (force->boltz * *target_temperature_tcp);
    double factor;
    if (reaction_distance < small) {
        factor = (1 + particle_number_B()) * volume * pow(10, -pKb)
                 * (1 + pow(10, 14 - pH - pI_minus)) / (particle_number_B_plus() * particle_number_S_minus());
    } else {
        factor = (1 + particle_number_B()) * (4.0 * PI * pow(reaction_distance, 3) / 3.0) * pow(10, -pKb)
                 * (1 + pow(10, 14 - pH - pI_minus)) / (particle_number_B_plus() * particle_number_S_minus());
    }
    double energy_before = energy_stored;
    int nlocal_temp = atom->nlocal;
    int *type_temp = atom->type;
    double *charge_temp = atom->q;
    double **x_temp = atom->x;
    double *position = new double[3];
    double *check1 = new double[7];
    double *check2 = new double[7];
    double check1_all = 0;
    double check2_all = 0;
    check1 = delete_particle(base_type, +1);
    MPI_Allreduce(&check1[0], &check1_all, 1, MPI_DOUBLE, MPI_SUM, world);
    successes_temp = 0;
    if (check1_all > 0) {
        position[0] = check1[2];
        position[1] = check1[3];
        position[2] = check1[4];
        check2 = delete_particle_cutoff(anion_type, -1, position);
        MPI_Allreduce(&check2[0], &check2_all, 1, MPI_DOUBLE, MPI_SUM, world);
        if (check2_all > 0) {
            successes_temp = 1;
            position[0] = check1[2];
            position[1] = check1[3];
            position[2] = check1[4];
            if (check1[5] > 0) {
                atom->q[int(check1[1])] = 0;
            }
        } else {
            if (check1[5] > 0)
                atom->q[int(check1[1])] = 1;
        }
    }
    double energy_after = energy_full();
    if (successes_temp > 0 && energy_after < MAXENERGYTEST &&
        random_equal->uniform() < (1.0 / factor) * exp(beta * (energy_before - energy_after))) {
        nbase_successes += 1;
        atom->natoms--;
        energy_stored = energy_after;
        if (check2[5] > 0) {
            atom->avec->copy(atom->nlocal - 1, int(check2[1]), 1);
            atom->nlocal--;
        }
    } else {
        energy_stored = energy_before;
        if (successes_temp > 0) {
            if (check1[5] >= 0) {
                atom->q[int(check1[1])] = 1;
                atom->mask[int(check1[1])] = int(check1[6]);
            }
            if (check2[5] >= 0) {
                atom->q[int(check2[1])] = -1;
                atom->mask[int(check2[1])] = int(check2[6]);
            }
        }
    }
}

void Fix_charge_regulation::forward_salt() {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
    volume = (xhi - xlo) * (yhi - ylo) * (zhi - zlo) * pow(lb, 3) / 1.66;
    beta = 1.0 / (force->boltz * *target_temperature_tcp);
    double energy_before = energy_stored;
    double factor = volume * volume * (pow(10, -pH) + pow(10, -pI_plus))
                    * (pow(10, -14 + pH) + pow(10, -pI_minus)) /
                    ((1 + particle_number_S_plus()) * (1 + particle_number_S_minus()));
    int nlocal_temp = atom->nlocal;
    int *type_temp = atom->type;
    double *charge_temp = atom->q;
    double **x_temp = atom->x;
    double flag1 = 0;
    double flag2 = 0;
    flag1 = insert_particle(cation_type, +1);
    flag2 = insert_particle(anion_type, -1);
    double energy_after = energy_full();
    if ( energy_after < MAXENERGYTEST && random_equal->uniform() < factor * exp(beta * (energy_before - energy_after))) {
        energy_stored = energy_after;
        nsalt_successes += 1;
    } else {
        energy_stored = energy_before;
        atom->natoms--;
        if (flag1 > 0) {
            atom->nlocal--;
        }
        atom->natoms--;
        if (flag2 > 0) {
            atom->nlocal--;
        }
        if (force->kspace) force->kspace->qsum_qsq();
        if (force->pair->tail_flag) force->pair->reinit();
    }
}

void Fix_charge_regulation::backward_salt() {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
    volume = (xhi - xlo) * (yhi - ylo) * (zhi - zlo) * pow(lb, 3) / 1.66;
    beta = 1.0 / (force->boltz * *target_temperature_tcp);
    double energy_before = energy_stored;
    double factor = volume * volume * (pow(10, -pH) + pow(10, -pI_plus))
                    * (pow(10, -14 + pH) + pow(10, -pI_minus)) /
                    ((1 + particle_number_S_plus()) * (1 + particle_number_S_minus()));
    int nlocal_temp = atom->nlocal;
    int *type_temp = atom->type;
    double *charge_temp = atom->q;
    double **x_temp = atom->x;
    double *position = new double[3];
    double *check1 = new double[6];
    double *check2 = new double[6];
    double check1_all = 0;
    double check2_all = 0;
    int maxmol, maxmol_all;
    check1 = delete_particle_salt(cation_type, +1);
    MPI_Allreduce(&check1[0], &check1_all, 1, MPI_DOUBLE, MPI_SUM, world);
    if (check1_all > 0) {
        check2 = delete_particle_salt(anion_type, -1);
        MPI_Allreduce(&check2[0], &check2_all, 1, MPI_DOUBLE, MPI_SUM, world);
        if (check2_all > 0) {
            double energy_after = energy_full();
            if (energy_after < MAXENERGYTEST &&
                random_equal->uniform() < (1.0 / factor) * exp(beta * (energy_before - energy_after))) {
                energy_stored = energy_after;
                nsalt_successes += 1;
            } else {
                energy_stored = energy_before;
                if (check1[5] > 0) {
                    position[0] = check1[2];
                    position[1] = check1[3];
                    position[2] = check1[4];
                    atom->avec->create_atom(2, position);
                    int m = atom->nlocal - 1;
                    atom->q[m] = 1;
                    modify->create_attribute(m);
                }
                atom->natoms++;
                if (check2[5] > 0) {
                    position[0] = check2[2];
                    position[1] = check2[3];
                    position[2] = check2[4];
                    atom->avec->create_atom(4, position);
                    int m = atom->nlocal - 1;
                    atom->q[m] = -1;
                    modify->create_attribute(m);
                }
                atom->natoms++;
                if (force->kspace) force->kspace->qsum_qsq();
                if (force->pair->tail_flag) force->pair->reinit();
            }
        } else {
            if (check1[5] > 0) {
                position[0] = check1[2];
                position[1] = check1[3];
                position[2] = check1[4];
                atom->avec->create_atom(2, position);
                int m = atom->nlocal - 1;
                atom->q[m] = 1;
                modify->create_attribute(m);
            }
            atom->natoms++;
        }
    }
    if (force->kspace) force->kspace->qsum_qsq();
    if (force->pair->tail_flag) force->pair->reinit();
}

double Fix_charge_regulation::insert_particle(int ptype, double charge) {
    double coord[3];
    double value;
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    value = 0;
    coord[0] = xlo + random_equal->uniform() * (xhi - xlo);
    coord[1] = ylo + random_equal->uniform() * (yhi - ylo);
    coord[2] = zlo + random_equal->uniform() * (zhi - zlo);
    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
        coord[1] >= sublo[1] && coord[1] < subhi[1] &&
        coord[2] >= sublo[2] && coord[2] < subhi[2]) {
        value = 1;
        atom->avec->create_atom(ptype, coord);
        int m = atom->nlocal - 1;
        atom->mask[m] = groupbitall;
        for (int igroup = 0; igroup < ngrouptypes; igroup++) {
            if (ngcmc_type == grouptypes[igroup])
                atom->mask[m] |= grouptypebits[igroup];
        }
        tfac_insert = 1.0;
        gas_mass = atom->mass[ptype];
        sigma = sqrt(force->boltz * *target_temperature_tcp * tfac_insert / gas_mass / force->mvv2e);
        atom->v[m][0] = random_unequal->gaussian() * sigma;
        atom->v[m][1] = random_unequal->gaussian() * sigma;
        atom->v[m][2] = random_unequal->gaussian() * sigma;
        atom->q[m] = charge;
        modify->create_attribute(m);
    }


    atom->nghost = 0;
    comm->borders();
    atom->natoms++;
    return (value);
}

double Fix_charge_regulation::insert_particle_cutoff(int ptype, double charge, double *target) {
    double coord[3];
    double value;

    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    value = 0;
    if (reaction_distance < small) {
        coord[0] = xlo + random_equal->uniform() * (xhi - xlo);
        coord[1] = ylo + random_equal->uniform() * (yhi - ylo);
        coord[2] = zlo + random_equal->uniform() * (zhi - zlo);
    } else {
        double radius = reaction_distance * random_equal->uniform();
        double theta = random_equal->uniform() * PI;
        double phi = random_equal->uniform() * 2 * PI;
        coord[0] = target[0] + radius * sin(theta) * cos(phi);
        coord[1] = target[1] + radius * sin(theta) * sin(phi);
        coord[2] = target[2] + radius * cos(theta);
        coord[0] = coord[0] - floor(1.0 * (coord[0] - xlo) / (xhi - xlo)) * (xhi - xlo);
        coord[1] = coord[1] - floor(1.0 * (coord[1] - ylo) / (yhi - ylo)) * (yhi - ylo);
        coord[2] = coord[2] - floor(1.0 * (coord[2] - zlo) / (zhi - zlo)) * (zhi - zlo);
        coord[0] = xlo + random_equal->uniform() * (xhi - xlo);
        coord[1] = ylo + random_equal->uniform() * (yhi - ylo);
        coord[2] = zlo + random_equal->uniform() * (zhi - zlo);
    }

    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
        coord[1] >= sublo[1] && coord[1] < subhi[1] &&
        coord[2] >= sublo[2] && coord[2] < subhi[2]) {
        value = 1;
        atom->avec->create_atom(ptype, coord);
        int m = atom->nlocal - 1;
        atom->mask[m] = groupbitall;
        for (int igroup = 0; igroup < ngrouptypes; igroup++) {
            if (ngcmc_type == grouptypes[igroup])
                atom->mask[m] |= grouptypebits[igroup];
        }
        tfac_insert = 1.0;
        gas_mass = atom->mass[ptype];
        sigma = sqrt(force->boltz * *target_temperature_tcp * tfac_insert / gas_mass / force->mvv2e);
        atom->v[m][0] = random_unequal->gaussian() * sigma;
        atom->v[m][1] = random_unequal->gaussian() * sigma;
        atom->v[m][2] = random_unequal->gaussian() * sigma;
        atom->q[m] = charge;
        modify->create_attribute(m);
    }
    atom->nghost = 0;
    comm->borders();
    atom->natoms++;
    return (value);
}

double *Fix_charge_regulation::delete_particle(int ptype, double charge) {
    int nlocal = atom->nlocal;
    ptype_ID = new int[nlocal];
    int count_local, count_global, count_before;
    int tmpmask;
    double *value = new double[7];
    count_local = 0;
    count_global = 0;
    count_before = 0;
    for (int i = 0; i < nlocal; i++) {
        if (atom->type[i] == ptype && fabs(atom->q[i] - charge) < small) {
            ptype_ID[count_local] = i;
            count_local = count_local + 1;
        }
    }
    MPI_Allreduce(&count_local, &count_global, 1, MPI_INT, MPI_SUM, world);
    MPI_Scan(&count_local, &count_before, 1, MPI_INT, MPI_SUM, world);
    count_before -= count_local;

    if (count_global > 0) {
        const int ID_global = floor(random_equal->uniform() * count_global);
        if ((ID_global >= count_before) && (ID_global < (count_before + count_local))) {
            const int ID_local = ID_global - count_before;
            const int m = ptype_ID[ID_local];
            atom->q[m] = 0.0;
            tmpmask = atom->mask[m];
            double *position = new double[3];
            position[0] = atom->x[m][0];
            position[1] = atom->x[m][1];
            position[2] = atom->x[m][2];
            value[0] = 1;//value[0]=1: The particle of given type and charge exists.
            value[1] = double(1.0 * m);//The local ID of the selected particle
            value[2] = position[0];
            value[3] = position[1];
            value[4] = position[2];
            value[5] = 1;//value[5] This labels the processor which contains the selected particle; otherwise it is -1
            value[6] = tmpmask;
        } else {
            value[0] = 0;//value[0]=0: The particle of given type and charge does not exist.
            value[5] = -1;
            value[6] = tmpmask;
        }
    } else {
        value[0] = 0;
        value[5] = -1;
        value[6] = tmpmask;
    }
    return (value);
}

double *Fix_charge_regulation::delete_particle_cutoff(int ptype, double charge, double *target) {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    int nlocal = atom->nlocal;
    ptype_ID = new int[nlocal];
    int count_local, count_global, count_before;
    int tmpmask;
    double *value = new double[7];
    count_local = 0;
    count_global = 0;
    count_before = 0;
    value[0] = 0;
    value[5] = -1;
    double distance_check, dx, dy, dz;
    for (int i = 0; i < nlocal; i++) {
        if (reaction_distance > small) {//reaction_distance<small: No geometry constraint of atom inserting and deleting
            dx = fabs(atom->x[i][0] - target[0]);
            dx -= static_cast<int>(1.0 * dx / (xhi - xlo) + 0.5) * (xhi - xlo);
            dy = fabs(atom->x[i][1] - target[1]);
            dy -= static_cast<int>(1.0 * dy / (yhi - ylo) + 0.5) * (yhi - ylo);
            dz = fabs(atom->x[i][2] - target[2]);
            dz -= static_cast<int>(1.0 * dz / (zhi - zlo) + 0.5) * (zhi - zlo);
            distance_check = sqrt(dx * dx + dy * dy + dz * dz);
        }
        if ((distance_check < reaction_distance || reaction_distance < small) && atom->type[i] == ptype &&
            fabs(atom->q[i] - charge) < small) {
            ptype_ID[count_local] = i;
            count_local = count_local + 1;
        }
    }
    MPI_Allreduce(&count_local, &count_global, 1, MPI_INT, MPI_SUM, world);
    MPI_Scan(&count_local, &count_before, 1, MPI_INT, MPI_SUM, world);
    count_before -= count_local;

    if (count_global > 0) {
        const int ID_global = floor(random_equal->uniform() * count_global);
        if ((ID_global >= count_before) && (ID_global < (count_before + count_local))) {
            const int ID_local = ID_global - count_before;
            const int m = ptype_ID[ID_local];
            atom->q[m] = 0.0;
            tmpmask = atom->mask[m];
            atom->mask[m] = exclusion_group_bit;
            double *position = new double[3];
            position[0] = atom->x[m][0];
            position[1] = atom->x[m][1];
            position[2] = atom->x[m][2];
            value[0] = 1;//value[0]=1: The particle of given type and charge exists.
            value[1] = double(1.0 * m);//The local ID of the selected particle
            value[2] = position[0];
            value[3] = position[1];
            value[4] = position[2];
            value[5] = 1;//value[5] This labels the processor which contains the selected particle; otherwise it is -1
            value[6] = tmpmask;
        } else {
            value[0] = 0;
            value[5] = -1;
            value[6] = tmpmask;
        }
    }
    return (value);
}

double *Fix_charge_regulation::delete_particle_salt(int ptype, double charge) {
    int nlocal = atom->nlocal;
    ptype_ID = new int[nlocal];
    int count_local, count_global, count_before;
    int tmpmask;
    double *value = new double[6];
    count_local = 0;
    count_global = 0;
    count_before = 0;
    for (int i = 0; i < nlocal; i++) {
        if (atom->type[i] == ptype && fabs(atom->q[i] - charge) < small) {
            ptype_ID[count_local] = i;
            count_local = count_local + 1;
        }
    }
    MPI_Allreduce(&count_local, &count_global, 1, MPI_INT, MPI_SUM, world);
    MPI_Scan(&count_local, &count_before, 1, MPI_INT, MPI_SUM, world);
    count_before -= count_local;
    if (count_global > 0) {
        atom->natoms--;
        const int ID_global = floor(random_equal->uniform() * count_global);
        if ((ID_global >= count_before) && (ID_global < (count_before + count_local))) {
            const int ID_local = ID_global - count_before;
            const int m = ptype_ID[ID_local];
            atom->q[m] = 0.0;
            atom->mask[m] = exclusion_group_bit;
            double *position = new double[3];
            position[0] = atom->x[m][0];
            position[1] = atom->x[m][1];
            position[2] = atom->x[m][2];
            atom->avec->copy(atom->nlocal - 1, m, 1);
            atom->nlocal--;
            value[0] = 1;
            value[1] = double(1.0 * m);
            value[2] = position[0];
            value[3] = position[1];
            value[4] = position[2];
            value[5] = 1;
        } else {
            value[0] = 0;
            value[5] = -1;
        }
        if (atom->tag_enable) {
            if (atom->map_style) atom->map_init();
        }
        atom->nghost = 0;
        comm->borders();
    } else {
        value[0] = 0;
        value[5] = -1;
    }
    return (value);
}

int Fix_charge_regulation::particle_number_A() {
    int nlocal = atom->nlocal;
    int *type = atom->type;
    int count = 0;
    double *charge = atom->q;
    for (int i = 0; i < nlocal; i++) {
        if (atom->type[i] == acid_type && fabs(charge[i]) < small)
            count = count + 1;
    }
    int count_sum = 0;
    MPI_Allreduce(&count, &count_sum, 1, MPI_INT, MPI_SUM, world);
    return (count_sum);
}

int Fix_charge_regulation::particle_number_A_minus() {
    int nlocal = atom->nlocal;
    int *type = atom->type;
    int count = 0;
    double *charge = atom->q;
    for (int i = 0; i < nlocal; i++) {
        if (atom->type[i] == acid_type && charge[i] < 0)
            count = count + 1;
    }
    int count_sum = 0;
    MPI_Allreduce(&count, &count_sum, 1, MPI_INT, MPI_SUM, world);
    return (count_sum);
}

int Fix_charge_regulation::particle_number_B() {
    int nlocal = atom->nlocal;
    int *type = atom->type;
    int count = 0;
    double *charge = atom->q;
    for (int i = 0; i < nlocal; i++) {
        if (atom->type[i] == base_type && fabs(charge[i]) < small)
            count = count + 1;
    }
    int count_sum = 0;
    MPI_Allreduce(&count, &count_sum, 1, MPI_INT, MPI_SUM, world);
    return (count_sum);
}

int Fix_charge_regulation::particle_number_B_plus() {
    int nlocal = atom->nlocal;
    int *type = atom->type;
    int count = 0;
    double *charge = atom->q;
    for (int i = 0; i < nlocal; i++) {
        if (atom->type[i] == base_type && charge[i] > 0)
            count = count + 1;
    }
    int count_sum = 0;
    MPI_Allreduce(&count, &count_sum, 1, MPI_INT, MPI_SUM, world);
    return (count_sum);
}

int Fix_charge_regulation::particle_number_S_plus() {
    int nlocal = atom->nlocal;
    int *type = atom->type;
    int count = 0;
    for (int i = 0; i < nlocal; i++) {
        if (atom->type[i] == cation_type)
            count = count + 1;
    }
    int count_sum = 0;
    MPI_Allreduce(&count, &count_sum, 1, MPI_INT, MPI_SUM, world);
    return (count_sum);
}

int Fix_charge_regulation::particle_number_S_minus() {
    int nlocal = atom->nlocal;
    int *type = atom->type;
    int count = 0;
    for (int i = 0; i < nlocal; i++) {
        if (atom->type[i] == anion_type)
            count = count + 1;
    }
    int count_sum = 0;
    MPI_Allreduce(&count, &count_sum, 1, MPI_INT, MPI_SUM, world);
    return (count_sum);
}

double Fix_charge_regulation::energy_full() {
    int imolecule;
    if (triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    comm->exchange();
    atom->nghost = 0;
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal + atom->nghost);
    if (modify->n_pre_neighbor) modify->pre_neighbor();
    neighbor->build(1);
    int eflag = 1;
    int vflag = 0;
    if (overlap_flag) {
        int overlaptestall;
        int overlaptest = 0;
        double delx, dely, delz, rsq;
        double **x = atom->x;
        tagint *molecule = atom->molecule;
        int nall = atom->nlocal + atom->nghost;
        for (int i = 0; i < atom->nlocal; i++) {
            for (int j = i + 1; j < nall; j++) {
                delx = x[i][0] - x[j][0];
                dely = x[i][1] - x[j][1];
                delz = x[i][2] - x[j][2];
                rsq = delx * delx + dely * dely + delz * delz;
                if (rsq < overlap_cutoffsq) {
                    overlaptest = 1;
                    break;
                }
            }
            if (overlaptest) break;
        }
        MPI_Allreduce(&overlaptest, &overlaptestall, 1,
                      MPI_INT, MPI_MAX, world);
        if (overlaptestall) return MAXENERGYSIGNAL;
    }
    size_t nbytes = sizeof(double) * (atom->nlocal + atom->nghost);
    if (nbytes) memset(&atom->f[0][0], 0, 3 * nbytes);

    if (modify->n_pre_force) modify->pre_force(vflag);

    if (force->pair) force->pair->compute(eflag, vflag);

    if (atom->molecular) {
        if (force->bond) force->bond->compute(eflag, vflag);
        if (force->angle) force->angle->compute(eflag, vflag);
        if (force->dihedral) force->dihedral->compute(eflag, vflag);
        if (force->improper) force->improper->compute(eflag, vflag);
    }

    if (force->kspace) force->kspace->compute(eflag, vflag);

    if (modify->n_post_force) modify->post_force(vflag);
    if (modify->n_end_of_step) modify->end_of_step();
    update->eflag_global = update->ntimestep;
    double total_energy = c_pe->compute_scalar();
    return total_energy;
}

double Fix_charge_regulation::compute_vector(int n) {
    double count_temp = 0;
    int nlocal = atom->nlocal;
    if (n == 0) {
        count_temp = nacid_attempts + nbase_attempts + nsalt_attempts;
        return count_temp;
    }
    if (n == 1) {
        count_temp = nacid_successes + nbase_successes + nsalt_successes;
        return count_temp;
    }
    if (n == 2) {
        count_temp = 0;
        for (int i = 0; i < nlocal; i++) {
            if (atom->type[i] == acid_type && fabs(atom->q[i]) < small) {
                count_temp = count_temp + 1;
            }
        }
        double count_temp_sum = 0;
        MPI_Allreduce(&count_temp, &count_temp_sum, 1, MPI_DOUBLE, MPI_SUM, world);
        return count_temp_sum;
    }
    if (n == 3) {
        count_temp = 0;
        for (int i = 0; i < nlocal; i++) {
            if (atom->type[i] == acid_type && fabs(atom->q[i] + 1) < small) {
                count_temp = count_temp + 1;
            }
        }
        double count_temp_sum = 0;
        MPI_Allreduce(&count_temp, &count_temp_sum, 1, MPI_DOUBLE, MPI_SUM, world);
        return count_temp_sum;
    }
    if (n == 4) {
        count_temp = 0;
        for (int i = 0; i < nlocal; i++) {
            if (atom->type[i] == base_type && fabs(atom->q[i]) < small) {
                count_temp = count_temp + 1;
            }
        }
        double count_temp_sum = 0;
        MPI_Allreduce(&count_temp, &count_temp_sum, 1, MPI_DOUBLE, MPI_SUM, world);
        return count_temp_sum;
    }
    if (n == 5) {
        count_temp = 0;
        for (int i = 0; i < nlocal; i++) {
            if (atom->type[i] == base_type && fabs(atom->q[i] - 1) < small) {
                count_temp = count_temp + 1;
            }
        }
        double count_temp_sum = 0;
        MPI_Allreduce(&count_temp, &count_temp_sum, 1, MPI_DOUBLE, MPI_SUM, world);
        return count_temp_sum;
    }
    if (n == 6) {
        count_temp = 0;
        for (int i = 0; i < nlocal; i++) {
            if (atom->type[i] == cation_type) {
                count_temp = count_temp + 1;
            }
        }
        double count_temp_sum = 0;
        MPI_Allreduce(&count_temp, &count_temp_sum, 1, MPI_DOUBLE, MPI_SUM, world);
        return count_temp_sum;
    }
    if (n == 7) {
        count_temp = 0;
        for (int i = 0; i < nlocal; i++) {
            if (atom->type[i] == anion_type) {
                count_temp = count_temp + 1;
            }
        }
        double count_temp_sum = 0;
        MPI_Allreduce(&count_temp, &count_temp_sum, 1, MPI_DOUBLE, MPI_SUM, world);
        return count_temp_sum;
    }
    return 0.0;
}

void Fix_charge_regulation::setThermoTemperaturePointer() {
    int ifix = -1;
    ifix = modify->find_fix(idftemp);
    if (ifix == -1) {
        error->all(FLERR, "Fix charge regulation regulation could not find a temperature fix id provided by tempfixid\n");
    }
    Fix *temperature_fix = modify->fix[ifix];
    int dim;
    target_temperature_tcp = (double *) temperature_fix->extract("t_target", dim);
   // if (fabs(reservoir_temperature - (*target_temperature_tcp)) > small) {
    //    error->all(FLERR, "Temperature of charge regulation MC moves is inconsistent with thermostat temperature");
   // }
}

void Fix_charge_regulation::assign_tags() {
    // Assign tags to ions with zero tags
    if (atom->tag_enable) {
        tagint *tag = atom->tag;
        tagint maxtag_all = 0;
        tagint maxtag = 0;
        for (int i = 0; i < atom->nlocal; i++) maxtag = MAX(maxtag, tag[i]);
        maxtag_all = maxtag;
        MPI_Allreduce(&maxtag, &maxtag_all, 1, MPI_LMP_TAGINT, MPI_MAX, world);
        if (maxtag_all >= MAXTAGINT)
            error->all(FLERR, "New atom IDs exceed maximum allowed ID");

        tagint notag = 0;
        tagint notag_all;
        for (int i = 0; i < atom->nlocal; i++)
            if (tag[i] == 0 && (atom->type[i] == cation_type || atom->type[i] == anion_type))notag++;
        notag_all = notag;
        MPI_Allreduce(&notag, &notag_all, 1, MPI_LMP_TAGINT, MPI_SUM, world);
        if (notag_all >= MAXTAGINT)
            error->all(FLERR, "New atom IDs exceed maximum allowed ID");

        tagint notag_sum = notag;
        MPI_Scan(&notag, &notag_sum, 1, MPI_LMP_TAGINT, MPI_SUM, world);
        // itag = 1st new tag that my untagged atoms should use

        tagint itag = maxtag_all + notag_sum - notag + 1;
        for (int i = 0; i < atom->nlocal; i++) {
            if (tag[i] == 0 && (atom->type[i] == cation_type || atom->type[i] == anion_type)) {
                tag[i] = itag++;
            }
        }
        if (atom->map_style) atom->map_init();
        atom->nghost=0;
        comm->borders();
    }
}


void Fix_charge_regulation::options(int narg, char **arg) {
    if (narg < 0) error->all(FLERR, "Illegal fix charge regulation command");

    // defaults
    add_tags_flag = false;
    nevery = 100;
    nexchanges = 100;
    lb = 0.72;
    reservoir_temperature = 1.0;
    reaction_distance = 0;
    seed = 12345;
    target_temperature_tcp = &reservoir_temperature;

    int iarg = 0;
    while (iarg < narg) {

        if (strcmp(arg[iarg], "lb") == 0) {
            if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge regulation command");
            lb = force->numeric(FLERR, arg[iarg+1]);
            iarg += 2;
        }
        else if (strcmp(arg[iarg], "temp") == 0) {
            if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge regulation command");
            reservoir_temperature = force->numeric(FLERR, arg[iarg+1]);
            iarg += 2;
        }
        else if (strcmp(arg[iarg], "tempfixid") == 0) {
            if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge regulation command");
            int n = strlen(arg[iarg+1]) + 1;
            delete [] idftemp;
            idftemp = new char[n];
            strcpy(idftemp,arg[iarg+1]);
            setThermoTemperaturePointer();
            iarg += 2;
        }
        else if (strcmp(arg[iarg], "rxd") == 0) {
            if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge regulation command");
            reaction_distance = force->numeric(FLERR, arg[iarg+1]);
            iarg += 2;
        }
        else if (strcmp(arg[iarg], "nevery") == 0) {
            if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge regulation command");
            nevery = force->numeric(FLERR, arg[iarg+1]);
            iarg += 2;
        }
        else if (strcmp(arg[iarg], "nexchange") == 0) {
            if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge regulation command");
            nexchanges = force->numeric(FLERR, arg[iarg+1]);
            iarg += 2;
        }
        else if (strcmp(arg[iarg], "seed") == 0) {
            if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge regulation command");
            seed = force->numeric(FLERR, arg[iarg+1]);
            iarg += 2;
        }
        else if (strcmp(arg[iarg], "tag") == 0) {
            if (iarg + 2 > narg) error->all(FLERR, "Illegal fix charge regulation command");
            if (strcmp(arg[iarg + 1], "yes") == 0) {
                add_tags_flag = true;
            } else if (strcmp(arg[iarg + 1], "no") == 0) {
                add_tags_flag = false;
            } else { error->all(FLERR, "Illegal fix charge regulation command"); }
            iarg += 2;
        }

        else { error->all(FLERR, "Illegal fix charge regulation command"); }
    }
}