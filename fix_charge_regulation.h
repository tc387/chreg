/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(charge_regulation,Fix_charge_regulation)

#else

#ifndef LMP_FIX_charge_regulation_H
#define LMP_FIX_charge_regulation_H

#include "fix.h"

namespace LAMMPS_NS {

    class Fix_charge_regulation : public Fix {
    public:
        Fix_charge_regulation(class LAMMPS *, int, char **);
        ~Fix_charge_regulation();
        int setmask();
        void init();
        void pre_exchange();
        void forward_acid();
        void backward_acid();
        void forward_base();
        void backward_base();
        void forward_ions();
        void forward_ions_multival();
        void backward_ions();
        void backward_ions_multival();
        int get_random_particle(int, double, double, double *);
        int insert_particle(int, double, double, double *);
        double energy_full();
        int particle_number(int, double);
        int particle_number_xrd(int, double, double, double *);
        double compute_vector(int n);
        void assign_tags();
        void options(int, char **);
        void setThermoTemperaturePointer();
        double memory_usage();

    private:
        int exclusion_group, exclusion_group_bit;
        int ngcmc_type, nevery, seed;
        int ncycles, nexchanges, nmcmoves;
        double lb, pH, pKa, pKb, pKs, pI_plus, pI_minus;
        int npart_xrd;            // # of particles (ions) within xrd
        int npart_xrd2;            // # of particles (ions) within xrd
        double vlocal_xrd;         // # local volume within xrd
        int exchmode;             // exchange ATOM or MOLECULE
        int movemode;             // move ATOM or MOLECULE
        int regionflag;           // 0 = anywhere in box, 1 = specific region
        int iregion;              // gcmc region
        char *idregion;           // gcmc region id
        bool only_salt_flag;      // true if performing only salt insertion/deletion, no acid/base dissociation.
        bool add_tags_flag;       // true if each inserted atom gets its unique atom tag

        int groupbitall;          // group bitmask for inserted atoms
        int ngroups;              // number of group-ids for inserted atoms
        char **groupstrings;      // list of group-ids for inserted atoms
        int ngrouptypes;          // number of type-based group-ids for inserted atoms
        char **grouptypestrings;  // list of type-based group-ids for inserted atoms
        int *grouptypebits;       // list of type-based group bitmasks
        int *grouptypes;          // list of type-based group types

        double nacid_attempts, nacid_successes, nbase_attempts, nbase_successes, nsalt_attempts, nsalt_successes;
        int nacid_neutral, nacid_charged, nbase_neutral, nbase_charged, ncation, nanion;

        int cr_nmax;              //  max number of local particles
        int max_region_attempts;
        double gas_mass;
        double reservoir_temperature;
        double beta, sigma, volume, volume_rx;
        int salt_charge[2];    // charge of salt ions: [0] - cation, [1] - anion
        int salt_charge_ratio ;
        double xlo, xhi, ylo, yhi, zlo, zhi;
        double region_xlo, region_xhi, region_ylo, region_yhi, region_zlo, region_zhi;
        double region_volume;
        double energy_stored;  // full energy of old/current configuration
        double *sublo, *subhi;
        double **cutsq;
        int *ptype_ID;
        double overlap_cutoffsq; // square distance cutoff for overlap
        int overlap_flag;
        int max_ngas;
        int min_ngas;
        int acid_type, cation_type, base_type, anion_type, reaction_distance_flag;
        double reaction_distance;



        double energy_intra;

        class Pair *pair;

        class RanPark *random_equal;

        class RanPark *random_unequal;

        class Atom *model_atom;

        char *idftemp; // pointer to the temperature fix
        int triclinic;                         // 0 = orthog box, 1 = triclinic

        class Compute *c_pe;

        double *target_temperature_tcp;  // current temperature of the thermostat

    };
}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

Self-explanatory.

*/
