# charge-regulation
fix_charge_regulation for LAMMPS, development in progress...

Syntax

fix ID group-ID charge_regulation cation_type anion_type keyword value(s)

ID, group-ID are documented in fix command

charge_regulation = style name of this fix command

cation_type = atom type of free cations

anion_type = atom type of free anions


zero or more keyword/value pairs may be appended

keyword = pH, pKa, pKb, pIp, pIm, acid_type, base_type, lb, temp, tempfixid, nevery, nexchange, xrd, seed, tag, group, onlysalt, 

pH value = pH of the solution (chemical potential of dissociated cations)

pKa value = dissociation constant of acid reaction 

pKb value = dissociation constant of base reaction

pIp value = chemical potential of salt cations

pIm value = chemical potential of salt anions

acid_type = atom type of acid groups

base_type = atom type of base groups

lb value = Bjerrum length (# in the unit of nm)

temp value = temperature 

tempfixid value = fix ID of temperature thermostat

nevery value = invoke this fix every nevery steps

nexchange value = number of charge regulation moves to attempt every nevery steps

xrd value = cutoff distance for acid/base reaction

seed value = random # seed (positive integer)

tag value = yes or no (yes: The code assign unique tags to inserted ions; no: The tag of all inserted ions is "0")

group value = group-ID, inserted ions are assigned to group group-ID. Can be used multiple times to assign inserted ions to multiple groups.

onlysalt values = flag charge_cation charge_anion. 
  flag = yes or no (yes: the fix performs only ion insertion/deletion, no: perform acid/base dissociation and ion insertion/deletion
  charge_cation, charge_anion = value of cation/anion charge, must be an integer (only specify if flag = yes)

Examples

fix 1 all charge_regulation 1 1 pIp 3 pIm 3 

fix chareg all charge_regulation 1 2 acid_type 3 base_type 4 pKa 5 pKb 7 lb 1.0 nevery 200 nexchange 200 seed 123 rxd 5 tempfixid flang tag yes 

fix chareg all charge_regulation 1 2 pIp 3 pIm 3 tempfixid flang tag yes onlysalt yes 2 -1

Description

This fix performs grand canonical Monte Carlo moves of acid and base reactions as well as the exchange of ion pairs with a reservoir…(Add more details)….

Default:

pH = 7; pKa = 100; pKb = 100; pIp = 100; pIm = 100; acid_type = -1; base_type = -1; lb = 0.72; temp = 1.0; nevery = 100; nexchange = 100; xrd = 0; seed = 2345; tag = no; onlysalt = no;
