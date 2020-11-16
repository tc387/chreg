# charge-regulation
fix_charge_regulation for LAMMPS, development in progress...

Syntax

fix ID group-ID charge_regulation pH pKa pKb pIp pIm acid_type base_type cation_type anion_type keyword value

ID, group-ID are documented in fix command

charge_regulation = style name of this fix command

pH = pH value of the solution 

pKa = dissociation constant of acid reaction 

pKb = dissociation constant of base reaction

pIp = chemical potential of positively charged species

pIm = chemical potential of negatively charged species

acid_type = atom type of acid groups

base_type = atom type of base groups

cation_type = atom type of free cations

anion_type = atom type of free anions


zero or more keyword/value pairs may be appended

keyword = lb, temp, tempfixid, nevery, nexchange, xrd, seed, tag, group, onlysalt, 

lb value = Bjerrum length (# in the unit of nm)

temp value = temperature 

tempfixid value = fix ID of temperature thermostat

nevery value = invoke this fix every nevery steps

nexchange value = number of charge regulation moves to attempt every nevery steps

xrd value = cutoff distance for acid/base reaction

seed value = random # seed (positive integer)

tag value = yes or no (yes: The code assign unique tags to inserted ions; no: The tag of all inserted ions is "0")

group value = group-ID, inserted ions are assigned to group group-ID. Can be used multiple times to assign inserted ions to multiple groups.

onlysalt values = flag qcation qanion. 
  flag = yes or no (yes: the fix performs only ion insertion/deletion, no: perform acid/base dissociation and ion insertion/deletion
  qcation, qanion = value of cation/anion charge, must be an integer (only provide if flag = yes)

Examples

fix 1 all charge_regulation 7 3 3 4 4 1 2 3 4 

fix chareg all charge_regulation 7 3 3 4 4 1 3 2 4 lb 1.0 nevery 200 nexchange 200 seed 123 rxd 5 tempfixid flang tag yes 

fix chareg all charge_regulation 0 0 0 4 4 0 0 2 4 tempfixid flang tag yes onlysalt yes 2 -1

Description

This fix performs grand canonical Monte Carlo moves of acid and base reactions as well as the exchange of ion pairs with a reservoir…(Add more details)….

Default:

lb = 0.72; temp = 1.0; nevery = 100; nexchange = 100; xrd = 0; seed = 2345; tag = no; onlysalt = no;
