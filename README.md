# Associative_Desorption
Continuing from the observation in [Improved Hydrogen Evolution Activity Descriptors from First-principles Electrochemical Kinetics](https://doi.org/10.26434/chemrxiv-2025-p65sw) that associative desorption of $\mathrm{H_2}$ does not follow BEP kinetics, we use [Taskblaster](https://taskblaster.readthedocs.io/en/latest/) to perform high-throughput calculations on the barriers of other diatomic species in vacuum.
Explanation of Scripts:
* materials.py : Workflow for unit cell minimizations
* surfaces.py : Workflow for adsorbate surface calculations
* endstates.py : Workflow for NEB initial,final states and MLNEB calculations
* build_db.py : Run to build the complete ase database of adsorbates across metals and sites
* tasks.py : Masterlist of tasks for each workflow
