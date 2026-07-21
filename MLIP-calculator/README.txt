========================================================================
MLIP INTERFACE USER GUIDE
========================================================================

This documentation provides a workflow for calculating third-order 
and fourth-order interatomic force constants (IFCs) using the 
MLIP_module.py interface with the Thirdorder and Fourthorder 
packages.
------------------------------------------------------------------------
SUPPORTED MODELS
------------------------------------------------------------------------
The following Machine Learning Interatomic Potentials are supported:
* Foundational Models: mattersim, mace, chgnet, m3gnet
* Traditional Models: gap, deepmd, mtp

------------------------------------------------------------------------
COMMAND LINE ARGUMENTS
------------------------------------------------------------------------
The script MLIP_module.py accepts the following flags:

-m, --model  : (Required) Specify the MLIP model to use. Choices include
               mattersim, mace, chgnet, m3gnet, gap, deepmd, mtp.
-p, --path   : (Optional) Absolute or relative path to the pre-trained 
               model file. If omitted, the script searches the current 
               directory for standard default filenames.
-o, --order  : (Optional) Specify the calculation order. Use 3 to target 
               3RD.POSCAR.*** files, or 4 to target 4TH.POSCAR.*** files. 
               If omitted, the script defaults to standard POSCAR-***.

------------------------------------------------------------------------
STEP-BY-STEP WORKFLOW
------------------------------------------------------------------------

------------------------------------------------------------------------
1. Generate Displaced Structures
------------------------------------------------------------------------
Create the required displaced supercell structures in VASP 
format using the 'sow' command.

Example for third-order IFCs:
thirdorder_vasp.py sow na nb nc cutoff

------------------------------------------------------------------------
2. Calculate Forces and Energies
------------------------------------------------------------------------
Execute the machine learning interatomic potential interface in the same 
directory containing your generated POSCAR files. Specify your desired 
model and the order of the force constants using the command line flags. 
This step automatically evaluates the structures and outputs VASP 
compatible XML files.

Example using MACE for third-order structures:
python MLIP_module.py --model mace --order 3

------------------------------------------------------------------------
3. Reap the Force Constants
------------------------------------------------------------------------
Once all vasprun-*.xml files are generated, extract the force constants 
using the "reap" command. The following command finds all matching XML 
files, sorts them numerically by version, and pipes them directly into 
the harvester.

Example command:
find . -maxdepth 1 -name "vasprun-*.xml" | sort -V | thirdorder_vasp.py reap na nb nc cutoff

========================================================================