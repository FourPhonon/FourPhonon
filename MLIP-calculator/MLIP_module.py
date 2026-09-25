import os
import glob
import argparse
import torch
from ase.io import read

def write_vasprun_xml(filepath, atoms, forces, energy):
    """
    Writes a VASP-compatible XML file containing the atominfo, 
    energy, and forces blocks required by Phonopy.
    """
    num_atoms = len(atoms)
    
    xml_lines = [
        '<?xml version="1.0" encoding="ISO-8859-1"?>',
        '<modeling>',
        ' <generator>',
        '  <i name="program" type="string">vasp</i>',
        '  <i name="version" type="string">5.4.4</i>',
        ' </generator>',
        ' <atominfo>',
        f'  <atoms>{num_atoms:5d}</atoms>',
        ' </atominfo>',
        ' <calculation>',
        '  <energy>',
        f'   <i name="e_fr_energy">{energy:18.8f}</i>',
        f'   <i name="e_wo_entrp">{energy:18.8f}</i>',
        f'   <i name="e_0_energy">{energy:18.8f}</i>',
        '  </energy>',
        '  <varray name="forces" >'
    ]
    
    for f in forces:
        xml_lines.append(f'   <v>{f[0]:15.8f} {f[1]:15.8f} {f[2]:15.8f} </v>')
        
    xml_lines.extend([
        '  </varray>',
        ' </calculation>',
        '</modeling>'
    ])
    
    with open(filepath, 'w') as f:
        f.write('\n'.join(xml_lines) + '\n')

def get_calculator(model_name, device, model_path):
    """
    Factory function to initialize and return the selected ASE calculator.
    Defaults to the current directory if no path is specified.
    """
    model_name = model_name.lower()
    
    if model_name == "mattersim":
        from mattersim.forcefield import MatterSimCalculator
        path = model_path if model_path else "mattersim-v1.0.0-1M.pth"
        print(f"Loading MatterSim model from: {path}")
        return MatterSimCalculator(load_path=path, device=device)
        
    elif model_name == "mace":
        from mace.calculators import mace_mp
        if model_path:
            print(f"Loading MACE model from: {model_path}")
            return mace_mp(model=model_path, device=device, default_dtype="float32")
        else:
            print("Loading default MACE-MP medium model.")
            return mace_mp(model="medium", device=device, default_dtype="float32")
        
    elif model_name == "chgnet":
        from chgnet.model.dynamics import CHGNetCalculator
        print("Loading default CHGNet model.")
        return CHGNetCalculator(use_device=device)
        
    elif model_name == "m3gnet":
        import matgl
        from matgl.ext.ase import M3GNetCalculator
        if model_path:
            print(f"Loading M3GNet model from: {model_path}")
            pot = matgl.load_model(model_path)
        else:
            print("Loading default M3GNet-MP-2021.2.8-PES model.")
            pot = matgl.load_model("M3GNet-MP-2021.2.8-PES")
        return M3GNetCalculator(potential=pot, state_attr=None, stress_weight=1.0)
        
    elif model_name == "gap":
        from quippy.potential import Potential
        path = model_path if model_path else "gap_model.xml"
        print(f"Loading GAP model from: {path}")
        return Potential(f"IP GAP param_filename={path}")
        
    elif model_name == "deepmd":
        from deepmd.calculator import DP
        path = model_path if model_path else "graph.pb"
        print(f"Loading DeepMD model from: {path}")
        return DP(model=path)
        
    elif model_name == "mtp":
        from mlip_ase import MLIP
        path = model_path if model_path else "pot.mtp"
        print(f"Loading MTP model from: {path}")
        return MLIP(path)
        
    else:
        raise ValueError(f"Model '{model_name}' is not supported.")

def main():
    parser = argparse.ArgumentParser(description="Calculate forces and energies using MLIPs.")
    parser.add_argument('-m', '--model', type=str, required=True, 
                        choices=['mattersim', 'mace', 'chgnet', 'm3gnet', 'gap', 'deepmd', 'mtp'],
                        help="Specify the MLIP model to use.")
    parser.add_argument('-p', '--path', type=str, default=None,
                        help="Path to the model file. Defaults to current directory if omitted.")
    parser.add_argument('-o', '--order', type=int, choices=[3, 4], default=None,
                        help="Specify 3 for thirdorder (3RD.POSCAR.*) or 4 for fourthorder (4TH.POSCAR.*).")
    args = parser.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Initializing {args.model.upper()} on {device}...")
    
    try:
        calculator = get_calculator(args.model, device, args.path)
    except ImportError as e:
        print(f"Error: Missing required library for {args.model}. Details: {e}")
        return
    except Exception as e:
        print(f"Error initializing {args.model}: {e}")
        return
        
    # Determine the correct file pattern based on the --order argument
    if args.order == 3:
        glob_pattern = "3RD.POSCAR.*"
    elif args.order == 4:
        glob_pattern = "4TH.POSCAR.*"
    else:
        glob_pattern = "POSCAR-*"
        
    poscar_files = sorted(glob.glob(glob_pattern))
    
    if not poscar_files:
        print(f"Error: No files matching '{glob_pattern}' were found in the current directory.")
        return
        
    print(f"Found {len(poscar_files)} structures matching {glob_pattern}. Beginning calculations...")

    for poscar in poscar_files:
        # Extract the sequence number correctly based on the naming convention
        if args.order in [3, 4]:
            sequence_num = poscar.split('.')[-1]
        else:
            sequence_num = poscar.split('-')[-1]
            
        xml_filename = f"vasprun-{sequence_num}.xml"
        
        atoms = read(poscar, format="vasp")
        atoms.calc = calculator
        forces = atoms.get_forces()
        energy = atoms.get_potential_energy()
        
        write_vasprun_xml(xml_filename, atoms, forces, energy)
        print(f"  [+] Processed {poscar} -> {xml_filename}")

    print("\nCalculations complete.")

if __name__ == "__main__":
    main()
