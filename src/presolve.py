from __future__ import print_function # Ensure print works on Python 2
import sys
import os
import gurobipy as gp
from gurobipy import GRB

def save_and_show_presolved_model(input_file_path):
    # Check if file exists
    if not os.path.exists(input_file_path):
        print("Error: File '{}' not found.".format(input_file_path))
        return

    try:
        # 1. Create a Gurobi environment and load the original model
        print("Loading model from {}...".format(input_file_path))
        model = gp.read(input_file_path)

        # 2. Run Presolve
        # The presolve() method returns a NEW model object
        print("Running presolve...")
        model.Params.Presolve = 2
        presolved_model = model.presolve()

        # 3. Construct output filename
        base_name = os.path.splitext(input_file_path)[0]
        output_file_path = "{}_presolved.lp".format(base_name)

        # 4. Write the presolved model to a .lp file
        print("Writing presolved model to {}...".format(output_file_path))
        presolved_model.write(output_file_path)

        print("-" * 30)
        print("Success!")
        print("Original variables: {}".format(model.NumVars))
        print("Original constraints: {}".format(model.NumConstrs))
        print("Presolved variables: {}".format(presolved_model.NumVars))
        print("Presolved constraints: {}".format(presolved_model.NumConstrs))
        print("-" * 30)
        print("--- Content of {} ---".format(output_file_path))
        
        # 5. Read and display the file content
        with open(output_file_path, 'r') as f:
            print(f.read())
            
        print("-" * 30)

    except gp.GurobiError as e:
        print("Gurobi Error: {}".format(e))
    except Exception as e:
        print("Error: {}".format(e))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python presolve.py <path_to_lp_file>")
    else:
        save_and_show_presolved_model(sys.argv[1])
