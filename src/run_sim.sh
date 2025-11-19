#!/bin/bash

# Ask for testcase name (without extension)
read -p "Enter testcase name (e.g. sample2, ami33): " tc

# Ask for x value
read -p "Enter x value: " x

echo "=== Step 1: Generate LP for ${tc} with x = ${x} ==="
../bin/hw3 "../testcase/${tc}.txt" "../output/${tc}.lp" "${x}"

echo "=== Step 2: Solve with Gurobi ==="
cd ../output || exit 1
gurobi_cl ResultFile="../visualize/${tc}.sol" "${tc}.lp"

echo "=== Step 3: Visualize result ==="
cd ../visualize || exit 1
../bin/visualize "../testcase/${tc}.txt" "${tc}.sol" "${tc}.svg"

echo "Done. Generated ../visualize/${tc}.svg"
