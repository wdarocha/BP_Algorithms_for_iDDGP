# BP_Algorithms_for_iDDGP

This code was developed based on the article:
"A Novel Branch-and-Prune Algorithmic Framework for the 3D Interval Discretizable Distance Geometry Problem: An Approach Based on Torsion Angles of Molecular Structures"
Preprint available at: https://arxiv.org/abs/2508.09143

---

This project implements the **Branch-and-Prune algorithms** for the **Interval Discretizable Distance Geometry Problem (iDDGP)**, with three variants:

- **iBP**  – interval Branch-and-Prune   
- **iABP** – angular Branch-and-Prune  
- **iTBP** – torsional Branch-and-Prune 

The algorithms embed molecular structures from interval distance data and torsional constraints, producing multiple possible conformations under geometric feasibility rules.

---

## 📦 Dependencies

The project is written in **C (C99)** and requires the following libraries:

- **LAPACKE** – C interface to LAPACK  
- **LAPACK** – linear algebra routines  
- **BLAS** – Basic Linear Algebra Subprograms  
- **Math library (libm)** – standard math functions  

### Install on Ubuntu/Debian

```bash
sudo apt update
sudo apt install build-essential cmake git liblapacke-dev liblapack-dev libblas-dev
```

Other systems may require package managers like `brew` (macOS) or `vcpkg` (Windows).

---

## ⚙️ Build Instructions

The repository already provides a `Makefile` that builds the entire project.  
To compile:

```bash
make
```

The executable will be generated in:

```
build/bin/main
```

---

## ▶️ Usage

Run the program with:

```bash
./build/bin/main <input_file> <output_folder>
```

Example:

```bash
./build/bin/main dataset/1A11/input.txt results/
```

### Input File Format

The input file is a **structured text file** with key-value pairs:

```
structure id: 1TOS
structure chain: A
method: iabp
distance constraints file: dataset/I_1TOS_model1_chainA_ddgpHCorder9.dat
cliques and given torsion angles file: dataset/T_1TOS_model1_chainA_ddgpHCorder9.dat
reference structure xyz file: dataset/X_1TOS_model1_chainA_ddgpHCorder9.dat
time limit (days-hours:minutes:seconds): 0-12:00:00
largest admissible distance deviation (in Ångströms): 0.01
angular resolution (in degrees): 5.00
number of solutions (set to 0 for all solutions): 25
sample size: 5
RMSD threshold to consider solutions different (in Ångströms): 3.0
```

- **structure id**: PDB identifier (e.g., `1TOS`)  
- **structure chain**: Chain label (e.g., `A`)  
- **method**: `ibp`, `iabp`, or `itbp`  
- **distance constraints file** (`I_*`): interval distance constraints between atoms  
- **cliques file** (`T_*`): torsional cliques and given torsion angles  
- **reference xyz file** (`X_*`): optional initial conformation in XYZ format  
- **time limit**: maximum execution time  
- **largest admissible distance deviation**: maximum allowed distance error (LDE)
- **angle resolution**: discretization precision  
- **number of solutions**: maximum number of conformations to generate  
- **sample size**: random sampling parameter  
- **RMSD threshold**: cutoff to consider solutions different  

---

## 📂 Input Files

### Distance constraints (`I_*`)
Example:
```
10     7      2      1   1.3423   1.3423    N    C ASN TRP
10     8      2      2   1.0264   1.0264    N    H ASN ASN
10     9      2      2   1.4840   1.4840    N   CA ASN ASN
...
```
Each line specifies atom indices, residue indices, lower/upper distance bounds, atom names, and residue names.

---

### Cliques and torsion angles (`T_*`)
Example:
```
36 34 35 33 0 105.740586 74.259414
37 36 34 35 -1 117.600372 0.000000
38 37 34 36 0 136.329988 43.670012
...
```
Each line encodes a clique and a torsion angle, containing:  
- **Atom indices** – four atoms defining the clique, which in turn determine the torsion angle;  
- **Torsion angle sign** – `+1` or `-1` specify the orientation, while `0` indicates that both orientations are considered;  
- **Torsion angle value** – torsion angle absolute value (in degrees);  
- **Torsion angle deviation** – allowed deviation from the absolute value (in degrees).  

---

### Reference structure (`X_*`)
Example:
```
-2.436 -1.273 5.440
-1.143 0.380 6.203
-2.099 -0.316 5.298
...
```
XYZ coordinates of the reference conformation.

---

## 📤 Output

The program writes results into the specified **output folder**:

- **`<structure_id>_<structure_chain>.pdb`** → reference structure in PDB format (if provided)  
- **`<structure_id>.pdb`** → all generated conformations in PDB format  
- **`results.txt`** → run metrics, including:  
  - CPU time  
  - Number of embedded vertices  
  - Number of solutions found  
  - Maximum deviation errors (MDE/LDE)  
  - Minimum RMSD  

---

## 🔬 Example Run

```bash
./build/bin/main dataset/1TOS/input.txt results/
```

Generates:

```
results/
 ├── 1TOS_0.pdb        # reference structure
 ├── 1TOS.pdb          # embedded conformations
 └── results.txt       # metrics summary
```

Excerpt of `results.txt`:
```
CPU time = 9553.93796900
Last embedded vertex = 52/52
Number of embedded vertices = 729132453
Number of solutions found = 134308150
Number of considered solutions = 8
maximum MDE = 0.00005065
maximum LDE = 0.04222249
minimum RMSD = 0.23850656
```

---

## 📖 Citation

If this code is useful in your research, please cite the preprint below (also available via the **“Cite this repository”** button on GitHub thanks to the included `CITATION.cff`):

**W. A. A. da Rocha, C. Lavor, L. Liberti, L. de Melo Costa, L. D. Secchin, T. E. Malliavin.**  
*A Novel Branch-and-Prune Algorithmic Framework for the 3D Interval Discretizable Distance Geometry Problem: An Approach Based on Torsion Angles of Molecular Structures.*  
**arXiv:2508.09143**, 2025.  
https://arxiv.org/abs/2508.09143

### BibTeX
```bibtex
@misc{darocha2025novelbranchandprunealgorithmicframework,
  title        = {A Novel Branch-and-Prune Algorithmic Framework for the 3D Interval Discretizable Distance Geometry Problem: An Approach Based on Torsion Angles of Molecular Structures},
  author       = {Wagner A. A. da Rocha and Carlile Lavor and Leo Liberti and Leticia de Melo Costa and Leonardo D. Secchin and Therese E. Malliavin},
  year         = {2025},
  eprint       = {2508.09143},
  archivePrefix= {arXiv},
  primaryClass = {q-bio.BM},
  url          = {https://arxiv.org/abs/2508.09143}
}
```

---

## 📜 License

This repository is licensed under the [MIT License](./LICENSE).  
© 2025 Wagner Alan Aparecido da Rocha

---

## 👤 Author

Developed and maintained by [Wagner Alan Aparecido da Rocha](https://github.com/wdarocha).  

