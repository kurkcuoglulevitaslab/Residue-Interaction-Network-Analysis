# **Residue Interaction Network (RIN) analysis**

This repository provides scripts for performing **Residue Interaction Network (RIN) analysis** of biomolecular structures.

The workflow involves two main steps:  

1. **Affinity Matrix Generation** using the Fortran script `affinity_matrix.f`.  
2. **Hub Residue Identification & Graph Spectral Analysis** using the MATLAB script `graph_RIN.m`.




## 📖 **Overview**  

**1. Affinity Matrix Generation**

The `affinity_matrix.f` script computes the **affinity (aᵢⱼ)**, or **local interaction strength**, of residue pairs from an input PDB file. aij=Nij/(sqrt(Ni*Nj)), where Nij is the number of residue-residue atom contacts within 4.5 Å; Ni and Nj are the number of heavy atoms in residues i and j.

**2. Hub Residue Identification & Graph Spectral Analysis**

The `graph_RIN.m` MATLAB script:  

- Represents the protein structure as a **residue interaction network (RIN)**.  
  - **Nodes:** Residue/nucleotide centers of mass.  
  - **Edges:** Weighted by the inverse of local interaction strength, \(1/aij\).  
- Calculates **betweenness centrality (CB)** for each node.  
- Identifies **hub residues** (residues in the top 5% of CB values).  
- Performs **graph spectral analysis** to partition the network into node clusters.



## 🛠 Dependencies

- **Fortran Compiler** (e.g., `gfortran`) → Required for `affinity_matrix.f`.  
- **MATLAB R2019a or later** → Required for `graph_RIN.m`.



## 🚀 **Usage**  

### Step 1: Affinity Matrix Generation (`affinity_matrix.f`)  

**Input:**  
- PDB file containing atomic coordinates (e.g., `7aku-apo.pdb`).

**Output:**  
- `coordmin1.pdb`: Residue centroid coordinates.  
- `connectivity.out`: Local interaction strengths (aij).  

**Key Parameters:**  
- `rcutt`: Distance cutoff for residue interactions (default: 4.5 Å).  
- `rno`: Number of residues in the system.  
- `nresidue`: Number of atoms in the system.  

**Command-line Usage:**  
```bash
# Compile
gfortran -o affinity_matrix affinity_matrix.f  

# Run
./affinity_matrix
```

(Ensure the input PDB file is in the same directory.)

### Step 2: Hub Residue Identification & Graph Spectral Analysis (`graph_RIN.m`)  

**Input:**  
- `connectivity.out` (from Step 1).
- `coordmin1.pdb` (from Step 1).

**Output:**  
- `betweenness.out`: Betweenness centrality (CB) values. 
- `top_5_percent.out`: Hub residues (top 5% of CB values).
- `10modes.txt`: First 10 non-zero eigenvectors of the Laplacian matrix.
- `10modes_normalized.txt`: Normalized eigenvectors (-1 or +1).

  **Figures**:
- `Figure 1`: Residue interaction network showing nodes colored by CB. The plot can be rotated in 3D.
- `Figure 2`: CB of the residues.
- `Figure 3`: CB distributions.
- `Figure 4`: Spectral analysis shown on the network (Fiedler vector).

**MATLAB Usage:**

- Open MATLAB R2019a or later

- Ensure that `connectivity.out` and `coordmin1.pdb` are present in the same directory

  run('graph_RIN.m')


## 📂 **Example Input & Output**

The example_data/ directory provides sample files:

 - **Input**
   * 7aku-apo.pdb (example PDB structure).
   * affinity_matrix.f (Fortran script).
   * graph_RIN.m (MATLAB script).
 - **Output**
   * coordmin1.pdb (centroid coordinates).
   * connectivity.out (affinity values).
   * betweenness.out (betweenness centrality values).
   * top_5_percent.out (hub residues).
   * 10modes.txt & 10modes_normalized.txt (spectral vectors).

  
## 📜 **License**

This project is licensed under the MIT License.
See the LICENSE file for details.

## 📌 **Citation**

If you use this work in your research, please cite this study and the associated publication:

Kurkcuoglu O. (2018) Exploring allosteric communication in multiple states of the bacterial ribosome using residue network analysis, Turkish Journal of Biology, 42(5): 392–404.
DOI: 10.3906/biy-1802-77
  
## ✉️ **Contact**

For questions or issues, please contact:

Dr. Ozge Kurkcuoglu — 📧 olevitas@itu.edu.tr
