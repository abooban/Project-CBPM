# Project CBPM
### Introduction to Computer-Based Physical Modeling — Final Project

**First-year coursework · International Physics Studies Program (IPSP) · Leipzig University**

This repository preserves my final project for **Introduction to Computer-Based Physical Modeling**, completed during the first year of the International Physics Studies Program (IPSP) at Leipzig University.

The original submission is a **Jupyter notebook**. This is a *completed course project*, not an ongoing research project or a polished software library.

## Project at a glance

| | |
| --- | --- |
| **Course** | Introduction to Computer-Based Physical Modeling |
| **Program** | International Physics Studies Program (IPSP) |
| **University** | Leipzig University |
| **Level** | First-year undergraduate |
| **Format** | Jupyter Notebook |
| **Status** | Completed coursework |
| **Original repository upload** | July 2024 |

## Explore the project

**[View the original Jupyter notebook](Project.ipynb)**

The notebook is the central project artifact. It is approximately **85 MB** and is tracked with **Git Large File Storage (Git LFS)**, so GitHub may not provide a convenient browser preview.

To examine the notebook in full, download the repository with Git LFS and open it with JupyterLab.

### Open on Ubuntu / Linux / WSL

Install the basic tools (if not already installed):

```bash
sudo apt update
sudo apt install git git-lfs python3 python3-venv
git lfs install
```

Clone the repository and open the notebook:

```bash
git clone https://github.com/abooban/Project-CBPM.git
cd Project-CBPM

python3 -m venv .venv
source .venv/bin/activate
python -m pip install jupyterlab
jupyter lab Project.ipynb
```

This installs **JupyterLab**, the notebook interface. Additional scientific Python packages may be required by individual notebook cells; the original repository does not include a verified dependency list. Consult the imports in the notebook rather than assuming a particular environment.

If the notebook is downloaded as a tiny text file beginning with `version https://git-lfs.github.com/spec/v1`, that is only the LFS pointer. Run `git lfs pull` from inside the cloned repository to retrieve the actual notebook.

## Repository contents

```text
Project-CBPM/
├── Project.ipynb       # Original course-project notebook (Git LFS)
├── .gitattributes      # Git LFS tracking configuration
└── README.md           # Project introduction and viewing instructions
```

## Context

I originally developed this work as part of my first-year physics coursework. I keep it in my GitHub portfolio as a record of my early experience with programming and computational approaches to physics.

**Scope note:** The notebook is preserved as coursework. No claim is made here that its results constitute peer-reviewed research, a production-ready simulation package, or a reproducible benchmark under current software versions.

---

*Original project upload: 2024 · Repository maintained as an academic portfolio archive.*
