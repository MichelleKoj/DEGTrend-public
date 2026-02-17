# DEGTrend — User Guide

**DEGTrend** (Differentially Expressed Genes Trending Analysis) is a program that helps you analyze gene expression data. You do **not** need to write code: you run the program, follow the menus, and provide your data files when asked.

**Author:** Michelle Kojekine  
Please feel free to contact me.
Email: michelle.kojekine@mail.huji.ac.il.

**Language:** R (version 4.4 or higher required)

---

## Table of Contents

1. [What You Need Before Starting](#what-you-need-before-starting)
2. [Example Data](#example-data)
3. [Required R Packages (Optional Pre-Install)](#required-r-packages-optional-pre-install)
4. [How to Run DEGTrend](#how-to-run-degtrend)
5. [First Time: Package Installation](#first-time-package-installation)
6. [Step-by-Step: Using the Program](#step-by-step-using-the-program)
7. [Main Menu Options Explained](#main-menu-options-explained)
8. [Input Files: What You Need](#input-files-what-you-need)
9. [Where Your Results Are Saved](#where-your-results-are-saved)
10. [Useful Commands While Running](#useful-commands-while-running)
11. [Troubleshooting](#troubleshooting)

---

## What You Need Before Starting

- **R** installed on your computer (version **4.4 or higher**).
  - If you don’t have R: download it from [https://cran.r-project.org/](https://cran.r-project.org/) (or [https://cloud.r-project.org](https://cloud.r-project.org)).
- **Your data files** ready (DEG tables, count matrix, sample info, etc.—see [Input Files](#input-files-what-you-need)).
- **A folder** where you want all results to be saved (the program will ask you for this).

You do **not** need to install R packages yourself; DEGTrend will install them the first time you run it (see [First Time: Package Installation](#first-time-package-installation)). If you prefer to install everything **before** running the script, use the section below.

---

## Example Data

Example data for DEGTrend is available in the following Google Drive folder:

**[Example data (Google Drive)](https://drive.google.com/drive/folders/1HgtBbxDIWFP47Bp6M3hT0WvLkLlPJlVE?usp=sharing)**

This folder is open for **HUJI (Hebrew University of Jerusalem) students**.

---

## Required R Packages (Optional Pre-Install)

If you want to install the required packages **before** running DEGTrend (e.g. to avoid waiting during the first run), open R or RStudio and run the commands below. You only need to do this once per R installation.

### Step 1 — Install BiocManager (needed for Bioconductor packages)

In the R console, run:

```r
install.packages("BiocManager")
```

### Step 2 — Install CRAN packages

These are the packages that come from the main R repository (CRAN):

| Package        | Package       | Package     | Package   | Package    |
|----------------|---------------|-------------|-----------|------------|
| shiny          | ggVennDiagram | plotly      | DT        | miniUI     |
| sf             | ggupset       | readxl      | ggplot2   | RColorBrewer |
| circlize       | viridis       | reshape2    | openxlsx  | VennDiagram |
| dplyr          | tidyr         | stringr     | matrixStats |           |

Run this in R (one line):

```r
install.packages(c("shiny", "ggVennDiagram", "plotly", "DT", "miniUI", "sf", "ggupset", "readxl", "ggplot2", "RColorBrewer", "circlize", "viridis", "reshape2", "openxlsx", "VennDiagram", "dplyr", "tidyr", "stringr", "matrixStats"), dependencies = TRUE)
```

### Step 3 — Install Bioconductor packages

These are the packages that come from Bioconductor:

| Package             | Package    | Package   | Package     | Package            |
|---------------------|------------|-----------|-------------|--------------------|
| DESeq2              | RUVSeq     | edgeR     | rtracklayer | ComplexHeatmap     |
| Biobase             | SummarizedExperiment |           |             |                    |

Run this in R (one line):

```r
BiocManager::install(c("DESeq2", "RUVSeq", "edgeR", "rtracklayer", "ComplexHeatmap", "Biobase", "SummarizedExperiment"), dependencies = TRUE, update = TRUE, ask = FALSE)
```

### Note

- **grid** is used by DEGTrend but is part of base R; you do **not** need to install it.
- If you are behind a firewall or have no internet when running DEGTrend, pre-installing with internet access ensures the script can run offline later.
- If any package fails to install, see [Troubleshooting](#troubleshooting).

---

## How to Run DEGTrend

1. **Option A — Using RStudio (recommended if you use RStudio)**  
   - Open RStudio.  
   - Go to **File → Open File** and open `DEGTrend.R`.  
   - Click **Source** (or press Ctrl+Shift+S) to run the whole script.  
   - The program will start in the **Console** (bottom of RStudio). Type your answers and press **Enter** when the program asks for input.

2. **Option B — Using R alone**  
   - Open **R** (RGui or R in a terminal).  
   - In the R console, run:
     ```r
     source("C:/path/to/your/DEGTrend.R")
     ```
     Replace `C:/path/to/your/` with the real folder where `DEGTrend.R` is saved.  
   - Example if the file is on your Desktop in a folder named "DEGTrend share":
     ```r
     source("C:/Users/YourName/Desktop/DEGTrend share/DEGTrend.R")
     ```
   - Then answer the questions in the console.

---

## First Time: Package Installation

The first time you run DEGTrend, it will automatically install the R packages it needs. This can take several minutes and may show a lot of text in the console—that is normal.

- You may see questions like “Do you want to install from sources?” You can type **n** (no) and press Enter if you prefer not to compile from source; the program will try to use pre-built packages when possible.
- If everything installs correctly, the program will then ask you for an **output directory** and show the **Main Menu**. If an installation fails, see [Troubleshooting](#troubleshooting).

---

## Step-by-Step: Using the Program

After the packages load, you will see:

### Step 1 — Choose an output folder

The program will ask:

**“Please provide a directory to save all exported files (or type 'quit' to exit):”**

- Type the **full path** to the folder where you want all results saved.  
- Examples (Windows):
  - `C:/Users/YourName/Documents/DEGTrend_Results`
  - `D:/MyProject/Output`
- You can type **quit** (or **exit**) to close the program.

The program will create this folder if it does not exist. All tables, plots, and exports will go inside it.

### Step 2 — Main menu

You will then see the **Main Menu**:

```
Main Menu
[1] Load DEGs tables and analyze
[2] Process count files and perform differential expression analysis
[3] Load genomic data tables and analyze
[4] Quit
Please select an option by entering the corresponding number (or type 'back' to return to the previous step):
```

- Type **1**, **2**, **3**, or **4** and press **Enter**.  
- Typing **back** takes you back to Step 1 (choosing the output folder).

What each option does is described in [Main Menu Options Explained](#main-menu-options-explained).

---

## Main Menu Options Explained

### Option 1 — Load DEGs tables and analyze

Use this when you **already have**:

- A **DEGs file** (differentially expressed genes from another tool or pipeline).  
- A **counts file** (gene counts per sample).  
- A **coldata file** (sample information: which sample belongs to which group/condition).

The program will then ask you for these three files one by one. After that you can run **Automatic** analysis (heatmaps, Venn diagrams, trending plots) or **Manual** analysis (more control over each step).

### Option 2 — Process count files and perform differential expression analysis

Use this when you want DEGTrend to **do the differential expression analysis** for you.

You will choose:

- **[1]** Process count files and perform differential expression analysis  
  - You can provide either a **count matrix** file or a **folder of count files**, plus a **sample info table** if you have one.  
- **[2]** Perform differential expression analysis without processing count files  
  - You provide a **count matrix** and a **sample info table** (no raw count files).

The program will run the analysis (e.g. DESeq2-style) and then you can export results and generate plots.

### Option 3 — Load genomic data tables and analyze

Use this when you have **genomic** (e.g. DNA-level) differential results in a table.

- You provide **one genomic data file** that contains columns such as **Gene**, **comparison**, **log2FoldChange**, **padj**.  
- You can then search for genes of interest, generate heatmaps, and create Venn diagrams.

### Option 4 — Quit

Type **4** to exit the program.

---

## Input Files: What You Need

### Supported file formats

For most inputs, DEGTrend accepts:

- **CSV** (`.csv`)  
- **Text / tab-separated** (`.txt`, `.tsv`)  
- **Excel** (`.xlsx`)

### Option 1 — DEGs file

- **Required columns:**  
  `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`, `comparison`, `Gene`
- The **Gene** column should contain gene names/IDs.  
- **comparison** should describe the comparison (e.g. “ConditionA vs ConditionB”).

### Option 1 — Counts data file

- **Rows:** genes (first column = gene ID/name, or row names).  
- **Columns:** sample IDs (one column per sample).  
- Values are **counts** (raw or normalized, as used by your DEG pipeline).  
- **Sample IDs** in the column headers must match the **row names** (or first column) of the **coldata** file.

### Option 1 — Coldata (sample info) file

- **Rows:** one row per sample; row names (or first column) = **sample ID**.  
- **Columns:** any information about samples (e.g. condition, batch, treatment).  
- Sample IDs must match the **column names** in the counts file.  
- You will later choose which column to use for **grouping** (e.g. “condition”).

### Option 2 — Count matrix

- Same idea as the counts file above: genes as rows, samples as columns, counts as values.  
- First column or row names = gene ID.

### Option 2 — Sample info table

- Same as coldata: sample IDs (row names or first column) and columns for group/condition, etc.

### Option 3 — Genomic data file

- **Required columns:**  
  `Gene`, `comparison`, `log2FoldChange`, `padj`  
- Optional: other columns (e.g. baseMean, stat) if you use them elsewhere.

### Typing file paths (Windows)

- You can use **forward slashes**: `C:/Users/YourName/Documents/my_counts.csv`  
- Or **backslashes**: `C:\Users\YourName\Documents\my_counts.csv`  
- Paths can be in quotes if needed: `"C:/My Folder/counts.csv"`

---

## Where Your Results Are Saved

All results are saved inside the **output folder** you chose at the start. Typical subfolders:

| Folder / file       | Contents |
|---------------------|----------|
| **results/**        | DEG tables, normalized counts, “mega” results (e.g. Excel), coldata exports. |
| **PCA_RLE/**        | PCA and RLE plots (e.g. raw, normalised, RUV, final). |
| **heatmaps/**       | Heatmaps (e.g. by scaling: log2, z-score; by type: combined, up/down). |
| **venn_upset/**      | Venn diagrams and UpSet plots (PNG, SVG) and count/gene tables. |
| **trending/**       | Trending pattern plots and data (e.g. CSV). |
| **genes_of_interest/** | Plots and tables for genes you looked up. |
| **genomic/**        | Genomic analysis outputs: genes of interest, heatmaps. |

Exact names may vary slightly; the program creates the folders it needs. You can open the output folder in Windows Explorer and browse by folder name.

---

## Useful Commands While Running

- **back** — Go back to the previous step or menu (e.g. from a file path question to the main menu).  
- **quit** or **exit** — Exit the program.  
- When asked for a **number** (e.g. menu option), type only the number (e.g. `1`) then Enter.  
- When asked for a **path**, type the full path to the file or folder, then Enter.

---

## Troubleshooting

### “The DEGs file is missing the following required columns…”

- Your table must contain the required column names exactly (see [Option 1 — DEGs file](#option-1--degs-file)).  
- Check spelling and spaces (e.g. `log2FoldChange` not `log2 Fold Change`).  
- Re-export from your analysis tool with these column names.

### “Samples are present in coldata but missing in counts_data”

- The **sample IDs** in the counts file (column names) must match the **sample IDs** in the coldata file (row names or first column).  
- Check for extra spaces, different spelling, or different order. Fix one of the files so the IDs match exactly.

### “Unsupported file format”

- Use **CSV**, **TXT/TSV**, or **XLSX** only.  
- If your file is Excel with another extension, save it as `.xlsx` or export as CSV.

### “Directory exists but is not writable” / “Parent directory is not writable”

- Choose a folder where you have permission to create and write files (e.g. Documents, Desktop, or a project folder).  
- Avoid system or program folders (e.g. “C:\Program Files”).

### Package installation fails (first run)

- Make sure you have **internet** access so R can download packages.  
- If you are behind a **firewall or proxy**, R may need to be configured to use it (your IT or R documentation).  
- Ensure **R version is 4.4 or higher** (in R or RStudio, type `R.version.string` and press Enter).  
- If one package fails, read the error message: it sometimes suggests installing a system dependency (e.g. Rtools on Windows) or trying again later.

### Invalid choice / invalid number

- When the menu says “enter 1, 2, 3, or 4”, type only that number and press Enter.  
- Don’t type the text in brackets (e.g. type `1` not “Load DEGs tables”).

### Program seems stuck

- It may be computing (e.g. differential expression or large heatmaps). Wait a few minutes.  
- If you want to stop, press **Ctrl+C** once (or use the Stop button in RStudio). You can then run the script again and choose **quit** from the menu for a clean exit next time.

---

## Summary

1. Install **R 4.4+** and put your data files in a known place.  
2. Run **DEGTrend.R** (e.g. via RStudio “Source” or `source(".../DEGTrend.R")` in R).  
3. First run: wait for automatic package installation.  
4. Enter the **output folder** path when asked.  
5. Choose **Main Menu** option **1**, **2**, or **3** and provide the requested file paths when asked.  
6. Follow the on-screen prompts; use **back** to go back and **quit** or **exit** to close.  
7. Find your results in the output folder and its subfolders (results, heatmaps, venn_upset, trending, etc.).

If you need to rerun an analysis, start DEGTrend again and point it to the same or a new output folder as you prefer.

