# Single Amino Acid Mutation Analysis and Visualization for Base editors

## Overview 
single_AAwindow_mut project consists of two Python scripts that perform comprehensive analysis and visualization of single amino acid mutations. The first script (`Analysis of BE mutation profiles.py`) generates mutation data based on codon mutations, and the second script (`Plot of BE mutation profiles.py`) visualizes these mutation data using heatmaps.

## Prerequisites
- **Python**: This project is developed and tested with Python 3.10.9. It is recommended to use this version or a compatible Python 3.x version to ensure smooth execution.
- **Libraries**: The following Python libraries are required:
  - `pandas`: Used for data manipulation and analysis.
  - `seaborn`: Used for creating statistical graphics, especially heatmaps.
  - `matplotlib`: Used for creating visualizations.
  - `itertools`: A built - in Python library used for generating combinations.

You can install the required libraries using `pip`:
```bash
pip install pandas seaborn matplotlib
```

## Scripts

### 1. `Analysis of BE mutation profiles.py`
#### Functionality
- **Codon Library Generation**: Generates a complete list of all possible codons (64 in total) using the four DNA bases (A, T, C, G).
- **Codon - Amino Acid Mapping**: Defines a codon table that maps each codon to its corresponding amino acid or 'Stop' signal.
- **Mutation Simulation**: Simulates all possible single, double, and triple base mutations for each codon in the codon library.
- **Data Output**: Stores the codon library, codon - amino acid mapping, and mutation results in separate Excel files (`codon_database.xlsx`, `amino_database.xlsx`, `mutation_results.xlsx`, `ori_mut.xlsx`).

#### Usage
Run the script in your Python environment:
```bash
python "Analysis of BE mutation profiles.py"
```
After running the script, you will see the following message indicating that the results have been saved:
```
results have been saved to codon_database.xlsx, amino_database.xlsx, mutation_results.xlsx,ori_mut.xlsx
```

### 2. `Plot of BE mutation profiles.py`
#### Functionality
- **Data Loading**: Reads the mutation data from the `ori_mut.xlsx` file generated by the previous script.
- **Mutation Counting**: Counts the number of mutations between different amino acids and stores the results in a DataFrame.
- **Heatmap Generation**: Creates two heatmaps to visualize the amino acid mutation profiles:
  - A heatmap showing the overall amino acid mutation counts.
  - A heatmap comparing synonymous and non - synonymous mutations at different mutation sites.
- **Plot Saving**: Saves the generated heatmaps as PDF files (`amino_acid_mutation_heatmap.pdf` and `Cluster_amino_acid_mutation_heatmap.pdf`).

#### Usage
Run the script in your Python environment:
```bash
python "Plot of BE mutation profiles.py"
```
After running the script, you will see the following messages indicating that the plots have been saved:
```
plot has been saved in amino_acid_mutation_heatmap.pdf
plot has been saved in Cluster_amino_acid_mutation_heatmap.pdf。
```

## Notes
- Make sure the `Analysis of BE mutation profiles.py` script is run before the `Plot of BE mutation profiles.py` script, as the latter depends on the `ori_mut.xlsx` file generated by the former.
- The code assumes that the Excel files are in the same directory as the Python scripts. If the files are located elsewhere, you need to modify the file paths in the code accordingly.

# BEs_MutMapper: From Base sequence Mutation to Amino Acid sequence Mutation

## Overview
- `BEs_MutMapper.python` is a Python - based tool developed to analyze the base mutation effects of single base editors (ABE, CBE) and dual - function base editors (DFBE) on base sequences of a specified window length. It also calculates the amino acid mutation effects. Users can input any base sequence (length ≤ 51 bases and a multiple of 3) and select the type of base editor. After that, the tool generates a data table of total mutation effects and statistical graphs of mutation scenarios for subsequent analysis.
- A software version of BEs_MutMapper is provided, and users can directly download the `BEs_MutMapper_exe.zip` and unzip it to use the `BEs_MutMapper.exe` tool (window system), without the need to deploy the related environment.
## Features
- **Multiple Editor Types**: Supports ABE (A to G), CBE (C to T), and DFBE (A to G and C to T) base editors.
- **Data Output**: Generates an Excel file (`user_input_mutation_results.xlsx`) containing detailed mutation information.
- **Graphical Visualization**: Creates a PDF file (`amino_acid_mutation_graphs.pdf`) with a distribution histogram and a proportion chart of amino acid change numbers.

## Prerequisites
- Python 3.x
- Required Python libraries:
  - `tkinter`: For creating the graphical user interface.
  - `pandas`: For data manipulation and generating Excel files.
  - `matplotlib`: For creating statistical graphs.
  - `matplotlib.backends.backend_pdf`: For saving graphs as PDF files.

You can install the required libraries using `pip`:
```bash
pip install pandas matplotlib
```

## Installation
1. Clone the repository or download the `BEs_MutMapper.py` file to your local machine.
```bash
git clone <repository_url>
```

## Usage
1. Run the Python script:
```bash
python BEs_MutMapper.py
```
2. A graphical user interface will appear:
   - **Select Base Editor**: Choose from ABE (A to G), CBE (C to T), or DFBE (A to G and C to T) using the radio buttons.
   - **Enter Sequence**: Input a base sequence with a length that is a multiple of 3 and does not exceed 51 bases.
   - **Submit**: Click the "Submit" button to start the analysis.
3. After the analysis is completed, you will see a message indicating that the results have been saved.
   - An Excel file named `user_input_mutation_results.xlsx` will be generated in the current directory, containing information such as the original sequence, mutated sequence, total mutation site number, and amino acid change number.
   - A PDF file named `amino_acid_mutation_graphs.pdf` will also be generated, showing the distribution histogram and proportion chart of amino acid change numbers.

## Code Structure
- **Base and Codon Definitions**: Defines the base list (`A`, `T`, `C`, `G`), codon table, and mapping from three - letter to single - letter amino acids.
- **Mutation Generation**: The `generate_mutations` function generates all possible mutations for a given codon based on the selected editor type.
- **Calculation Function**: The `calculate` function performs the main analysis, including calculating mutation results, saving data to an Excel file, and generating statistical graphs.
- **Graphical User Interface**: Uses `tkinter` to create a simple GUI for user input and interaction.


## Contributing
Contributions are welcome! If you have any suggestions, bug reports, or feature requests, please open an issue or submit a pull request.