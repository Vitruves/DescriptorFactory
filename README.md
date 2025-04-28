# desfact: Descriptor Factory

`desfact` is a high-performance command-line tool for calculating a diverse set of chemical descriptors from molecules provided in a CSV file. Built with C++, it leverages RDKit for cheminformatics operations and optionally uses Intel TBB (Threading Building Blocks) for parallel processing to accelerate calculations.

## Features

*   **CSV Input/Output:** Reads SMILES strings and optional properties from a CSV file and writes calculated descriptors to a new CSV file.
*   **RDKit Integration:** Utilizes the RDKit library for molecule parsing, sanitization, and structural calculations.
*   **Diverse Descriptor Set:** Includes molecular descriptors, as well as novel string-based, sum-based, fractional, vague, electronic, and graph-based eigenvalue descriptors derived from molecular properties and structure.
*   **Parallel Processing:** Supports multi-threaded computation using Intel TBB to significantly speed up descriptor calculation on multi-core processors (optional, falls back to single-threaded if TBB is not found).
*   **Efficient I/O:** Implements custom CSV reading and writing with buffering and optional header handling.
*   **Progress Bar:** Provides a visual progress indicator during computation.
*   **Configurable:** Allows specification of input/output paths, descriptors to calculate, batch size, number of threads, SMILES column name, and delimiter via command-line arguments.
*   **Invalid Molecule Handling:** Gracefully handles unparsable SMILES strings, logs errors, and skips calculation for those rows while retaining original input data if possible.

## Dependencies

*   **C++17 Compiler:** GCC, Clang, or MSVC supporting C++17.
*   **CMake:** Version 3.16 or higher.
*   **RDKit:** A recent version of the RDKit C++ library. RDKit must be installed and discoverable by CMake (e.g., installed system-wide or available via `find_package`). Using a package manager like `conda` or `mamba` is a common way to install RDKit dependencies.
*   **Intel TBB (Optional):** Required for multi-threaded parallel processing. If not found by CMake, the project will build for single-threaded execution.

## Building

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/Vitruves/DescriptorFactory
    cd DescriptorFactory
    ```

2.  **Ensure RDKit and TBB are installed:** Follow RDKit and TBB installation instructions for your system. Conda/mamba is often the easiest way:
    ```bash
    # Using mamba (recommended)
    mamba create -n desfact-env -c conda-forge rdkit tbb cmake cxxopts rapidjson
    mamba activate desfact-env
    ```

3.  **Create build directory and run CMake:**
    ```bash
    mkdir build
    cd build
    cmake ..
    ```
    *   CMake will search for RDKit and TBB. Check the CMake output to confirm they are found (`TBB_FOUND`, RDKit components).
    *   If RDKit is not found, check your `CMAKE_PREFIX_PATH` or how RDKit was installed.
    *   If TBB is not found, `WITH_TBB` will not be defined, and the build will proceed in single-threaded mode.
    *   (Optional) To force building without TBB even if found: `cmake .. -DWITH_TBB=OFF`
    *   (Optional) The `CMakeLists.txt` includes a placeholder for CUDA, but descriptor implementations currently don't use it.

4.  **Compile the project:**

    ```bash
    make
    ```

    This will build the shared library (`libdesfact.so` or `.dylib` or `.dll`) and the executable (`desfact`).

    For maximum performance, consider using Profile-Guided Optimization (PGO):

    ```bash
    # From the repository root
    python3 INSTALL.py
    ```
    PGO can improve execution speed by up to 1000% compared to a standard build by optimizing frequently executed code paths based on actual usage patterns. The script automatically:
    - Builds an instrumented binary
    - Runs benchmarks with sample or provided SMILES data
    - Builds the final optimized binary
    
    Run `python3 INSTALL.py --help` for more options, including custom benchmark data.

5.  **Install (Optional):**

    ```bash
    make install
    ```
    By default, this installs the executable to `CMAKE_INSTALL_PREFIX/bin` and the library to `CMAKE_INSTALL_PREFIX/lib`.

## Usage

The main executable `desfact` is run from the command line:

```bash
./build/desfact [options]

# or, if installed:

desfact [options]
```

### Command Line Options

**Required Options:**

*   `-i, --input <path>`: Path to the input CSV file.
*   `-o, --output <path>`: Path to the output CSV file.
*   `-d, --descriptors <list>`: Comma-separated list of descriptor names to calculate, or `all` to calculate all available descriptors.

**Optional Options:**

*   `-l, --list`: List all available descriptors and their descriptions.
*   `-b, --batch-size <size>`: Number of lines to process per batch (default: 64). Larger batches might use more memory but can be faster.
*   `-t, --threads <num>`: Maximum number of parallel threads to use (default: auto-detected, typically CPU cores - 1). Set to 1 for single-threaded execution.
*   `--verbose`: Enable verbose output (DEBUG and INFO level messages). By default, only WARNING, ERROR, and FATAL messages are shown.
*   `-s, --smiles-column <name>`: Name of the column containing SMILES strings (default: SMILES). Case-sensitive.
*   `--delimiter <char>`: CSV delimiter character (default: ,). Supports single characters.
*   `--no-header`: Specify that the input CSV file has no header row. If this flag is used, the first column (index 0) is assumed to be the SMILES column.

### Examples

**List available descriptors:**

```bash
./build/desfact --list
```

**Calculate Molecular Weight and LogP for all molecules in input.csv:**

```bash
./build/desfact -i input.csv -o output.csv -d MolWt,LogP
```

**Calculate all descriptors for molecules in molecules.csv with custom SMILES column:**

```bash
./build/desfact -i molecules.csv -o results.csv -d all -s Canonical_SMILES -t 8 -b 512 --verbose
```

**Process a file with no header, using tab as delimiter:**

```bash
./build/desfact -i data.txt -o data_results.csv -d NumAtoms --no-header --delimiter $'\t'
```
(Note: Using $'\t' for tab might require a shell like Bash)

## Available Descriptors

Here is a list of descriptors implemented in desfact. You can also get this list using the `--list` option.

### Core Descriptors

*   **MolWt** - Molecular Weight
*   **LogP** - Octanol/Water Partition Coefficient (Crippen method)
*   **NumAtoms** - Total number of atoms

### Sum Descriptors

*   **SumEN** - Sum of atom electronegativities
*   **SumCovalentRadii** - Sum of atomic covalent radii
*   **SumVdWVolume** - Sum of Van der Waals atomic volumes
*   **SumPolarizability** - Sum of atomic polarizabilities
*   **SumIonizationEnergy** - Sum of first ionization energies
*   **SumElectronAffinity** - Sum of electron affinities
*   **SumAtoms** - Total number of atoms
*   **SumHeavyAtoms** - Total number of atoms with Z > 1
*   **SumHeteroatoms** - Total number of non-C, non-H atoms
*   **SumHalogens** - Count of F, Cl, Br, I atoms
*   **SumChargedAtoms** - Count of atoms with formal charge ≠ 0
*   **SumBonds** - Total number of bonds
*   **SumDoubleBonds** - Count of double bonds
*   **SumTripleBonds** - Count of triple bonds
*   **SumAromaticBonds** - Count of aromatic bonds
*   **SumPolarBonds** - Count of bonds with significant EN difference
*   **SumUnpolarBonds** - Count of bonds with small or zero EN difference
*   **SumSp3Bonds** - Count of sp3-hybridized central atom bonds
*   **SumSp2Bonds** - Count of sp2-hybridized central atom bonds
*   **SumRings** - Total number of rings
*   **SumAromaticRings** - Total number of aromatic rings
*   **SumRotatableBonds** - Total number of rotatable bonds
*   **SumBridgeAtoms** - Number of atoms connecting multiple rings
*   **SumRingAtoms** - Number of atoms part of any ring
*   **SumENMW** - Sum of (electronegativity × atomic weight)
*   **SumPolMW** - Sum of (polarizability × atomic weight)
*   **SumENRcov** - Sum of (electronegativity × covalent radius)
*   **SumENPol** - Sum of (electronegativity × polarizability)
*   **SumSp3** - Total count of sp3-hybridized atoms
*   **SumSp2** - Total count of sp2-hybridized atoms
*   **SumHBDonors** - Total number of hydrogen bond donors
*   **SumHBAcceptors** - Total number of hydrogen bond acceptors
*   **SumAmideGroups** - Count of amide functional groups (via SMARTS)
*   **SumCarboxylGroups** - Count of carboxyl functional groups (via SMARTS)
*   **SumElectronegativityRing** - EN sum over atoms in rings
*   **SumElectronegativityBonded** - EN sum for atoms with ≥1 bond
*   **SumAtomicRadiusHeavy** - Covalent radius sum for heavy atoms
*   **SumPolzENRatio** - Sum of polarizability / EN
*   **SumIEENRatio** - Sum of IE / EN
*   **SumEAMWEN** - Sum of EA × MW × EN
*   **SumRadicalCount** - Total number of radical atoms
*   **SumRadicalWeight** - Sum of radical_count × MW
*   **SumSp2CAtoms** - Count of sp2-hybridized carbon atoms
*   **SumSpAtomsEN** - Sum of sp-hybridized atoms × EN
*   **SumFormalChargeAbs** - Sum of absolute formal charges
*   **SumFormalChargeHeavy** - Formal charge sum over heavy atoms
*   **SumFormalChargePositive** - Sum of positive formal charges
*   **SumFormalChargeNegative** - Sum of negative formal charges
*   **SumENMWRatio** - Sum of EN / MW
*   **SumPolzRing** - Polarizability sum for atoms in ring systems
*   **SumIEBonded** - IE sum for atoms with ≥1 bond
*   **SumEAHeavyBonded** - EA sum over heavy atoms with ≥1 bond
*   **SumSp3ENWeighted** - Sum of EN × sp3 atoms
*   **SumHeteroMWEN** - Sum of MW × EN for heteroatoms
*   **SumENSp2Heavy** - EN of heavy sp2 atoms
*   **SumENRadicalAtoms** - EN of radical atoms
*   **SumOxHeteroMW** - Sum of abs(ox state) × MW for heteroatoms (approx)
*   **SumENRingHeavy** - EN of heavy atoms in rings
*   **SumPolENBonded** - Sum of polarizability × EN for bonded atoms

### Fractional Descriptors

*   **FcC** - Fractional contribution of Carbon to molecular weight
*   **FcF** - Fluorine contribution to molecular weight
*   **FcO**, **FcCl**, **FcBr**, **FcI**, **FcS**, **FcN** - Other element contributions to molecular weight
*   **FcHalo** - Combined halogen contribution: F, Cl, Br, I
*   **FcHetero** - Combined heteroatom contribution: all non-C, non-H atoms
*   **FcPolar** - Fraction of atoms with electronegativity > C
*   **FcApolar** - Fraction of atoms with electronegativity ≤ C
*   **FcCSp3** - Proportion of carbon-carbon sp3 bonds
*   **FcCSp2** - Proportion of carbon-carbon sp2 bonds
*   **FcUnpol** - Proportion of bonds between atoms of equal electronegativity
*   **FcPol** - Proportion of polar bonds (atoms with different electronegativity)
*   **FcSumPolAt** - Average electronegativity per atom
*   **FcSumPolMW** - Average electronegativity normalized by molecular weight
*   **FcBondAt** - Total bond valence (sum of bond orders) divided by atom count
*   **FcBondN** - Fraction of non-single bonds on nitrogen atoms
*   **FcBondO** - Fraction of non-single bonds on oxygen atoms
*   **FcBondC** - Fraction of non-single bonds on carbon atoms
*   **FcBondS** - Fraction of non-single bonds on sulfur atoms
*   **FcBondP** - Fraction of non-single bonds on phosphorus atoms
*   **FcHBDonors** - Fraction of atoms that are H-bond donors
*   **FcHBAcceptors** - Fraction of atoms that are H-bond acceptors
*   **FcAromaticAtoms** - Fraction of atoms part of aromatic systems
*   **FcRingAtoms** - Fraction of atoms included in rings
*   **FcBridgeAtoms** - Fraction of atoms acting as ring junctions or bridges
*   **FcChargedAtoms** - Fraction of atoms carrying formal charge
*   **FcENAboveAvg** - Fraction of atoms with electronegativity above molecule average
*   **FcENBelowAvg** - Fraction of atoms with electronegativity below molecule average
*   **FcENHigh** - Fraction of atoms with electronegativity > 3.5
*   **FcSmallR** - Fraction of atoms with covalent radius < 0.77 Å
*   **FcLargeR** - Fraction of atoms with covalent radius > 1.1 Å
*   **FcLowPolz** - Fraction of atoms with atomic polarizability < 3.0 Å³
*   **FcHighPolz** - Fraction of atoms with polarizability > 3.5 Å³
*   **FcHighEA** - Fraction of atoms with electron affinity > 2 eV
*   **FcLowEA** - Fraction with EA < 0.5 eV
*   **FcSmallVdW** - Fraction of atoms with VdW volume < 15 Å³
*   **FcLargeVdW** - Fraction of atoms with VdW volume > 25 Å³
*   **FcENMW** - Sum(atom electronegativity × atomic MW) / total MW
*   **FcENBonded** - Fraction of bonds where both atoms have EN > 3.0
*   **FcVdWMW** - Sum(atom VdW volume × MW) / total MW
*   **FcRcovMW** - Sum(atom covalent radius × MW) / total MW
*   **FcHighOxState** - Fractionc of atoms with common oxidation state ≥ +4 (approx)
*   **FcMetalloid** - Fraction of metalloids
*   **FcHETpol** - Fraction of heteroatoms that are also polar
*   **FcHALpol** - Fraction of halogen atoms with polar bonds
*   **FcHeavyPol** - Fraction of heavy atoms (Z > 10) in polar bonds
*   **FcSp3HeavyAtoms** - sp3 heavy atoms / total heavy atoms
*   **FcSp2ENAboveAvg** - sp2 atoms with EN > avg EN / total sp2 atoms
*   **FcIEOdd** - Atoms where int(IE) is odd / total atoms
*   **FcEvenValenceAtoms** - Even valence electron count / total atoms
*   **FcNonHeteroHighEA** - Non-heteroatoms with EA > 1.0 eV / total atoms
*   **FcHeavyLowIE** - Heavy atoms with IE < 11.0 / total heavy atoms
*   **FcHeteroLowPolz** - Heteroatoms with polarizability < 3.0 / total heteroatoms
*   **FcOxENAboveThreshold** - Atoms with abs(ox state) × EN > 10.0 / total atoms (approx)
*   **FcFormalChargeNonZero** - Non-zero formal charge atoms / total atoms
*   **FcFormalChargePositive** - Atoms with formal charge > 0 / total atoms
*   **FcFormalChargeNegative** - Atoms with formal charge < 0 / total atoms
*   **FcRadiusMWRatioAbove1** - Atoms where radius/MW > 0.07 / total atoms (approx)
*   **FcGroup16Atoms** - Atoms in group 16 / total atoms
*   **FcSpAtomsHighEN** - sp-hybridized atoms with EN > 2.5 / total sp atoms
*   **FcRingOxHigh** - Atoms in ring with approx ox state > 2 / total ring atoms
*   **FcHeavyFormalCharge** - Heavy atoms with non-zero formal charge / heavy atom count
*   **FcSp3PolarizabilityAbove10** - sp3 atoms with polarizability > 3.5 / sp3 atom count (approx)
*   **FcHeavyOxNegative** - Heavy atoms with approx ox state < 0 / total heavy atoms
*   **FcENAboveMoleculeAvg** - Atoms with EN > molecular average / total atoms

### String-Based Descriptors

*   **ascii_sum** - Sum of ASCII values of all characters
*   **ascii_average** - Average of ASCII values
*   **ascii_variance** - Variance of ASCII values
*   **ascii_range** - Range between highest and lowest ASCII values
*   **bit_count** - Count of set bits in byte representation
*   **bit_density** - Density of set bits (bit count / total bits)
*   **bit_entropy** - Shannon entropy of bit representation
*   **bit_transition_rate** - Rate of 0->1 and 1->0 transitions
*   **odd_bit_ratio** - Ratio of bits in odd positions that are set
*   **run_length_encoding_size** - Size after simple run-length encoding
*   **character_repetition_score** - Score based on character repetition patterns
*   **nibble_pattern_count** - Count of unique 4-bit patterns
*   **byte_pattern_count** - Count of unique byte patterns
*   **entropy_per_byte** - Shannon entropy normalized by string length
*   **ascii_histogram_uniformity** - Uniformity of ASCII character distribution
*   **bit_clustering_coefficient** - Measure of how bits cluster together (approx)
*   **longest_character_run** - Length of longest run of same character
*   **special_character_density** - Density of non-alphanumeric characters
*   **case_transition_rate** - Rate of transitions between uppercase and lowercase
*   **alternating_character_score** - Score for alternating character patterns (e.g., ABAB)
*   **character_trigram_count** - Count of unique 3-character sequences
*   **letter_digit_ratio** - Ratio of letters to digits
*   **character_bigram_entropy** - Entropy of character bigrams
*   **byte_palindrome_score** - Score for palindromic byte patterns (direct string symmetry)
*   **bit_palindrome_score** - Score for palindromic bit patterns
*   **byte_mirror_symmetry** - Symmetry score for byte patterns considering bit reversal
*   **bit_rotation_symmetry** - Rotation symmetry of bit patterns (0 or 1)
*   **bit_periodicity** - Highest match rate for bit patterns across different periods
*   **bit_autocorrelation** - Autocorrelation of bit patterns (lag 1)
*   **hamming_weight_distribution** - Variance of per-byte Hamming weights (normalized)
*   **byte_parity_distribution** - Fraction of bytes with even parity
*   **bit_run_entropy** - Entropy of bit run lengths
*   **position_weight** - Weighted average of ASCII values by position (normalized)
*   **front_back_ratio** - Ratio of ASCII sums of the front vs back halves
*   **bracket_asymmetry** - Difference in average position of '(' vs ')' (normalized)
*   **char_dispersion** - Average standard deviation of character positions (normalized)
*   **digit_position_mean** - Mean position of digits (normalized)
*   **symmetry_score** - Overall string symmetry score (palindrome check)
*   **aromatic_ratio** - Ratio of lowercase characters (approximate aromatic count) to total letters
*   **heteroatom_ratio** - Ratio of non-C,H characters based on basic parsing
*   **element_diversity** - Diversity of elements based on basic parsing
*   **saturation_index** - Index based on explicit bond symbols ('-', '=', '#', ':')
*   **pattern_entropy** - Shannon entropy of simple SMILES tokens
*   **branch_complexity** - Complexity based on bracket nesting depth and total branches
*   **ring_position_variance** - Variance in position of ring closing digits (normalized)
*   **bond_position_entropy** - Entropy of distances between explicit bond symbols
*   **element_position_bias** - Average standard deviation of normalized element positions
*   **bond_complexity** - Complexity based on explicit bond types
*   **ring_complexity** - Complexity based on ring closure digits (unique indices * total closures)

### Vague3 Descriptors (Connectivity, Shape, Atom/Bond Types)

*   **AtomCentralityVar** - Variance of atom closeness centrality in molecular graph
*   **PeripheryCoreRatio** - Peripheral atoms (degree ≤1) / core atoms (degree ≥3)
*   **MSTreeDepth** - Depth of the molecular minimum spanning tree
*   **TerminalBranchRatio** - Terminal branches / total branches
*   **AvgAtomPathRedundancy** - Mean number of alternative shortest paths between atom pairs
*   **CarbonClusterDensity** - Clusters of ≥3 contiguous carbon atoms / total carbons
*   **HeteroatomClusterDensity** - Clusters of ≥2 contiguous heteroatoms / total heteroatoms
*   **TerminalHeavyAtomDiversity** - Diversity of terminal atom types
*   **CentralHeteroatomRatio** - Central heteroatoms (degree ≥3) / total heteroatoms
*   **ChainEndElementDiversity** - Unique element types at chain ends
*   **HighENNeighborDensity** - Atoms bonded to atoms with EN > 3.0 / total atoms
*   **ElectronDeficientCRatio** - Carbons bonded directly to electronegative substituents / total carbons
*   **PolarizableAtomDispersion** - Variance of distances between atoms with high polarizability
*   **DipoleMomentProxy** - Fraction of bonds between atoms with EN difference > 1.0
*   **LongChainFragmentDensity** - Number of chain fragments ≥5 atoms / total chain fragments
*   **ShortChainDensity** - Chain fragments ≤3 atoms / total chain fragments
*   **SubstitutionDensityPerRing** - Mean substituents per ring atom
*   **LinearVsBranchedCRatio** - Linear-chain carbons / branched carbons
*   **RingToBranchConnectivityRatio** - Ring atoms connected directly to branches / total ring atoms
*   **IsolatedHBondSiteDensity** - H-bond capable atoms with no neighboring H-bond capable atoms
*   **FunctionalGroupSeparationIndex** - Mean shortest path between different functional group atoms
*   **PeripheralHBondDensity** - Peripheral atoms capable of hydrogen bonding / total peripheral atoms
*   **AromaticCoreRatio** - Aromatic atoms in ring cores / total aromatic atoms
*   **ConjugationGapRatio** - Number of conjugation breaks / total conjugated fragments
*   **NonAromaticRingSubstitution** - Substituents on non-aromatic rings / total ring substituents
*   **AromaticToNonAromaticRingRatio** - Aromatic rings / total rings
*   **ConjugationTerminalDensity** - Terminal atoms involved in conjugation / total conjugated atoms
*   **ChargePairDensity** - Count of pairs of oppositely charged atoms separated by ≤4 bonds
*   **FormalChargeAccessibility** - Charged atoms at molecular periphery / total charged atoms
*   **ChargeSeparationIndex** - Mean graph distance between charged atoms
*   **NeutralAtomRatio** - Atoms with zero formal charge / total atoms
*   **FunctionalGroupChargeBalance** - Charged functional groups / total functional groups
*   **StereocenterDensity** - Chiral centers / total heavy atoms
*   **DoubleBondConfigurationDensity** - Double bonds with distinct substituents / total double bonds
*   **RingStereogenicAtomDensity** - Chiral centers within rings / total ring atoms
*   **TerminalStereocenterRatio** - Chiral terminal atoms / total terminal atoms
*   **AdjacentStereocenterDensity** - Pairs of directly bonded stereocenters / total stereocenters
*   **GroupPeriodicDiversity** - Diversity index across periodic groups
*   **PeriodDiversityIndex** - Diversity index across periodic periods
*   **AtomicMassVariance** - Variance of atomic masses of constituent atoms
*   **RareElementDensity** - Fraction of atoms from uncommon periodic groups (group number ≥15, excluding common heteroatoms)
*   **AlkaliAlkalineEarthRatio** - Alkali and alkaline earth atoms / total atoms
*   **UnsaturationClusters** - Clusters of ≥2 adjacent unsaturated bonds (double or triple) / total unsaturated bonds
*   **RingSpanningBondDensity** - Bonds connecting different rings / total bonds
*   **TerminalUnsaturatedBondRatio** - Terminal unsaturated bonds / total unsaturated bonds
*   **HeavyAtomBondOrderVariance** - Variance in bond orders per heavy atom
*   **HeteroatomBondingDiversity** - Average distinct atom types bonded to each heteroatom

### Vague4 Descriptors (More Diverse Structural/Electronic/Topological)

*   **AtomicConnectivityImbalance** - Average absolute difference in degrees between bonded heavy atoms
*   **RingBridgeRatio** - Fraction of bonds bridging two different SSSR ring systems
*   **SubstitutionPatternComplexity** - Number of unique substitution patterns on aromatic rings (heuristic)
*   **OpenChainSaturationRatio** - Fraction of saturated bonds in open chains
*   **HeteroatomNeighborhoodDiversity** - Avg variance of atomic numbers of heteroatom neighbors per atom
*   **CarbonNeighborhoodUniformity** - Fraction of carbon atoms bonded exclusively to carbons
*   **PolarAtomNeighborhoodRatio** - Fraction of polar atoms adjacent only to other polar atoms
*   **AtomDegreeRange** - Difference between max and min atom degrees
*   **RingJunctionComplexity** - Number of atoms connecting >= 2 distinct SSSR ring systems
*   **NonFusedRingDensity** - Fraction of rings that are non-fused (isolated)
*   **RingChainAttachmentDensity** - Fraction of chains attached to rings
*   **RingTerminalSubstituentRatio** - Fraction of ring substituents that are terminal atoms
*   **RingSaturationBalance** - Ratio of saturated (sp3) to unsaturated ring atoms
*   **PolarizabilityGradient** - Average absolute polarizability difference between directly bonded atoms
*   **ElectronWithdrawingAtomDensity** - Fraction of atoms bonded directly to strong EWGs (NO2, CN, CF3)
*   **ElectronDonatingAtomDensity** - Fraction of atoms bonded directly to strong EDGs (alkoxy, amino)
*   **ElectronegativityGradientDensity** - Fraction of bonds with electronegativity difference > 1.2
*   **PeripheralElectronRichAtomRatio** - Fraction of peripheral heavy atoms with lone electron pairs
*   **FunctionalGroupIsolationIndex** - Fraction of functional groups separated by >= 4 non-functional atoms from others
*   **HydroxylGroupDispersion** - Variance in distances among oxygen atoms of hydroxyl groups
*   **AlkylChainDiversity** - Number of distinct alkyl chain lengths per molecule
*   **SubstituentPositionVariance** - Variance of attachment positions of substituents on aromatic rings
*   **TerminalFunctionalGroupClustering** - Fraction of terminal functional groups clustered within 2 bonds of another
*   **ChainSaturationVariability** - Variance of saturated bond fraction across different chains
*   **TripleBondTerminalRatio** - Fraction of triple bonds involving a terminal atom
*   **AdjacentHybridizationTransition** - Fraction of bonds connecting atoms with different hybridization
*   **CyclicHybridizationHomogeneity** - Average variance of hybridization states within each ring
*   **InternalChainUnsaturationDensity** - Fraction of unsaturated bonds strictly internal within chains
*   **HydrogenBondDonorClustering** - Fraction of H-bond donors clustered within <= 2 bonds of another donor
*   **AcceptorDonorRatioImbalance** - Absolute difference between H-bond donors and acceptors normalized by heavy atoms
*   **PeripheralDonorAcceptorBalance** - Ratio of peripheral H-bond donors to peripheral acceptors
*   **IntraringHBondPotential** - Potential H-bond pairs within rings normalized by total ring atoms
*   **ChainEndHydrogenBondDensity** - Fraction of chain-terminal heavy atoms capable of H-bonding
*   **FormalChargeNeighborhoodVariance** - Average variance of formal charges among neighbors for each atom
*   **OppositeChargeNeighborRatio** - Fraction of bonded charged atom pairs with opposite charges
*   **ChargeGradientAlongChains** - Count of formal charge sign changes along linear chain sequences
*   **PeripheralChargeNeutrality** - Fraction of peripheral heavy atoms with zero charge
*   **InterRingConjugationRatio** - Fraction of conjugated bonds connecting different aromatic rings
*   **AromaticChainDensity** - Fraction of aromatic atoms having at least one non-ring neighbor
*   **ConjugationLengthVariance** - Variance in lengths (number of bonds) of conjugated segments
*   **TerminalAromaticSubstitution** - Fraction of substituents on aromatic rings that are terminal
*   **CyclicVsChainAromaticRatio** - Ratio of aromatic atoms in SSSR rings vs aromatic atoms not in SSSR rings
*   **UniqueElementPairRatio** - Fraction of unique bonded element pairs relative to total bonds
*   **HeavyAtomDegreeDiversity** - Shannon diversity index of heavy atom degrees
*   **AtomPathDiversity** - Number of unique shortest path lengths between distinct heavy element types
*   **InternalAtomComplexity** - Number of internal heavy atoms bonded to atoms of >= 3 different elements
*   **HeteroatomBondOrderVariability** - Average variance in bond orders around each heteroatom

### Electronic Descriptors (VdW Radius & Gasteiger Charge Based)

*   **SumRadiusCharge** - Sum of (VdW Radius * Gasteiger Partial Charge)
*   **AvgRadiusCharge** - Average (VdW Radius * Gasteiger Partial Charge)
*   **StdDevRadiusCharge** - Std Dev of (VdW Radius * Gasteiger Partial Charge)
*   **RangeRadiusCharge** - Range of (VdW Radius * Gasteiger Partial Charge)
*   **SkewRadiusCharge** - Skewness of (VdW Radius * Gasteiger Partial Charge)
*   **KurtosisRadiusCharge** - Kurtosis of (VdW Radius * Gasteiger Partial Charge)
*   **SumRadiusPerAbsCharge** - Sum of VdW Radius / (|Gasteiger Charge| + epsilon)
*   **AvgRadiusPositiveCharge** - Average VdW Radius of positively charged atoms
*   **AvgRadiusNegativeCharge** - Average VdW Radius of negatively charged atoms
*   **RatioSumRadiusCharged** - Ratio: Sum VdW Radius (Pos Charge) / Sum VdW Radius (Neg Charge)
*   **WeightedAvgElectronegativity** - Weighted Average Electronegativity (by VdW Radius)
*   **WeightedStdDevElectronegativity** - Std Dev of Electronegativity (weighted by VdW Radius)
*   **SumRadiusSqAbsCharge** - Sum of VdW Radius^2 * abs(Gasteiger Charge)
*   **WeightedStdDevCharge** - Std Dev of Gasteiger Charge (weighted by VdW Radius)
*   **TotalAbsChargePerRadius** - Total (Abs(Charge) / VdW Radius)
*   **SumInvRadiusAbsCharge** - Sum of (1 / VdW Radius * Abs(Gasteiger Charge))
*   **SumEnRadiusSq** - Sum of Electronegativity * (VdW Radius)^2
*   **RatioAvgEnRadiusRingChain** - Ratio of Avg(EN*Radius) for ring atoms vs chain atoms
*   **AvgRatioChargeRadius** - Average Ratio (Gasteiger Charge / VdW Radius) over all atoms
*   **RatioAvgRChargeToAvgRPlusCharge** - Ratio of Avg (R * |q|) to Avg (R + |q|)
*   **AvgRadiusLonePairAtoms** - Avg VdW Radius for atoms with lone pairs (heuristic)
*   **AvgRadiusPiSystemAtoms** - Avg VdW Radius for atoms in pi systems
*   **StdDevRatioChargeRadius** - Std Dev of Ratios (|Gasteiger Charge| / VdW Radius)
*   **SumRadiusPositiveChargeThresh** - Sum of VdW Radius for atoms with Gasteiger Charge > 0.1
*   **SumRadiusNegativeChargeThresh** - Sum of VdW Radius for atoms with Gasteiger Charge < -0.1
*   **WeightedPathCount** - Weighted path count (R*|q|)
*   **AvgWeightedPathLength** - Average weighted path length (R*|q|)
*   **BalabanLikeIndexRCharge** - Balaban-like index (R+|q|)
*   **WienerLikeIndexRCharge** - Wiener-like index (R*|q|)
*   **EccentricityMaxRadiusAtomChargeWeighted** - Eccentricity of max R atom, weighted by its charge
*   **TopoDistMaxRMaxPosCharge** - Topological distance between max R atom and max positive charge atom
*   **TopoDistMaxRMinNegCharge** - Topological distance between max R atom and min negative charge atom
*   **AvgNeighborWeightRq** - Avg neighbor weight (VdW Radius * abs(Charge)) per heavy atom
*   **TopoAutocorrRqDist2** - Topological Autocorrelation (R*|q|) at distance 2
*   **TopoAutocorrRqDist3** - Topological Autocorrelation (R*|q|) at distance 3
*   **TopoAutocorrRqDist4** - Topological Autocorrelation (R*|q|) at distance 4
*   **AvgTopoDistPiWeightedRq** - Avg Topo Dist between Pi atoms, weighted by sum(R*|q|)
*   **AvgTopoDistLPWeightedRq** - Avg Topo Dist between Lone Pair atoms, weighted by sum(R*|q|)
*   **AvgTopoDistPiLpWeightedRq** - Avg Topo Dist between Pi and Lone Pair atoms, weighted by sum(R*|q|)
*   **PathCountAlternatingRCharge** - Count paths length 5 with alternating R/|q|
*   **RatioLongestWeightedPathRingChain** - Ratio of longest weighted paths (R*|q|) ring/chain
*   **SumWeightRChargeDegreeGt3** - Sum of (R+|q|) for atoms with degree > 3
*   **RatioAvgRWeightedChargeTerminal** - Ratio of Avg(R*|q|) for terminal vs non-terminal heavy atoms
*   **EigenvalueWeightedConnectivity** - Leading eigenvalue of weighted connectivity matrix
*   **WeightedBranchPointComplexity** - Sum of branch point R*|q| * Avg(branch R*|q|)
*   **ShannonEntropyRCharge** - Shannon Entropy of normalized (R*|q|) distribution
*   **AvgRadiusPlusChargePerDegree** - Average (VdW Radius + abs(Charge)) / Degree for heavy atoms
*   **SumDeltaRadiusChargeWeighted** - Sum of (VdW Radius - Covalent Radius) * abs(Charge)
*   **WeightedBuriedAtomCount** - Count of 'buried' atoms weighted by R*|q|
*   **RatioSumRqFormalCharge** - Ratio of Sum(R*|q|) for formally charged vs uncharged heavy atoms
*   **AvgRadiusHBDonorWeightedCharge** - Avg VdW Radius of HBD heavy atoms, weighted by their charge
*   **AvgRadiusHBAcceptorWeightedCharge** - Avg VdW Radius of HBA heavy atoms, weighted by their charge
*   **RatioSumRPolarNonpolarFrag** - Ratio Sum(R) polar frags / Sum(R) nonpolar frags
*   **SumRadiusSp2CWeightedCharge** - Sum of VdW Radius for sp2 carbons, weighted by Gasteiger charge
*   **SumRadiusSp3CWeightedCharge** - Sum of VdW Radius for sp3 carbons, weighted by Gasteiger charge
*   **AvgNeighborChargeRadiusPerDegree** - Avg Ratio: (R_heavy * Avg |q|_neighbors) / Degree
*   **KierShapeIndexVariant3** - Kier Kappa3 index variant using R*|q|

### Eigenvalue-Based Graph Descriptors

#### Adjacency Matrix Descriptors
*   **AdjacencyNonZeroEntries** - Number of non-zero entries in adjacency matrix (edges)
*   **AdjacencyMatrixTrace** - Trace of adjacency matrix
*   **AdjacencyFrobeniusNorm** - Frobenius norm of adjacency matrix
*   **AdjacencySpectralRadius** - Spectral radius (largest eigenvalue) of adjacency matrix
*   **AdjacencySmallestEigenvalue** - Smallest eigenvalue of adjacency matrix
*   **AdjacencyEigenvalueKurtosis** - Kurtosis of eigenvalues of adjacency matrix
*   **AdjacencyPositiveEigenvalues** - Number of positive eigenvalues of adjacency matrix
*   **AdjacencyNegativeEigenvalues** - Number of negative eigenvalues of adjacency matrix
*   **AdjacencyZeroEigenvalues** - Number of zero eigenvalues of adjacency matrix
*   **AdjacencyMatrixRank** - Rank of adjacency matrix
*   **GraphEnergy** - Energy of the graph (sum of absolute eigenvalues)

#### Laplacian Matrix Descriptors
*   **LaplacianSpectralRadius** - Spectral radius of Laplacian matrix
*   **LaplacianAlgebraicConnectivity** - Algebraic connectivity (second-smallest eigenvalue of Laplacian)
*   **LaplacianZeroEigenvalues** - Number of zero eigenvalues in Laplacian (number of components)
*   **LaplacianEnergy** - Laplacian energy (sum of squares of eigenvalues)
*   **LaplacianMatrixTrace** - Trace of Laplacian matrix
*   **LaplacianTotalEffectiveResistance** - Total effective resistance (sum of inverse nonzero Laplacian eigenvalues)
*   **LaplacianKirchhoffIndex** - Kirchhoff index (sum of reciprocal Laplacian eigenvalues)
*   **LaplacianEigenvalueVariance** - Variance of Laplacian eigenvalues
*   **LaplacianEigenvalueSkewness** - Skewness of Laplacian eigenvalues

#### Normalized Laplacian Descriptors
*   **NormalizedLaplacianSpectralRadius** - Spectral radius of normalized Laplacian
*   **NormalizedLaplacianSmallestNonzero** - Smallest nonzero normalized Laplacian eigenvalue
*   **NormalizedLaplacianLargestEigenvalue** - Largest normalized Laplacian eigenvalue
*   **NormalizedLaplacianEnergy** - Normalized Laplacian energy
*   **NormalizedLaplacianTrace** - Normalized Laplacian trace

#### Walk-Based Graph Descriptors
*   **NumberOf2Walks** - Number of 2-walks (Tr(A²))
*   **NumberOf3Walks** - Number of 3-walks (Tr(A³))
*   **NumberOf4Walks** - Number of 4-walks (Tr(A⁴))
*   **MeanClosed3WalksPerNode** - Mean number of closed 3-walks per node
*   **MeanClosed4WalksPerNode** - Mean number of closed 4-walks per node
*   **Walk2Energy** - 2-walk energy (sum of singular values of A²)
*   **Walk3Energy** - 3-walk energy (sum of singular values of A³)
*   **Walk4Energy** - 4-walk energy (sum of singular values of A⁴)
*   **GraphIrregularityWalkCount** - Graph irregularity based on walk counts
*   **TrianglesToPathsRatio** - Ratio Tr(A³)/Tr(A²) (triangles to paths)
*   **NumberTriangles** - Number of triangles in the graph (from trace of A³/6)

#### Singular Value Decomposition Descriptors
*   **MaxSingularValue** - Maximum singular value of adjacency matrix
*   **MinNonZeroSingularValue** - Minimum non-zero singular value of adjacency matrix
*   **ConditionNumber** - Condition number (max/min singular value)
*   **SumSingularValues** - Sum of singular values of adjacency matrix
*   **FrobeniusNormSingularValues** - Frobenius norm from singular values
*   **SingularValueEntropy** - Entropy of singular values after normalization
*   **SingularValueVariance** - Variance of singular values
*   **SingularValueSkewness** - Skewness of singular values
*   **SpectralEffectiveRank** - Spectral effective rank (Shannon entropy of normalized singular values)
*   **NuclearNorm** - Nuclear norm (sum of singular values)

#### Normalized Adjacency Matrix Descriptors
*   **NormalizedAdjacencySpectralRadius** - Normalized adjacency spectral radius
*   **NormalizedEigenvalueGap** - Gap between first two normalized eigenvalues
*   **SumNormalizedEigenvalues** - Sum of normalized eigenvalues
*   **VarianceNormalizedEigenvalues** - Variance of normalized eigenvalues
*   **CountNormalizedEigenvaluesAboveHalf** - Number of normalized eigenvalues > 0.5
*   **NormalizedEnergy** - Normalized energy (sum of absolute normalized eigenvalues)
*   **LargestNormalizedEigenvectorCentrality** - Largest normalized eigenvector centrality
*   **AverageNormalizedEigenvectorCentrality** - Average normalized eigenvector centrality
*   **NormalizedAdjacencyMatrixRank** - Normalized adjacency matrix rank
*   **NormalizedAdjacencyEntropy** - Normalized adjacency entropy

#### Signless Laplacian Matrix Descriptors
*   **SignlessLaplacianSpectralRadius** - Signless Laplacian spectral radius
*   **SignlessLaplacianSmallestEigenvalue** - Smallest eigenvalue of signless Laplacian
*   **SignlessLaplacianEnergy** - Energy of signless Laplacian (sum of absolute eigenvalues)
*   **SignlessLaplacianTrace** - Trace of signless Laplacian
*   **SignlessLaplacianDeterminant** - Determinant of signless Laplacian
*   **SignlessLaplacianZeroEigenvalues** - Number of zero eigenvalues of signless Laplacian
*   **SignlessLaplacianPositiveEigenvalues** - Number of positive eigenvalues of signless Laplacian
*   **SignlessLaplacianNegativeEigenvalues** - Number of negative eigenvalues of signless Laplacian
*   **SignlessLaplacianEigenvalueVariance** - Signless Laplacian eigenvalue variance
*   **SignlessLaplacianEigenvalueSkewness** - Signless Laplacian eigenvalue skewness

#### Graph Topology and Structure Descriptors
*   **WeightedAdjacencySpectralRadius** - Weighted adjacency spectral radius
*   **MeanFirstPassageTime** - Mean first passage time approximated by pseudo-inverse of Laplacian
*   **CommuteTimeDistance** - Commute time distance (using pseudo-inverse of Laplacian)
*   **KirchhoffIndexVariance** - Kirchhoff index variance among subgraphs
*   **EffectiveGraphResistanceDistribution** - Effective graph resistance distribution (subgraphs)
*   **LocalClusteringCoefficientDistribution** - Local clustering coefficient distribution
*   **GraphRobustnessIndex** - Graph robustness index (based on spectral gap)
*   **NormalizedEstradaIndex** - Estrada index from normalized adjacency matrix
*   **GraphBipartivityIndex** - Graph bipartivity index (based on eigenvalues)
*   **SpanningTreeEntropy** - Spanning tree entropy (related to Laplacian eigenvalues)
*   **WienerIndex** - Wiener index (sum of all shortest-path distances)
*   **EstradaIndex** - Estrada index (sum of exp(eigenvalues of adjacency matrix))
*   **NumberSpanningTrees** - Number of spanning trees (via Laplacian minor)
*   **GraphEccentricity** - Graph eccentricity (maximum eccentricity of any vertex)
*   **SpectralGap** - Spectral gap (difference between largest and second-largest eigenvalues)
*   **GraphDiameter** - Graph diameter (maximum distance between any two vertices)

## Handling Invalid Molecules

desfact attempts to parse each SMILES string using RDKit. If RDKit fails to parse or sanitize a molecule, it is marked as invalid. Descriptors will not be calculated for invalid molecules.

*   By default (without `--verbose`), parsing errors from RDKit might not be explicitly printed line-by-line to the console, but they are often captured internally.
*   If `--verbose` is enabled, specific parsing/sanitization errors for individual SMILES strings are logged.
*   The output CSV will include the original row data for invalid molecules, but the descriptor columns for that row will contain placeholder values (0.0 for doubles, -1 for integers, "Error: ..." for strings).

## License

This project is likely licensed under a BSD-style license, similar to RDKit. Refer to the LICENSE file in the repository root for details.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.