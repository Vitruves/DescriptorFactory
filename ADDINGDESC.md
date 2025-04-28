MEMORANDUM

TO: me
FROM: me 
DATE: 2025
SUBJECT: Guidelines for Adding New Descriptors Efficiently

This memo describes the standard procedure for implementing new molecular descriptors within the `desfact` project architecture. By following this process, we can add new functionality with minimal changes to core logic, primarily modifying only the new descriptor file(s) and the central registration point. This approach reduces development overhead and maintains code organization.

**Core Architecture for Extensibility**

The `desfact` project uses a **Descriptor Factory** pattern (`DescriptorFactory` class in `include/descriptors.hpp` and `src/descriptors.cpp`) to manage descriptors. Each specific descriptor is implemented as a class inheriting from a base `Descriptor` class (also in `include/descriptors.hpp`). This factory pattern allows us to register new descriptor classes without modifying the main descriptor calculation loop or file processing logic.

Furthermore, common descriptor categories (like sum-based or fractional descriptors) utilize base classes (`SumDescriptor`, `FractionalDescriptor`, `StringDescriptor`) which provide shared helper functions. This reduces code duplication and promotes consistency.

**Shared Functions and Resources**

New descriptor implementations should utilize the following shared resources where appropriate:

1.  **`Molecule` Object:**
    *   Input to the primary `calculate` method (`std::variant<double, int, std::string> calculate(const Molecule& mol) const override;`).
    *   Provides access to the underlying RDKit `ROMol` object (`mol.getMolecule()`) and basic information (`mol.isValid()`, `mol.getSmiles()`, `mol.getNumAtoms()`, `mol.getNumBonds()`, `mol.getOriginalSmiles()`, `mol.getErrorMessage()`).
    *   Always check `mol.isValid()` at the beginning of your `calculate` method and return an appropriate default/error value if false (e.g., `0.0`, `-1`, `"Error: Invalid Molecule"`).

2.  **`Descriptor` Base Class:**
    *   You must inherit from this class (or a derived base like `SumDescriptor`).
    *   Defines the virtual interface including `getName()`, `getDescription()`, and the crucial `calculate()` method.

3.  **Category Base Classes (e.g., `SumDescriptor`, `FractionalDescriptor`, `StringDescriptor`):**
    *   Inherit from these if your descriptor fits the category, allowing you to reuse common helper logic.
    *   **`SumDescriptor`:** Provides static helpers for iterating atoms/bonds and summing values:
        *   `static double calcAtomicSum(const Molecule& mol, const std::function<double(const RDKit::Atom*)>& valueFunction, bool normalize = false);` - Iterates atoms, applies `valueFunction`, and sums results. Can normalize by atom count.
        *   `static double calcBondSum(const Molecule& mol, const std::function<double(const RDKit::Bond*)>& valueFunction, bool normalize = false);` - Iterates bonds, applies `valueFunction`, and sums results. Can normalize by bond count.
        *   Also provides static helpers for accessing common atomic properties (electronegativity, radii, polarizability, etc.) via `getAtomElectronegativity`, `getAtomCovalentRadius`, etc. and basic atom/bond checks (`isHeteroatom`, `isPolarBond`).
    *   **`FractionalDescriptor`:** Provides static helpers for calculating fractions of atoms or bonds matching a predicate:
        *   `static double calcAtomicFraction(const Molecule& mol, const std::function<bool(const RDKit::Atom*)>& predicate);` - Counts atoms where `predicate` is true, returns count / total atoms.
        *   `static double calcBondFraction(const Molecule& mol, const std::function<bool(const RDKit::Bond*)>& predicate);` - Counts bonds where `predicate` is true, returns count / total bonds.
        *   Also provides static helpers for common atom checks (`isElement`, `isHalogen`, `isHeteroatom`, `isMetalAtom`, etc.) and properties, similar to `SumDescriptor`.
    *   **`StringDescriptor`:** Designed for descriptors calculated solely from the SMILES string.
        *   Requires implementing `virtual std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const = 0;`. The base class `calculate` calls this.
        *   Provides static helpers for string analysis (`calcAsciiSum`, `calcEntropy`, `getLongestRun`, `getSymmetryScore`, `computeStringStats`, etc.).

4.  **RDKit API:** Access the underlying RDKit `ROMol` object via `mol.getMolecule().get()`. Use RDKit functions (e.g., `mol->getNumAtoms()`, `atom->getAtomicNum()`, `bond->getBondType()`, `RDKit::MolOps`, `RDKit::Descriptors`, `RDKit::SubstructMatch`, `RDKit::SmartsToMol`, `RDKit::PeriodicTable`).

5.  **Element Properties:** The `ElementProperties` namespace in `src/descriptors/fractional.cpp` and the anonymous namespace maps in `src/descriptors/sum.cpp` provide common physical properties for elements. Access is typically via the static helper functions in the base classes.

**Step-by-Step Implementation of a New Descriptor**

Let's say we want to add a new descriptor called "ExampleDescriptor" that calculates the sum of the number of explicit hydrogens on heavy (non-Hydrogen) atoms.

1.  **Choose a Base Class:** This descriptor sums a property over atoms, specifically heavy atoms. `SumDescriptor` is a good fit as it provides `calcAtomicSum` and potentially atom-checking helpers.

2.  **Create the Header File:**
    *   Create a new file: `include/descriptors/example.hpp`
    *   Add standard include guards (`#pragma once`).
    *   Include the base class header (`#include "descriptors/sum.hpp"` in this case).
    *   Include any necessary RDKit headers for types used in the class definition (e.g., `<GraphMol/Atom.h>`, `<GraphMol/ROMol.h>`, etc. - often covered by the base header).
    *   Define your new class, inheriting from the chosen base.
    *   Declare the constructor and the `calculate` method override.

    ```cpp
    // include/descriptors/example.hpp
    #pragma once

    #include "descriptors/sum.hpp"
    #include <GraphMol/Atom.h> // Might be needed if using Atom pointers directly in signature, though calcAtomicSum hides this

    namespace desfact {
    namespace descriptors {

    class ExampleDescriptor : public SumDescriptor {
    public:
        ExampleDescriptor(); // Declare constructor
        // Override the virtual calculate method from the base Descriptor class
        // The return type MUST be std::variant<double, int, std::string>
        std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
    };

    } // namespace descriptors
    } // namespace desfact
    ```

3.  **Create the Source File:**
    *   Create a new file: `src/descriptors/example.cpp`
    *   Include your new header (`#include "descriptors/example.hpp"`).
    *   Include any other necessary RDKit headers for implementation details (e.g., `<GraphMol/GraphMol.h>`, `<GraphMol/AtomIterators.h>`).
    *   Implement the constructor and the `calculate` method.
    *   Use the `Molecule` object and appropriate helper functions or RDKit calls. Handle the `mol.isValid()` case.

    ```cpp
    // src/descriptors/example.cpp
    #include "descriptors/example.hpp"
    #include <GraphMol/GraphMol.h>
    #include <GraphMol/AtomIterators.h>
    #include <numeric> // For std::accumulate, maybe
    #include "utils.hpp" // For globalLogger, etc.

    namespace desfact {
    namespace descriptors {

    // Implement the constructor
    ExampleDescriptor::ExampleDescriptor() 
        : SumDescriptor("SumHHeavy", "Sum of explicit hydrogens on heavy atoms") 
    {
        // The base class constructor sets the name and description.
    }

    // Implement the calculate method
    std::variant<double, int, std::string> ExampleDescriptor::calculate(const Molecule& mol) const {
        // 1. Handle invalid molecule input
        if (!mol.isValid() || !mol.getMolecule()) {
            globalLogger.debug("ExampleDescriptor: Invalid molecule input.");
            return 0; // Or std::string("Error: Invalid Input")
        }

        // 2. Access RDKit molecule
        const RDKit::ROMol* rdkMol = mol.getMolecule().get();
        if (!rdkMol) {
             globalLogger.error("ExampleDescriptor: RDKit molecule is null despite Molecule being valid?");
             return std::string("Error: Internal RDKit Null");
        }

        // 3. Implement calculation logic
        // Option A: Manual iteration
        int sumHydrogens = 0;
        for (const RDKit::Atom* atom : rdkMol->atoms()) {
            // Check if atom is heavy (AtomicNum > 1)
            if (atom->getAtomicNum() > 1) {
                sumHydrogens += atom->getNumExplicitHs();
            }
        }
        return sumHydrogens;

        // Option B: Using SumDescriptor::calcAtomicSum helper (less direct here,
        // as it sums the *result* of the function, not a count, but possible)
        // This helper is usually for summing a property value (like EN, radius).
        // For simple counts or sums of integer properties like H count,
        // manual iteration is often clearer or necessary depending on the
        // exact calculation. Let's stick to the manual one for clarity here.

        // Option C: Using FractionalDescriptor::calcAtomicFraction (if calculating a fraction)
        // Example: Fraction of heavy atoms with >0 explicit hydrogens
        // return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        //     return atom->getAtomicNum() > 1 && atom->getNumExplicitHs() > 0;
        // });
    }

    } // namespace descriptors
    } // namespace desfact
    ```
    *   **Note on Error Handling:** Within the `calculate` method, focus on handling RDKit-specific errors during the calculation if they occur (e.g., trying to access properties that require sanitization if the molecule wasn't sanitized, or issues with specific RDKit functions). Wrap potential problematic RDKit calls in `try...catch`. Return `"Error: ..."` strings or default numerical values (`0.0`, `-1`, `NaN`) for errors. The `Molecule::isValid()` check handles basic parsing errors before your calculation starts.

4.  **Register the New Descriptor:**
    *   Open `src/descriptors.cpp`.
    *   Include your new header at the top: `#include "descriptors/example.hpp"`.
    *   Find the `DescriptorFactory::DescriptorFactory()` constructor.
    *   Add a line to register an instance of your new descriptor using `registerDescriptor` and `std::make_unique`.

    ```cpp
    // src/descriptors.cpp
    #include "descriptors.hpp"
    #include "descriptors/fractional.hpp"
    #include "descriptors/sum.hpp"
    #include "descriptors/strings.hpp"
    #include "descriptors/example.hpp" // <--- Add this line
    // ... other includes ...

    namespace desfact {
    // ... other descriptor implementations ...

    DescriptorFactory::DescriptorFactory() {
        // Basic descriptors (defined in descriptors.hpp)
        registerDescriptor(std::make_unique<MolecularWeightDescriptor>());
        registerDescriptor(std::make_unique<LogPDescriptor>());
        registerDescriptor(std::make_unique<NumAtomsDescriptor>());

        // Fractional descriptors (from fractional.hpp/cpp)
        registerDescriptor(std::make_unique<descriptors::FcCDescriptor>());
        registerDescriptor(std::make_unique<descriptors::FcFDescriptor>());
        // ... register all other FractionalDescriptors ...
        registerDescriptor(std::make_unique<descriptors::FcENAboveMoleculeAvgDescriptor>());

        // Sum descriptors (from sum.hpp/cpp)
        registerDescriptor(std::make_unique<descriptors::SumENDescriptor>());
        registerDescriptor(std::make_unique<descriptors::SumCovalentRadiiDescriptor>());
        // ... register all other SumDescriptors ...
        registerDescriptor(std::make_unique<descriptors::SumPolENBondedDescriptor>());

        // String descriptors (from strings.hpp/cpp)
        registerDescriptor(std::make_unique<descriptors::AsciiSumDescriptor>());
        registerDescriptor(std::make_unique<descriptors::AsciiAverageDescriptor>());
        // ... register all other StringDescriptors ...
        registerDescriptor(std::make_unique<descriptors::RingComplexityDescriptor>());

        // *** Add your new descriptor here ***
        registerDescriptor(std::make_unique<descriptors::ExampleDescriptor>()); // <--- Add this line

        globalLogger.info("Initialized descriptor factory with " +
                         std::to_string(descriptors.size()) + " descriptors");
    }

    // ... calculateBatch implementation ...

    } // namespace desfact
    ```

5.  **Update CMakeLists.txt:**
    *   Open the project's `CMakeLists.txt`.
    *   Find the `set(SOURCES ...)` block (around line 156).
    *   Add your new source file (`src/descriptors/example.cpp`) to this list.

    ```cmake
    # CMakeLists.txt
    # ... other CMake configuration ...

    # Updated source list
    set(SOURCES
        src/utils.cpp
        src/io.cpp
        src/descriptors.cpp
        src/descriptors/fractional.cpp
        src/descriptors/sum.cpp
        src/descriptors/strings.cpp
        src/descriptors/example.cpp # <--- Add this line
    )

    # ... rest of CMakeLists.txt ...
    ```
    *   You do *not* need to add the new header file to CMake; `#include` directives handle that, and the `include_directories` command ensures the compiler can find files in the `include/` path.

6.  **Build and Test:** Compile the project. Your new descriptor ("SumHHeavy" in this example) should now be available when you run `desfact_main`. You can verify it by listing descriptors (`desfact_main --list`) and then calculating it on a test CSV file.

**Summary of Changes:**

To add a new descriptor ("ExampleDescriptor"), the changes required are localized to:

*   `include/descriptors/example.hpp` (new file)
*   `src/descriptors/example.cpp` (new file)
*   `src/descriptors.cpp` (one line to register)
*   `CMakeLists.txt` (one line to add the new `.cpp` file)

This pattern effectively adds new functionality ("without introducing overhead") compared to modifying fundamental parts of the calculation loop or data handling for each new descriptor. The core processing pipeline (`main.cpp`, `CsvIO`, `DescriptorFactory::calculateBatch`) remains largely untouched, automatically discovering and using the newly registered descriptor.

Please follow these guidelines for future descriptor additions.
