#pragma once

#include "utils.hpp"
#include <string>
#include <vector>
#include <variant>
#include <memory>
#include <unordered_map>


namespace desfact {

using DescriptorValue = std::variant<double, int, std::string>;
using DescriptorMap = std::unordered_map<std::string, DescriptorValue>;

class IDescriptorCalculator {
public:
    virtual ~IDescriptorCalculator() = default;
    virtual std::string getName() const = 0;
    virtual DescriptorValue calculate(const Molecule& mol) = 0;
};

class DescriptorCalculator {
private:
    std::vector<std::shared_ptr<IDescriptorCalculator>> descriptors;
    std::unordered_map<std::string, std::shared_ptr<IDescriptorCalculator>> descriptorRegistry;

    void registerDefaultDescriptors();

public:
    DescriptorCalculator();
    void registerDescriptor(std::shared_ptr<IDescriptorCalculator> descriptor);
    std::vector<std::string> getAvailableDescriptors() const;
    void selectDescriptors(const std::vector<std::string>& names);
    std::vector<std::string> getSelectedDescriptors() const;

    DescriptorMap calculateDescriptors(const Molecule& mol);

    std::unordered_map<std::string, std::vector<DescriptorValue>>
    calculateDescriptorsForBatch(const std::vector<Molecule>& batch);
};

// Base descriptor class
class Descriptor {
protected:
    std::string name;
    std::string description;

public:
    Descriptor(const std::string& name, const std::string& description)
        : name(name), description(description) {}
    virtual ~Descriptor() = default;
    
    const std::string& getName() const { return name; }
    const std::string& getDescription() const { return description; }
    
    virtual std::variant<double, int, std::string> calculate(const Molecule& mol) const = 0;
    virtual bool isCudaSupported() const { return false; }
};

// Factory for creating descriptors
class DescriptorFactory {
private:
    std::unordered_map<std::string, std::unique_ptr<Descriptor>> descriptors;
    
public:
    DescriptorFactory();
    ~DescriptorFactory() = default;
    
    void registerDescriptor(std::unique_ptr<Descriptor> descriptor);
    const Descriptor* getDescriptor(const std::string& name) const;
    
    std::vector<std::string> getAvailableDescriptors() const;
    std::vector<std::string> getCudaEnabledDescriptors() const;
    
    // Calculate a single descriptor for a molecule
    std::variant<double, int, std::string> calculate(const std::string& descriptorName, const Molecule& mol) const;
    
    // Calculate multiple descriptors for a batch of molecules
    std::unordered_map<std::string, std::vector<std::variant<double, int, std::string>>> 
    calculateBatch(const std::vector<std::string>& descriptorNames, const MoleculeBatch& batch) const;

    // Add declaration for registerAllDescriptors
    void registerAllDescriptors();
};

// Implementation of common descriptors
class MolecularWeightDescriptor : public Descriptor {
public:
    MolecularWeightDescriptor() : Descriptor("MolWt", "Molecular Weight") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class LogPDescriptor : public Descriptor {
public:
    LogPDescriptor() : Descriptor("rdkit_LogP", "Octanol/Water Partition Coefficient") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class NumAtomsDescriptor : public Descriptor {
public:
    NumAtomsDescriptor() : Descriptor("NumAtoms", "Number of Atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

}  // namespace desfact