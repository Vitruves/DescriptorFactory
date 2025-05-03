#pragma once

#include "descriptors.hpp"
#include <GraphMol/Descriptors/MolDescriptors.h>

namespace desfact {
namespace descriptors {

class RDKitDescriptor : public Descriptor {
public:
    using DescriptorFunction = std::function<double(const RDKit::ROMol&)>;

    RDKitDescriptor(const std::string& name, const std::string& description, 
                    DescriptorFunction calcFunction);
    
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;

private:
    DescriptorFunction calcFunction_;
};

// Molecular Property Descriptors
class RDKitMWDescriptor : public RDKitDescriptor {
public:
    RDKitMWDescriptor();
};

class RDKitTPSADescriptor : public RDKitDescriptor {
public:
    RDKitTPSADescriptor();
};

class RDKitLogPDescriptor : public RDKitDescriptor {
public:
    RDKitLogPDescriptor();
};

class RDKitMRDescriptor : public RDKitDescriptor {
public:
    RDKitMRDescriptor();
};

class RDKitLipinskiHBADescriptor : public RDKitDescriptor {
public:
    RDKitLipinskiHBADescriptor();
};

class RDKitLipinskiHBDDescriptor : public RDKitDescriptor {
public:
    RDKitLipinskiHBDDescriptor();
};

class RDKitNumRotatableBondsDescriptor : public RDKitDescriptor {
public:
    RDKitNumRotatableBondsDescriptor();
};

class RDKitNumHBondAcceptorsDescriptor : public RDKitDescriptor {
public:
    RDKitNumHBondAcceptorsDescriptor();
};

class RDKitNumHBondDonorsDescriptor : public RDKitDescriptor {
public:
    RDKitNumHBondDonorsDescriptor();
};

class RDKitNumHeavyAtomsDescriptor : public RDKitDescriptor {
public:
    RDKitNumHeavyAtomsDescriptor();
};

class RDKitNumRingsDescriptor : public RDKitDescriptor {
public:
    RDKitNumRingsDescriptor();
};

class RDKitNumAromaticRingsDescriptor : public RDKitDescriptor {
public:
    RDKitNumAromaticRingsDescriptor();
};

class RDKitNumAliphaticRingsDescriptor : public RDKitDescriptor {
public:
    RDKitNumAliphaticRingsDescriptor();
};

class RDKitFractionCsp3Descriptor : public RDKitDescriptor {
public:
    RDKitFractionCsp3Descriptor();
};

class RDKitChiralCentersDescriptor : public RDKitDescriptor {
public:
    RDKitChiralCentersDescriptor();
};

class RDKitLabuteASADescriptor : public RDKitDescriptor {
public:
    RDKitLabuteASADescriptor();
};

class RDKitBalabanJDescriptor : public RDKitDescriptor {
public:
    RDKitBalabanJDescriptor();
};

class RDKitBertzCTDescriptor : public RDKitDescriptor {
public:
    RDKitBertzCTDescriptor();
};

class RDKitHallKierAlphaDescriptor : public RDKitDescriptor {
public:
    RDKitHallKierAlphaDescriptor();
};

class RDKitKappaDescriptor1 : public RDKitDescriptor {
public:
    RDKitKappaDescriptor1();
};

class RDKitKappaDescriptor2 : public RDKitDescriptor {
public:
    RDKitKappaDescriptor2();
};

class RDKitKappaDescriptor3 : public RDKitDescriptor {
public:
    RDKitKappaDescriptor3();
};

// Chi descriptors
class RDKitChi0vDescriptor : public RDKitDescriptor {
public:
    RDKitChi0vDescriptor();
};

class RDKitChi1vDescriptor : public RDKitDescriptor {
public:
    RDKitChi1vDescriptor();
};

class RDKitChi2vDescriptor : public RDKitDescriptor {
public:
    RDKitChi2vDescriptor();
};

class RDKitChi3vDescriptor : public RDKitDescriptor {
public:
    RDKitChi3vDescriptor();
};

class RDKitChi4vDescriptor : public RDKitDescriptor {
public:
    RDKitChi4vDescriptor();
};

// QED descriptor
class RDKitQEDDescriptor : public RDKitDescriptor {
public:
    RDKitQEDDescriptor();
};

// SASA descriptor
class RDKitSASADescriptor : public RDKitDescriptor {
public:
    RDKitSASADescriptor();
};

// PEOE VSA descriptors
class RDKitPEOE_VSA1Descriptor : public RDKitDescriptor {
public:
    RDKitPEOE_VSA1Descriptor();
};

// Exact MW from Formula
class RDKitExactMWDescriptor : public RDKitDescriptor {
public:
    RDKitExactMWDescriptor();
};

// MURCKO Framework Count
class RDKitMurckoScaffoldCountDescriptor : public RDKitDescriptor {
public:
    RDKitMurckoScaffoldCountDescriptor();
};

// Ring Count by Size
class RDKitRingCount5Descriptor : public RDKitDescriptor {
public:
    RDKitRingCount5Descriptor();
};

// Spiro atom count
class RDKitSpiroAtomCountDescriptor : public RDKitDescriptor {
public:
    RDKitSpiroAtomCountDescriptor();
};

// Bridgehead atom count
class RDKitBridgeheadAtomCountDescriptor : public RDKitDescriptor {
public:
    RDKitBridgeheadAtomCountDescriptor();
};

} // namespace descriptors
} // namespace desfact
