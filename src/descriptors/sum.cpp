#include "descriptors/sum.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/Descriptors/MolSurf.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <mutex>
#include <unordered_map>
#include <memory>
#include <cmath>

namespace desfact {
namespace descriptors {

namespace {

struct MolData {
    std::shared_ptr<const RDKit::ROMol> mol;
    std::vector<const RDKit::Atom*> atoms;
    std::vector<const RDKit::Bond*> bonds;
    RDKit::RingInfo* ringInfo;
};

static std::unordered_map<std::shared_ptr<const RDKit::ROMol>,
                          std::unique_ptr<MolData>> s_molDataCache;
static std::mutex s_molDataMutex;

static const MolData* getMolData(const Molecule& mol) {
    if(!mol.isValid() || !mol.getMolecule()) return nullptr;
    auto rdkMol = mol.getMolecule();
    {
        std::lock_guard<std::mutex> lock(s_molDataMutex);
        auto it = s_molDataCache.find(rdkMol);
        if(it != s_molDataCache.end()) {
            return it->second.get();
        }
    }
    auto md = std::make_unique<MolData>();
    md->mol = rdkMol;
    for(auto a : md->mol->atoms()) md->atoms.push_back(a);
    for(auto b : md->mol->bonds()) md->bonds.push_back(b);
    md->ringInfo = md->mol->getRingInfo();
    {
        std::lock_guard<std::mutex> lock(s_molDataMutex);
        if(!md->ringInfo->isInitialized()) {
            RDKit::MolOps::findSSSR(*const_cast<RDKit::ROMol*>(md->mol.get()));
        }
        if(!md->ringInfo->isInitialized()) {
            return nullptr;
        }
        s_molDataCache[md->mol] = std::move(md);
        return s_molDataCache[rdkMol].get();
    }
}

static const std::unordered_map<int, double> s_enMap = {
    {1, 2.20}, {6, 2.55}, {7, 3.04}, {8, 3.44}, {9, 3.98},
    {15, 2.19}, {16, 2.58}, {17, 3.16}, {35, 2.96}, {53, 2.66}
};

static const std::unordered_map<int, double> s_covalentRadius = {
    {1, 0.31}, {6, 0.76}, {7, 0.71}, {8, 0.66}, {9, 0.57},
    {15, 1.07}, {16, 1.05}, {17, 1.02}, {35, 1.20}, {53, 1.39}
};

static const std::unordered_map<int, double> s_vdwVolume = {
    {1, 7.2}, {6, 20.6}, {7, 15.6}, {8, 14.0}, {9, 13.3},
    {15, 24.4}, {16, 24.4}, {17, 22.5}, {35, 26.5}, {53, 32.9}
};

static const std::unordered_map<int, double> s_polarizability = {
    {1, 0.667}, {6, 1.76}, {7, 1.10}, {8, 0.802}, {9, 0.557},
    {15, 3.63}, {16, 2.90}, {17, 2.18}, {35, 3.05}, {53, 5.35}
};

static const std::unordered_map<int, double> s_ionizationEnergy = {
    {1, 13.6}, {6, 11.3}, {7, 14.5}, {8, 13.6}, {9, 17.4},
    {15, 10.5}, {16, 10.4}, {17, 13.0}, {35, 11.8}, {53, 10.5}
};

static const std::unordered_map<int, double> s_electronAffinity = {
    {1, 0.75}, {6, 1.26}, {7, 0.07}, {8, 1.46}, {9, 3.40},
    {15, 0.75}, {16, 2.08}, {17, 3.62}, {35, 3.36}, {53, 3.06}
};

static const std::unordered_set<int> s_metals = {
    3,4,11,12,13,19,20,21,22,23,24,25,26,27,28,29,30,31,
    37,38,39,40,41,42,43,44,45,46,47,48,49,50,
    55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,
    72,73,74,75,76,77,78,79,80,81,82,83,84,
    87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103
};

static bool isHeteroatom(const RDKit::Atom* atom) {
    if(!atom) return false;
    int n = atom->getAtomicNum();
    return n!=6 && n!=1;
}

static bool isHalogen(const RDKit::Atom* atom) {
    if(!atom) return false;
    int n = atom->getAtomicNum();
    return n==9||n==17||n==35||n==53||n==85;
}

static bool isMetalAtom(const RDKit::Atom* atom) {
    if(!atom) return false;
    return s_metals.find(atom->getAtomicNum())!=s_metals.end();
}

static double getAtomEN(const RDKit::Atom* atom) {
    if(!atom) return 0.0;
    auto it = s_enMap.find(atom->getAtomicNum());
    return it!=s_enMap.end()?it->second:2.20;
}

static double getAtomRcov(const RDKit::Atom* atom) {
    if(!atom) return 0.0;
    auto it = s_covalentRadius.find(atom->getAtomicNum());
    return it!=s_covalentRadius.end()?it->second:0.75;
}

static double getAtomVdW(const RDKit::Atom* atom) {
    if(!atom) return 0.0;
    auto it = s_vdwVolume.find(atom->getAtomicNum());
    return it!=s_vdwVolume.end()?it->second:20.0;
}

static double getAtomPol(const RDKit::Atom* atom) {
    if(!atom) return 0.0;
    auto it = s_polarizability.find(atom->getAtomicNum());
    return it!=s_polarizability.end()?it->second:1.0;
}

static double getAtomIE(const RDKit::Atom* atom) {
    if(!atom) return 0.0;
    auto it = s_ionizationEnergy.find(atom->getAtomicNum());
    return it!=s_ionizationEnergy.end()?it->second:12.0;
}

static double getAtomEA(const RDKit::Atom* atom) {
    if(!atom) return 0.0;
    auto it = s_electronAffinity.find(atom->getAtomicNum());
    return it!=s_electronAffinity.end()?it->second:1.0;
}

static bool isPolarBond(const RDKit::Bond* bond) {
    if(!bond) return false;
    double en1 = getAtomEN(bond->getBeginAtom());
    double en2 = getAtomEN(bond->getEndAtom());
    return std::fabs(en1 - en2) >= 0.5;
}

} // end anonymous namespace

double SumDescriptor::calcAtomicSum(const Molecule& mol,
    const std::function<double(const RDKit::Atom*)>& f, bool norm) {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s = 0.0;
    for(auto a : data->atoms) s += f(a);
    if(norm && !data->atoms.empty()) s /= data->atoms.size();
    return s;
}

double SumDescriptor::calcBondSum(const Molecule& mol,
    const std::function<double(const RDKit::Bond*)>& f, bool norm) {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s = 0.0;
    for(auto b : data->bonds) s += f(b);
    if(norm && !data->bonds.empty()) s /= data->bonds.size();
    return s;
}

double SumDescriptor::getAtomElectronegativity(const RDKit::Atom* a) {
    return getAtomEN(a);
}

double SumDescriptor::getAtomCovalentRadius(const RDKit::Atom* a) {
    return getAtomRcov(a);
}

double SumDescriptor::getAtomVdWVolume(const RDKit::Atom* a) {
    return getAtomVdW(a);
}

double SumDescriptor::getAtomPolarizability(const RDKit::Atom* a) {
    return getAtomPol(a);
}

double SumDescriptor::getAtomIonizationEnergy(const RDKit::Atom* a) {
    return getAtomIE(a);
}

double SumDescriptor::getAtomElectronAffinity(const RDKit::Atom* a) {
    return getAtomEA(a);
}

bool SumDescriptor::isHeteroatom(const RDKit::Atom* a) {
    return desfact::descriptors::isHeteroatom(a);
}

bool SumDescriptor::isHalogen(const RDKit::Atom* a) {
    return desfact::descriptors::isHalogen(a);
}

bool SumDescriptor::isMetalAtom(const RDKit::Atom* a) {
    return desfact::descriptors::isMetalAtom(a);
}

bool SumDescriptor::isPolarBond(const RDKit::Bond* b) {
    return desfact::descriptors::isPolarBond(b);
}

std::variant<double,int,std::string> SumENDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicSum(mol, [&](auto a){return getAtomEN(a);});
}

std::variant<double,int,std::string> SumCovalentRadiiDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicSum(mol, [&](auto a){return getAtomRcov(a);});
}

std::variant<double,int,std::string> SumVdWVolumeDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicSum(mol, [&](auto a){return getAtomVdW(a);});
}

std::variant<double,int,std::string> SumPolarizabilityDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicSum(mol, [&](auto a){return getAtomPol(a);});
}

std::variant<double,int,std::string> SumIonizationEnergyDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicSum(mol, [&](auto a){return getAtomIE(a);});
}

std::variant<double,int,std::string> SumElectronAffinityDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicSum(mol, [&](auto a){return getAtomEA(a);});
}

std::variant<double,int,std::string> SumAtomsDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0;
    return (int)data->atoms.size();
}

std::variant<double,int,std::string> SumHeavyAtomsDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicSum(mol, [](auto a){return a->getAtomicNum()>1?1.0:0.0;});
}

std::variant<double,int,std::string> SumHeteroatomsDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicSum(mol, [&](auto a){return isHeteroatom(a)?1.0:0.0;});
}

std::variant<double,int,std::string> SumHalogensDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicSum(mol, [&](auto a){return isHalogen(a)?1.0:0.0;});
}

std::variant<double,int,std::string> SumChargedAtomsDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicSum(mol, [](auto a){return a->getFormalCharge()!=0?1.0:0.0;});
}

std::variant<double,int,std::string> SumBondsDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0;
    return (int)data->bonds.size();
}

std::variant<double,int,std::string> SumDoubleBondsDescriptor::calculate(const Molecule& mol) const {
    return calcBondSum(mol, [&](auto b){return b->getBondType()==RDKit::Bond::DOUBLE?1.0:0.0;});
}

std::variant<double,int,std::string> SumTripleBondsDescriptor::calculate(const Molecule& mol) const {
    return calcBondSum(mol, [&](auto b){return b->getBondType()==RDKit::Bond::TRIPLE?1.0:0.0;});
}

std::variant<double,int,std::string> SumAromaticBondsDescriptor::calculate(const Molecule& mol) const {
    return calcBondSum(mol, [&](auto b){return b->getIsAromatic()?1.0:0.0;});
}

std::variant<double,int,std::string> SumPolarBondsDescriptor::calculate(const Molecule& mol) const {
    return calcBondSum(mol, [&](auto b){return isPolarBond(b)?1.0:0.0;});
}

std::variant<double,int,std::string> SumUnpolarBondsDescriptor::calculate(const Molecule& mol) const {
    return calcBondSum(mol, [&](auto b){return !isPolarBond(b)?1.0:0.0;});
}

std::variant<double,int,std::string> SumSp3BondsDescriptor::calculate(const Molecule& mol) const {
    return calcBondSum(mol, [](auto b){
        auto a1 = b->getBeginAtom();
        auto a2 = b->getEndAtom();
        return (a1->getHybridization()==RDKit::Atom::SP3 && a2->getHybridization()==RDKit::Atom::SP3)?1.0:0.0;
    });
}

std::variant<double,int,std::string> SumSp2BondsDescriptor::calculate(const Molecule& mol) const {
    return calcBondSum(mol, [](auto b){
        auto a1 = b->getBeginAtom();
        auto a2 = b->getEndAtom();
        return (a1->getHybridization()==RDKit::Atom::SP2 && a2->getHybridization()==RDKit::Atom::SP2)?1.0:0.0;
    });
}

std::variant<double,int,std::string> SumRingsDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0;
    return (int)data->ringInfo->numRings();
}

std::variant<double,int,std::string> SumAromaticRingsDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0;
    int c=0;
    for(unsigned int i=0;i<data->ringInfo->numRings();++i) {
        const auto& ring = data->ringInfo->atomRings()[i];
        bool aro=true;
        for(auto idx : ring) {
            if(!data->mol->getAtomWithIdx(idx)->getIsAromatic()){
                aro=false;break;
            }
        }
        if(aro) c++;
    }
    return c;
}

std::variant<double,int,std::string> SumRotatableBondsDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0;
    return (int)RDKit::Descriptors::calcNumRotatableBonds(*data->mol);
}

std::variant<double,int,std::string> SumBridgeAtomsDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0;
    int c=0;
    for(auto a : data->atoms) {
        if(data->ringInfo->numAtomRings(a->getIdx())>=2) c++;
    }
    return c;
}

std::variant<double,int,std::string> SumRingAtomsDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0;
    int c=0;
    for(auto a : data->atoms) {
        if(data->ringInfo->numAtomRings(a->getIdx())>0) c++;
    }
    return c;
}

std::variant<double,int,std::string> SumENMWDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto a : data->atoms) {
        double en=getAtomEN(a);
        double mw=RDKit::PeriodicTable::getTable()->getAtomicWeight(a->getAtomicNum());
        s+=en*mw;
    }
    return s;
}

std::variant<double,int,std::string> SumPolMWDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto a : data->atoms) {
        double p=getAtomPol(a);
        double mw=RDKit::PeriodicTable::getTable()->getAtomicWeight(a->getAtomicNum());
        s+=p*mw;
    }
    return s;
}

std::variant<double,int,std::string> SumENRcovDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto a : data->atoms) {
        double en = getAtomEN(a);
        double rc = getAtomRcov(a);
        s += en*rc;
    }
    return s;
}

std::variant<double,int,std::string> SumENPolDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto a : data->atoms) {
        double en=getAtomEN(a);
        double p=getAtomPol(a);
        s+=en*p;
    }
    return s;
}

std::variant<double,int,std::string> SumSp3AtomsDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicSum(mol, [](auto a){return a->getHybridization()==RDKit::Atom::SP3?1.0:0.0;});
}

std::variant<double,int,std::string> SumSp2AtomsDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicSum(mol, [](auto a){return a->getHybridization()==RDKit::Atom::SP2?1.0:0.0;});
}

std::variant<double,int,std::string> SumHBDonorsDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0;
    return (int)RDKit::Descriptors::calcNumHBD(*data->mol);
}

std::variant<double,int,std::string> SumHBAcceptorsDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0;
    return (int)RDKit::Descriptors::calcNumHBA(*data->mol);
}

std::variant<double,int,std::string> SumAmideGroupsDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0;
    std::string pat="[NX3][CX3](=[OX1])";
    RDKit::ROMol* query=RDKit::SmartsToMol(pat);
    if(!query) return 0;
    std::vector<std::vector<std::pair<int,int>>> matches;
    int n=RDKit::SubstructMatch(*data->mol,*query,matches);
    delete query;
    return n;
}

std::variant<double,int,std::string> SumCarboxylGroupsDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0;
    std::string pat="[CX3](=O)[OX2H1]";
    RDKit::ROMol* query=RDKit::SmartsToMol(pat);
    if(!query) return 0;
    std::vector<std::vector<std::pair<int,int>>> matches;
    int n=RDKit::SubstructMatch(*data->mol,*query,matches);
    delete query;
    return n;
}

std::variant<double,int,std::string> SumElectronegativityRingDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto a : data->atoms) {
        if(data->ringInfo->numAtomRings(a->getIdx())>0) {
            s+=getAtomEN(a);
        }
    }
    return s;
}

std::variant<double,int,std::string> SumElectronegativityBondedDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto b : data->bonds) {
        s += getAtomEN(b->getBeginAtom()) + getAtomEN(b->getEndAtom());
    }
    return s;
}

std::variant<double,int,std::string> SumAtomicRadiusHeavyDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicSum(mol, [&](auto a){
        return a->getAtomicNum()>1?getAtomRcov(a):0.0;
    });
}

std::variant<double,int,std::string> SumPolzENRatioDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto a : data->atoms) {
        double en=getAtomEN(a);
        double p=getAtomPol(a);
        if(std::abs(en) > 1e-9) s+=(p/en);
        else s+= (p > 1e-9 ? std::numeric_limits<double>::quiet_NaN() : 0.0);
    }
    if (std::isnan(s)) return 0.0;
    return s;
}

std::variant<double,int,std::string> SumIEENRatioDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto a : data->atoms) {
        double en=getAtomEN(a);
        double ie=getAtomIE(a);
        if(std::abs(en) > 1e-9) s+=(ie/en);
        else s+= (ie > 1e-9 ? std::numeric_limits<double>::quiet_NaN() : 0.0);
    }
    if (std::isnan(s)) return 0.0;
    return s;
}

std::variant<double,int,std::string> SumEAMWENDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto a : data->atoms) {
        double ea=getAtomEA(a);
        double en=getAtomEN(a);
        double mw=RDKit::PeriodicTable::getTable()->getAtomicWeight(a->getAtomicNum());
        s+=ea*mw*en;
    }
    return s;
}

std::variant<double,int,std::string> SumSp2CAtomsDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicSum(mol, [](auto a){
        return (a->getAtomicNum()==6 && a->getHybridization()==RDKit::Atom::SP2)?1.0:0.0;
    });
}

std::variant<double,int,std::string> SumSpAtomsENDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto a : data->atoms) {
        if(a->getHybridization()==RDKit::Atom::SP) {
            s+=getAtomEN(a);
        }
    }
    return s;
}

std::variant<double,int,std::string> SumFormalChargeAbsDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicSum(mol, [](auto a){
        return std::fabs((double)a->getFormalCharge());
    });
}

std::variant<double,int,std::string> SumFormalChargeHeavyDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicSum(mol, [](auto a){
        return a->getAtomicNum()>1?(double)a->getFormalCharge():0.0;
    });
}

std::variant<double,int,std::string> SumFormalChargePositiveDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicSum(mol, [](auto a){
        int c=a->getFormalCharge();
        return c>0?(double)c:0.0;
    });
}

std::variant<double,int,std::string> SumFormalChargeNegativeDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicSum(mol, [](auto a){
        int c=a->getFormalCharge();
        return c<0?std::fabs((double)c):0.0;
    });
}

std::variant<double,int,std::string> SumENMWRatioDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto a : data->atoms) {
        double en=getAtomEN(a);
        double mw=RDKit::PeriodicTable::getTable()->getAtomicWeight(a->getAtomicNum());
        if(mw > 1e-9) s+=(en/mw);
        else s+= (std::abs(en) > 1e-9 ? std::numeric_limits<double>::quiet_NaN() : 0.0);
    }
    if (std::isnan(s)) return 0.0;
    return s;
}

std::variant<double,int,std::string> SumPolzRingDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto a : data->atoms) {
        if(data->ringInfo->numAtomRings(a->getIdx())>0) {
            s+=getAtomPol(a);
        }
    }
    return s;
}

std::variant<double,int,std::string> SumIEBondedDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto a : data->atoms) {
        if(a->getDegree()>=1) {
            s+=getAtomIE(a);
        }
    }
    return s;
}

std::variant<double,int,std::string> SumEAHeavyBondedDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto a : data->atoms) {
        if(a->getAtomicNum()>1 && a->getDegree()>=1) {
            s+=getAtomEA(a);
        }
    }
    return s;
}

std::variant<double,int,std::string> SumSp3ENWeightedDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto a : data->atoms) {
        if(a->getHybridization()==RDKit::Atom::SP3) {
            s+=getAtomEN(a);
        }
    }
    return s;
}

std::variant<double,int,std::string> SumHeteroMWENDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto a : data->atoms) {
        if(isHeteroatom(a)) {
            double en=getAtomEN(a);
            double mw=RDKit::PeriodicTable::getTable()->getAtomicWeight(a->getAtomicNum());
            s+=mw*en;
        }
    }
    return s;
}

std::variant<double,int,std::string> SumENSp2HeavyDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto a : data->atoms) {
        if(a->getAtomicNum()>1 && a->getHybridization()==RDKit::Atom::SP2) {
            s+=getAtomEN(a);
        }
    }
    return s;
}

std::variant<double,int,std::string> SumENRadicalAtomsDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto a : data->atoms) {
        if(a->getNumRadicalElectrons()>0) {
            s+=getAtomEN(a);
        }
    }
    return s;
}

std::variant<double,int,std::string> SumOxHeteroMWDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto a : data->atoms) {
        if(isHeteroatom(a)) {
            double ox=0.0;
            int an=a->getAtomicNum();
            int val=a->getTotalValence();
            int fc=a->getFormalCharge();
            if(an==7) ox=-3+val+fc;
            else if(an==8) ox=-2+val+fc;
            else if(an==16) ox=-2+val+fc;
            else if(an==9) ox=-1+val+fc;
            else if(an==17) ox=-1+val+fc;
            else if(an==35) ox=-1+val+fc;
            else if(an==53) ox=-1+val+fc;
            double mw=RDKit::PeriodicTable::getTable()->getAtomicWeight(an);
            s+=std::fabs(ox)*mw;
        }
    }
    return s;
}

std::variant<double,int,std::string> SumENRingHeavyDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto a : data->atoms) {
        if(a->getAtomicNum()>1 && data->ringInfo->numAtomRings(a->getIdx())>0) {
            s+=getAtomEN(a);
        }
    }
    return s;
}

std::variant<double,int,std::string> SumPolENBondedDescriptor::calculate(const Molecule& mol) const {
    auto data = getMolData(mol);
    if(!data) return 0.0;
    double s=0.0;
    for(auto a : data->atoms) {
        if(a->getDegree()>=1) {
            s += getAtomPol(a)*getAtomEN(a);
        }
    }
    return s;
}

} // namespace descriptors
} // namespace desfact 