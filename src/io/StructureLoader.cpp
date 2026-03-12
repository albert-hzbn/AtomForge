#include "StructureLoader.h"

#include <openbabel3/openbabel/obconversion.h>
#include <openbabel3/openbabel/mol.h>
#include <openbabel3/openbabel/atom.h>
#include <openbabel3/openbabel/elements.h>

#include <iostream>

Structure loadStructure(const std::string& filename)
{
    Structure structure;

    OpenBabel::OBMol mol;
    OpenBabel::OBConversion conv;

    // detect format automatically
    if(!conv.SetInFormat(conv.FormatFromExt(filename.c_str())))
    {
        std::cerr<<"Unsupported format\n";
        return structure;
    }

    if(!conv.ReadFile(&mol, filename))
    {
        std::cerr<<"Failed to read file\n";
        return structure;
    }

    // iterator required in OpenBabel 3
    OpenBabel::OBAtomIterator ai;

    for(OpenBabel::OBAtom* atom = mol.BeginAtom(ai);
        atom;
        atom = mol.NextAtom(ai))
    {
        AtomSite site;

        site.symbol =
            OpenBabel::OBElements::GetSymbol(atom->GetAtomicNum());

        site.x = atom->GetX();
        site.y = atom->GetY();
        site.z = atom->GetZ();

        structure.atoms.push_back(site);
    }

    return structure;
}