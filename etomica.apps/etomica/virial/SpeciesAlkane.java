package etomica.virial;

import etomica.api.IAtomTypeLeaf;
import etomica.api.IMolecule;
import etomica.api.ISimulation;
import etomica.atom.Molecule;
import etomica.chem.elements.Element;
import etomica.chem.elements.ElementSimple;
import etomica.space.Space;
import etomica.species.SpeciesSpheresHetero;

public class SpeciesAlkane extends SpeciesSpheresHetero {

    public SpeciesAlkane(ISimulation sim, Space _space, int numCarbons) {
        super(_space, sim.isDynamic(), makeAtomTypeSpheres(new Element[]{new ElementSimple("CH3", 15), new ElementSimple("CH2", 14)}));
        setTotalChildren(numCarbons);
    }

    public IMolecule makeMolecule() {
        Molecule group = new Molecule(this);
        //make straight alkane CH3-CH2-...-CH2-CH3
        group.addChildAtom(makeLeafAtom(childTypes[0]));
        for(int j = 0; j < childCount[1]; j++) {
            group.addChildAtom(makeLeafAtom(childTypes[1]));
        }
        if (childCount[0] > 1) {
            group.addChildAtom(makeLeafAtom(childTypes[0]));
        }
        conformation.initializePositions(group.getChildList());
        return group;
    }
    
    public IAtomTypeLeaf getCH2Type() {
        return childTypes[0];
    }
    
    public IAtomTypeLeaf getCH3Type() {
        return childTypes[1];
    }
    
    public void setTotalChildren(int newTotalChildren) {
        if (newTotalChildren > 1) {
            childCount[0] = 2;
            childCount[1] = newTotalChildren - 2;
        }
        else {
            childCount[0] = 1;
            childCount[1] = 0;
        }
    }

    private static final long serialVersionUID = 1L;
}
