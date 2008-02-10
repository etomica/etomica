package etomica.virial;

import etomica.atom.IMolecule;
import etomica.atom.Molecule;
import etomica.chem.elements.Element;
import etomica.chem.elements.ElementSimple;
import etomica.simulation.ISimulation;
import etomica.species.SpeciesSpheresHetero;

public class SpeciesAlkane extends SpeciesSpheresHetero {

    public SpeciesAlkane(ISimulation sim, int numCarbons) {
        super(sim.getSpace(), sim.isDynamic(), makeAtomTypeSpheres(new Element[]{new ElementSimple("CH3", 15), new ElementSimple("CH2", 14)}));
        setTotalChildren(numCarbons);
    }

    public IMolecule makeMolecule() {
        isMutable = false;
        Molecule group = new Molecule(atomType);
        //make straight alkane CH3-CH2-...-CH2-CH3
        group.addChildAtom(makeLeafAtom(leafTypes[0]));
        for(int j = 0; j < childCount[1]; j++) {
            group.addChildAtom(makeLeafAtom(leafTypes[1]));
        }
        if (childCount[0] > 1) {
            group.addChildAtom(makeLeafAtom(leafTypes[0]));
        }
        return group;
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

}
