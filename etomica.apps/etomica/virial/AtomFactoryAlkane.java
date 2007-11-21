package etomica.virial;

import etomica.atom.AtomFactoryHetero;
import etomica.atom.AtomGroup;
import etomica.atom.IAtom;
import etomica.config.Conformation;
import etomica.space.Space;

public class AtomFactoryAlkane extends AtomFactoryHetero {
    public AtomFactoryAlkane(Space space, Conformation conformation) {
        super(space, conformation);
    }
    
    public IAtom makeAtom() {
        isMutable = false;
        AtomGroup group = new AtomGroup(atomType);
        //make straight alkane CH3-CH2-...-CH2-CH3
        group.addChildAtom(childFactory[0].makeAtom());
        for(int j = 0; j < childCount[1]; j++) {
            group.addChildAtom(childFactory[1].makeAtom());
        }
        if (childCount[0] > 1) {
            group.addChildAtom(childFactory[0].makeAtom());
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
