package etomica.species;

import etomica.api.IAtom;
import etomica.api.ISimulation;
import etomica.atom.AtomOriented;
import etomica.atom.AtomOrientedDynamic;
import etomica.atom.AtomTypeOrientedSphere;
import etomica.chem.elements.ElementSimple;
import etomica.space.ISpace;

/**
 * Species in which molecules are made of a single atom of type OrientedSphere
 *
 * @author David Kofke
 * @see AtomTypeOrientedSphere
 */
public class SpeciesSpheresRotating extends SpeciesSpheresMono {
    
    public SpeciesSpheresRotating(ISimulation sim, ISpace _space) {
        super(_space, new AtomTypeOrientedSphere(new ElementSimple(sim), 1.0, _space));
    }

    protected IAtom makeLeafAtom() {
        return isDynamic ? new AtomOrientedDynamic(space, leafAtomType)
                         : new AtomOriented(space, leafAtomType);
    }
    
    private static final long serialVersionUID = 1L;
}
