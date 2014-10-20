package etomica.species;

import etomica.api.IAtom;
import etomica.api.IElement;
import etomica.atom.AtomOriented;
import etomica.atom.AtomOrientedDynamic;
import etomica.atom.AtomTypeOrientedSphere;
import etomica.chem.elements.ElementSimple;
import etomica.simulation.Simulation;
import etomica.space.ISpace;

/**
 * Species in which molecules are made of a single atom of type OrientedSphere
 *
 * @author David Kofke
 * @see AtomTypeOrientedSphere
 */
public class SpeciesSpheresRotating extends SpeciesSpheresMono {
    
    public SpeciesSpheresRotating(Simulation sim, ISpace _space) {
        this(_space, new ElementSimple(sim));
    }
    
    public SpeciesSpheresRotating(ISpace _space, IElement element) {
        super(_space, new AtomTypeOrientedSphere(element, _space));
    }
    
    public void setAxisSymmetric(boolean isAxisSymmetric) {
        this.isAxisSymmetric = isAxisSymmetric;
    }

    protected IAtom makeLeafAtom() {
        return isDynamic ? new AtomOrientedDynamic(space, leafAtomType, isAxisSymmetric)
                         : new AtomOriented(space, leafAtomType, isAxisSymmetric);
    }

    protected boolean isAxisSymmetric;
}
