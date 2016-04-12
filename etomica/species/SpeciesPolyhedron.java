package etomica.species;

import java.util.List;

import etomica.api.IAtom;
import etomica.api.IElement;
import etomica.api.IVector;
import etomica.atom.AtomOrientedQuaternion;
import etomica.atom.AtomTypeOrientedSphere;
import etomica.atom.AtomTypeSpheroPolyhedron;
import etomica.chem.elements.ElementSimple;
import etomica.simulation.Simulation;
import etomica.space.ISpace;

/**
 * Species in which molecules are made of a single atom of type SpheroPolyhedron
 *
 * @author Andrew Schultz
 * @see AtomTypeOrientedSphere
 */
public class SpeciesPolyhedron extends SpeciesSpheresMono {
    
    public SpeciesPolyhedron(Simulation sim, ISpace _space, List<IVector> vertices, double sweepRadius) {
        this(_space, vertices, sweepRadius, new ElementSimple(sim));
    }
    
    public SpeciesPolyhedron(ISpace _space, List<IVector> vertices, double sweepRadius, IElement element) {
        super(_space, new AtomTypeSpheroPolyhedron(element, _space, vertices, sweepRadius));
    }

    protected IAtom makeLeafAtom() {
        return new AtomOrientedQuaternion(space, leafAtomType);
    }
}
