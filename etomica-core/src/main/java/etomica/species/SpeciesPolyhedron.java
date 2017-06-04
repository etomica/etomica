package etomica.species;

import java.util.List;

import etomica.atom.IAtom;
import etomica.api.IElement;
import etomica.space.Vector;
import etomica.atom.AtomOrientedQuaternion;
import etomica.atom.AtomTypeOrientedSphere;
import etomica.atom.AtomTypeSpheroPolyhedron;
import etomica.chem.elements.ElementSimple;
import etomica.simulation.Simulation;
import etomica.space.Space;

/**
 * Species in which molecules are made of a single atom of type SpheroPolyhedron
 *
 * @author Andrew Schultz
 * @see AtomTypeOrientedSphere
 */
public class SpeciesPolyhedron extends SpeciesSpheresMono {
    
    public SpeciesPolyhedron(Simulation sim, Space _space, List<Vector> vertices, double sweepRadius) {
        this(_space, vertices, sweepRadius, new ElementSimple(sim));
    }
    
    public SpeciesPolyhedron(Space _space, List<Vector> vertices, double sweepRadius, IElement element) {
        super(_space, new AtomTypeSpheroPolyhedron(element, _space, vertices, sweepRadius));
    }

    protected IAtom makeLeafAtom() {
        return new AtomOrientedQuaternion(space, leafAtomType);
    }
}
