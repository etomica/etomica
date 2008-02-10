package etomica.species;

import java.lang.reflect.Constructor;

import etomica.atom.AtomLeafAngular;
import etomica.atom.AtomLeafAngularDynamic;
import etomica.atom.AtomTypeOrientedSphere;
import etomica.atom.IAtomLeaf;
import etomica.chem.elements.ElementSimple;
import etomica.simulation.ISimulation;

/**
 * Species in which molecules are made of a single atom of type OrientedSphere
 *
 * @author David Kofke
 * @see AtomTypeOrientedSphere
 * 
 */
public class SpeciesSpheresRotating extends SpeciesSpheresMono {
    
    public SpeciesSpheresRotating(ISimulation sim) {
        super(sim.getSpace(), sim.isDynamic(), new AtomTypeOrientedSphere(new ElementSimple(sim), 1.0));
    }

    protected IAtomLeaf makeLeafAtom() {
        return isDynamic ? new AtomLeafAngularDynamic(space, leafAtomType)
                         : new AtomLeafAngular(space, leafAtomType);
    }
    
              
    public SpeciesSignature getSpeciesSignature() {
        Constructor constructor = null;
        try {
            constructor = this.getClass().getConstructor(new Class[]{ISimulation.class});
        }
        catch(NoSuchMethodException e) {
            System.err.println("you have no constructor.  be afraid");
        }
        return new SpeciesSignature(constructor,new Object[]{});
    }

    private static final long serialVersionUID = 1L;
}
