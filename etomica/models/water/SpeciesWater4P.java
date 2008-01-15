package etomica.models.water;
import java.lang.reflect.Constructor;

import etomica.simulation.ISimulation;
import etomica.species.Species;
import etomica.species.SpeciesSignature;

/**
 * Species for 4-point water molecule.
 */
public class SpeciesWater4P extends Species {
    
    private static final long serialVersionUID = 1L;

    public SpeciesWater4P(ISimulation sim) {
       super();
       setMoleculeFactory(new AtomFactoryWater4P(sim, this));
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
}