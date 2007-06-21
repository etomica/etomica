package etomica.models.propane;
import java.lang.reflect.Constructor;

import etomica.simulation.ISimulation;
import etomica.species.Species;
import etomica.species.SpeciesSignature;

/**
 * Species for 3-point water molecule.
 */
public class SpeciesUAPropane extends Species {
    
    private static final long serialVersionUID = 1L;

    public SpeciesUAPropane(ISimulation sim) {
       super(new AtomFactoryUAPropane(sim));
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