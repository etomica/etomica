package etomica.models.water;
import java.lang.reflect.Constructor;

import etomica.EtomicaElement;
import etomica.simulation.ISimulation;
import etomica.species.Species;
import etomica.species.SpeciesSignature;

/**
 * Species for 4-point water molecule.
 */
public class SpeciesWater4P extends Species implements EtomicaElement {
    
    private static final long serialVersionUID = 1L;

    public SpeciesWater4P(ISimulation sim) {
       super(new AtomFactoryWater4P(sim));
    }
    
    public SpeciesSignature getSpeciesSignature() {
        Constructor constructor = null;
        try {
            constructor = this.getClass().getConstructor(new Class[]{ISimulation.class});
        }
        catch(NoSuchMethodException e) {
            System.err.println("you have no constructor.  be afraid");
        }
        return new SpeciesSignature(getName(),constructor,new Object[]{});
    }
}