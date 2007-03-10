package etomica.models.water;
import java.lang.reflect.Constructor;

import etomica.EtomicaElement;
import etomica.atom.AtomTypeGroup;
import etomica.simulation.Simulation;
import etomica.species.Species;
import etomica.species.SpeciesSignature;

/**
 * Species for 3-point water molecule.
 */
public class SpeciesWater3P extends Species implements EtomicaElement {
    
    private static final long serialVersionUID = 1L;

    public SpeciesWater3P(Simulation sim) {
        this(sim, Species.makeAgentType(sim));
    }
    private SpeciesWater3P(Simulation sim, AtomTypeGroup agentType) {
       super(new AtomFactoryWater3P(sim, agentType),
               agentType);
    }
    
    public SpeciesSignature getSpeciesSignature() {
        Constructor constructor = null;
        try {
            constructor = this.getClass().getConstructor(new Class[]{Simulation.class});
        }
        catch(NoSuchMethodException e) {
            System.err.println("you have no constructor.  be afraid");
        }
        return new SpeciesSignature(getName(),constructor,new Object[]{});
    }
}