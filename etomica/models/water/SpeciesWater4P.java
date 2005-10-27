package etomica.models.water;
import java.lang.reflect.Constructor;

import etomica.EtomicaElement;
import etomica.atom.AtomTypeGroup;
import etomica.simulation.Simulation;
import etomica.species.Species;
import etomica.species.SpeciesSignature;

/**
 * Species for 4-point water molecule.
 */
public class SpeciesWater4P extends Species implements EtomicaElement {
    
    public SpeciesWater4P(Simulation sim) {
        this(sim, sim.getDefaults().moleculeCount);
    }
    public SpeciesWater4P(Simulation sim, int nM) {
        this(sim, nM, Species.makeAgentType(sim));
    }
    private SpeciesWater4P(Simulation sim, int nM, AtomTypeGroup agentType) {
       super(sim, new AtomFactoryWater4P(sim, agentType),
               agentType);
       nMolecules = nM;
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