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
    
    public SpeciesWater3P(Simulation sim) {
        this(sim, sim.getDefaults().moleculeCount);
    }
    public SpeciesWater3P(Simulation sim, int nM) {
        this(sim, nM, Species.makeAgentType(sim));
    }
    private SpeciesWater3P(Simulation sim, int nM, AtomTypeGroup agentType) {
       super(sim, new AtomFactoryWater3P(sim, agentType),
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