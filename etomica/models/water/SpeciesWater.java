package etomica.models.water;
import java.lang.reflect.Constructor;

import etomica.EtomicaElement;
import etomica.Simulation;
import etomica.atom.AtomTypeGroup;
import etomica.species.Species;
import etomica.species.SpeciesSignature;
import etomica.util.Default;

public class SpeciesWater extends Species implements EtomicaElement {
    
    public SpeciesWater(Simulation sim) {
        this(sim, Default.MOLECULE_COUNT);
    }
    public SpeciesWater(Simulation sim, int nM) {
        this(sim, nM, Species.makeAgentType(sim));
    }
    private SpeciesWater(Simulation sim, int nM, AtomTypeGroup agentType) {
       super(sim, new AtomFactoryWater(sim, agentType),
               agentType);
       factory.setSpecies(this);
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