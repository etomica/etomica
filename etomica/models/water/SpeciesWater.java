package etomica.models.water;
import etomica.AtomType;
import etomica.Default;
import etomica.EtomicaElement;
import etomica.Simulation;
import etomica.Species;

public class SpeciesWater extends Species implements EtomicaElement {
    
    public SpeciesWater(Simulation sim) {
        this(sim, Default.MOLECULE_COUNT);
    }
    public SpeciesWater(Simulation sim, int nM) {
        this(sim, nM, Species.makeAgentType(sim));
    }
    private SpeciesWater(Simulation sim, int nM, AtomType agentType) {
       super(sim, new AtomFactoryWater(sim.space, agentType.getIndexManager().makeChildManager()),
               agentType);
       factory.setSpecies(this);
       nMolecules = nM;
    }
    
}