package etomica.models.water;
import etomica.Default;
import etomica.EtomicaElement;
import etomica.Simulation;
import etomica.Species;

public class SpeciesWater extends Species implements EtomicaElement {
    
    public SpeciesWater(Simulation sim) {
        this(sim, Default.MOLECULE_COUNT);
    }
    public SpeciesWater(Simulation sim, int nM) {
       super(sim, new AtomFactoryWater(sim.space, sim.speciesRoot.type.getIndexManager().makeMoleculeIndexManager()));
       factory.setSpecies(this);
       nMolecules = nM;
    }
    
}