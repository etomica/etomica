package etomica.models.water;
import etomica.*;

public class SpeciesWater extends Species implements EtomicaElement {
    
   public SpeciesWater() {
        this(Simulation.instance);
    }     
  
    public SpeciesWater(Simulation sim) {
        this(sim, Default.MOLECULE_COUNT);
    }
    public SpeciesWater(Simulation sim, int nM) {
       super(sim, new AtomFactoryWater(sim));
       factory.setSpecies(this);
       nMolecules = nM;
    }
    
}