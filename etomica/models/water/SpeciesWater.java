package etomica.models.water;
import etomica.*;

public class SpeciesWater extends Species implements EtomicaElement {
    
   public SpeciesWater() {
        this(Simulation.instance.space);
    }     
  
    public SpeciesWater(Space space) {
        this(space, Default.MOLECULE_COUNT);
    }
    public SpeciesWater(Space space, int nM) {
       super(new AtomFactoryWater(space));
       factory.setSpecies(this);
       nMolecules = nM;
    }
    
}