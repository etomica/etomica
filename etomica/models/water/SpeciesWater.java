package etomica.models.water;
import etomica.Default;
import etomica.EtomicaElement;
import etomica.Simulation;
import etomica.Space;
import etomica.Species;

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