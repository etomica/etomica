package etomica.models.water;
import etomica.graphics.*;
import etomica.*;
import java.awt.Color;

public class SpeciesWater extends Species implements EtomicaElement {
    
    public AtomType.Sphere[] protoType;
    public SpeciesAgent sa;
    private double mass;
    
    private static AtomFactoryHetero makeFactory(Simulation sim) {
        int i = ColorSchemeByType.colorIndex;
        AtomFactoryMono hFactory = new AtomFactoryMono(sim, sim.getIteratorFactory().simpleSequencerFactory());
        AtomFactoryMono oFactory = new AtomFactoryMono(sim, sim.getIteratorFactory().simpleSequencerFactory());
        AtomType hType = new AtomType.Sphere(hFactory, 1.0, /*Electron.UNIT.toSim(0.4238),*/ 2.0/2);
        AtomType oType = new AtomType.Sphere(oFactory, 16.0, /*Electron.UNIT.toSim(-0.8476),*/ 3.167/2);
        
        hFactory.setType(hType);
        oFactory.setType(oType);
        ColorSchemeByType.setColor(hType, Color.green);
        ColorSchemeByType.setColor(oType, Color.pink);
        
        AtomFactory[] childFactory = new AtomFactory[3];
        childFactory[0] = oFactory;
        childFactory[1] = hFactory;
        childFactory[2] = hFactory;
        ConfigurationWater config = new ConfigurationWater(sim); 
        AtomSequencer.Factory sequencerFactory = sim.getIteratorFactory().neighborSequencerFactory();
        AtomFactoryHetero f = new AtomFactoryHetero(sim,sequencerFactory,childFactory,config);
        return f;
    }   
   
   public SpeciesWater() {
        this(Simulation.instance);
    }     
  
    public SpeciesWater(Simulation sim) {
        this(sim, Default.MOLECULE_COUNT);
    }
    public SpeciesWater(Simulation sim, int nM) {
       super(sim, makeFactory(sim));
       factory.setSpecies(this);
       nMolecules = nM;
    }
    
    
}