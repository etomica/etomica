package etomica;
import etomica.lattice.*;

public class ConfigurationFcc extends Configuration {
    
    LatticeCubicFcc fcc;
    
    public ConfigurationFcc(Simulation sim) {
        super(sim);
    }
    
    public double latticeConstant() {
        if(fcc != null) return fcc.getLatticeConstant();
        else return Double.NaN;
    }
    
    public void initializePositions(AtomIterator[] iterators){
        if(iterators == null || iterators.length == 0) return;
        AtomIterator iterator;
        if(iterators.length == 1) iterator = iterators[0];
        else iterator = new AtomIteratorCompound(iterators);//lump 'em all together
        if(iterator.size() == 0) {return;}
        if(iterator.size() == 1) {
            iterator.reset();
            iterator.nextAtom().coord.translateTo(space.origin());
            return;
        }
    
    // Count number of molecules
        int sumOfMolecules = iterator.size();
        if(sumOfMolecules == 0) {return;}
        
        fcc = LatticeCubicFcc.makeLattice(space, sumOfMolecules, dimensions[0]);
        fcc.shiftFirstToOrigin();

   // Place molecules  
        Space.Vector[] rLat = fcc.positions();
        int i = 0;
        iterator.reset();
        while(iterator.hasNext()) {
            Atom a = iterator.nextAtom();
            //initialize coordinates of child atoms
            try {//may get null pointer exception when beginning simulation
                a.creator().getConfiguration().initializeCoordinates(a);
            } catch(NullPointerException e) {}
            a.coord.translateTo(rLat[i++]);//use translateTo instead of E because atom might be a group
        }
    }//end of initializePositions
    
    
    //used only by ConfigurationSequential
    /**
     * @deprecated
     */
    public static Space3D.Vector[] lattice(int n) { 
        Space3D.Vector[] r = new Space3D.Vector[n];
        for(int i=0; i<n; i++) {r[i] = new Space3D.Vector();}
        LatticeCubicFcc fcc = LatticeCubicFcc.makeLattice(Simulation.instance.space, n, Default.BOX_SIZE);//need to extend--assumes cubic box
        AtomIteratorList iteratorsites = new AtomIteratorList(fcc.siteList());
        iteratorsites.reset();
        int i = 0;
        while (iteratorsites.hasNext()&& i < n){
            Atom site = iteratorsites.nextAtom();
            r[i].E(site.coord.position());
            i++ ;
        }
        return r;
    }//end of lattice
    
  public static void main(String[] args) {
    etomica.graphics.SimulationGraphic sim = new etomica.graphics.SimulationGraphic(new etomica.Space3D());
    Default.ATOM_SIZE = 6.6;
//    Default.DISPLAY_USE_OPENGL = false;
    etomica.Phase phase0  = new etomica.Phase();
    phase0.setConfiguration(new ConfigurationFcc(sim));
    etomica.SpeciesSpheresMono speciesSpheres0  = new etomica.SpeciesSpheresMono();
//    etomica.SpeciesSpheres speciesSpheres0  = new etomica.SpeciesSpheres();
      speciesSpheres0.setNMolecules(32);
    etomica.graphics.DisplayPhase displayPhase0  = new etomica.graphics.DisplayPhase();
    displayPhase0.setColorScheme(new etomica.graphics.ColorSchemeByType());
    etomica.graphics.ColorSchemeByType.setColor(speciesSpheres0, new java.awt.Color(0,255,0));
    sim.mediator().go(); 
    etomica.graphics.SimulationGraphic.makeAndDisplayFrame(sim);
  }//end of main
  
}//end of ConfigurationFcc