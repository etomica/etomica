package etomica;

/**
 * Sets the configuration to the zincblende structure, which consists
 * of two fcc lattices, with one shifted in each direction by one-quarter
 * of the lattice constant.
 */
public class ConfigurationZincblende extends Configuration {
    
    private ConfigurationFcc fcc;
    
    public ConfigurationZincblende(Space space) {
        super(space);
        if(space.D() != 3) {//need an exception for this
            throw new IllegalArgumentException("Illegal dimension for ConfigurationZincBlende");
        }
        fcc = new ConfigurationFcc(space);
    }
    
    /**
     * Initializes positions of atoms to the zincblende structure.  The given
     * array should hold exactly two iterators, each with the same number of iterates.
     */
    public void initializePositions(AtomIterator[] iterators){
        if(iterators == null || iterators.length != 2) {//need an exception for this
            System.out.println("inappropriate argument to ConfigurationZincBlende");
            return;
        }
        if(iterators[0].size() != iterators[1].size()) {
            System.out.println("Warning: different numbers of molecules for two species in ConfigurationZincBlende");
        }
        
        //create fcc lattice each species at same positions
        fcc.initializePositions(iterators[0]);
        fcc.initializePositions(iterators[1]);
        
        //shift lattice in all three directions by one-quarter the lattice constant
        Space3D.Vector shift = new Space3D.Vector();
        shift.E(0.25*fcc.latticeConstant());
        
        iterators[1].reset();
        while(iterators[1].hasNext()) {
            Atom atom = iterators[1].nextAtom();
            atom.coord.translateBy(shift);
        }
    }

    /**
     * Displays configuration without setting up full simulation.
     */
//    public static void main(String[] args) {
//        etomica.graphics.SimulationGraphic sim = new etomica.graphics.SimulationGraphic(new etomica.Space3D());
//        Simulation.instance = sim;
//        Default.ATOM_SIZE = 5.0;
//    //    Default.DISPLAY_USE_OPENGL = false;
//        etomica.Phase phase0  = new etomica.Phase(sim.space);
//        phase0.setConfiguration(new ConfigurationZincblende(sim.space));
//        etomica.SpeciesSpheresMono speciesSpheres0  = new etomica.SpeciesSpheresMono(sim);
//        etomica.SpeciesSpheresMono speciesSpheres1  = new etomica.SpeciesSpheresMono(sim);
//        speciesSpheres0.setNMolecules(32);
//        speciesSpheres1.setNMolecules(32);
//        etomica.graphics.ColorSchemeByType.setColor(speciesSpheres0,new java.awt.Color(0,255,0));
//        etomica.graphics.ColorSchemeByType.setColor(speciesSpheres1, java.awt.Color.red);
//        etomica.graphics.DisplayPhase displayPhase0  = new etomica.graphics.DisplayPhase(phase0);
//        sim.mediator().go(); 
//        etomica.graphics.SimulationGraphic.makeAndDisplayFrame(sim);
//    }//end of main
    
}//end of ConfigurationZincBlende