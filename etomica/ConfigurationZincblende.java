package etomica;

import etomica.lattice.LatticeCubicFcc;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;

/**
 * Sets the configuration to the zincblende structure, which consists
 * of two fcc lattices, with one shifted in each direction by one-quarter
 * of the lattice constant.
 */
public class ConfigurationZincblende extends Configuration {
    
    private ConfigurationLattice fcc;
    
    public ConfigurationZincblende() {
        super();
        fcc = new ConfigurationLattice(new LatticeCubicFcc());
    }
    
    /**
     * Initializes positions of atoms to the zincblende structure.  The given
     * array should hold exactly two iterators, each with the same number of iterates.
     */
    public void initializePositions(AtomIterator[] iterators){
        if(iterators == null || iterators.length != 2) {//need an exception for this
            System.err.println("inappropriate argument to ConfigurationZincBlende");
            return;
        }
        if(iterators[0].size() != iterators[1].size()) {
            System.err.println("Warning: different numbers of molecules for two species in ConfigurationZincBlende");
        }
        
        //create fcc lattice each species at same positions
        fcc.initializePositions(iterators[0]);
        fcc.initializePositions(iterators[1]);
        
        //shift lattice in all three directions by one-quarter the lattice constant
        Vector3D shift = new Vector3D();
        shift.E(0.125*fcc.getLatticeConstant());
        
        iterators[0].reset();
        while(iterators[0].hasNext()) {
            iterators[0].nextAtom().coord.translateBy(shift);
        }
        shift.TE(-1.0);
        iterators[1].reset();
        while(iterators[1].hasNext()) {
            iterators[1].nextAtom().coord.translateBy(shift);
        }
    }
    
    public void setDimensions(double[] dimensions) {
        super.setDimensions(dimensions);
        fcc.setDimensions(dimensions);
    }

    /**
     * Displays configuration without setting up full simulation.
     */
    public static void main(String[] args) {
        Simulation sim = new Simulation(Space3D.INSTANCE);
        Default.ATOM_SIZE = 5.0;
        Space space = sim.space;
        Phase phase = new Phase(space);
        etomica.SpeciesSpheresMono speciesSpheres0  = new etomica.SpeciesSpheresMono(sim);
        etomica.SpeciesSpheresMono speciesSpheres1  = new etomica.SpeciesSpheresMono(sim);
        speciesSpheres0.setNMolecules(32);
        speciesSpheres1.setNMolecules(32);
        etomica.graphics.ColorSchemeByType.setColor(speciesSpheres0,new java.awt.Color(0,255,0));
        etomica.graphics.ColorSchemeByType.setColor(speciesSpheres1, java.awt.Color.red);
        phase.setConfiguration(new ConfigurationZincblende());
        phase.speciesMaster.addSpecies(speciesSpheres0);
        phase.speciesMaster.addSpecies(speciesSpheres1);

        etomica.graphics.SimulationGraphic simGraphic = new etomica.graphics.SimulationGraphic(sim);
        simGraphic.makeAndDisplayFrame();
    }//end of main
    
}//end of ConfigurationZincBlende