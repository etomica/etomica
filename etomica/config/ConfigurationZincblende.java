package etomica.config;

import etomica.action.AtomActionTranslateBy;
import etomica.action.AtomGroupAction;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.lattice.LatticeCubicFcc;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;

/**
 * Sets the configuration to the zincblende structure, which consists
 * of two fcc lattices, with one shifted in each direction by one-quarter
 * of the lattice constant.
 */
public class ConfigurationZincblende extends Configuration {
    
    private final ConfigurationLattice fcc;
    private final AtomGroupAction translator;
    private final AtomIteratorListSimple iterator0, iterator1;
    private final LatticeCubicFcc lattice;
    
    public ConfigurationZincblende(Space space) {
        super(space);
        lattice = new LatticeCubicFcc();
        fcc = new ConfigurationLattice(lattice);
        translator = new AtomGroupAction(new AtomActionTranslateBy(space));
        iterator0 = new AtomIteratorListSimple();
        iterator1 = new AtomIteratorListSimple();
    }
    
    /**
     * Initializes positions of atoms to the zincblende structure.  The given
     * array should hold exactly two AtomLists, each with the same number of atoms.
     */
    public void initializePositions(AtomList[] lists) {
        if(lists == null || lists.length != 2) {//need an exception for this
            throw new IllegalArgumentException("inappropriate argument to ConfigurationZincBlende");
        }
        if(lists[0].size() != lists[1].size()) {
            System.err.println("Warning: different numbers of molecules for two species in ConfigurationZincBlende");
        }
        
        //create fcc lattice each species at same positions
        fcc.initializePositions(new AtomList[]{lists[0]});
        fcc.initializePositions(new AtomList[]{lists[1]});
        
        //shift lattice in all three directions by one-quarter the lattice constant
        Vector3D shift = new Vector3D();
        shift.E(0.125*lattice.getLatticeConstant());
        
        ((AtomActionTranslateBy)translator.getAction()).setTranslationVector(shift);

        iterator0.setList(lists[0]);
        iterator0.reset();
        while(iterator0.hasNext()) {
            translator.actionPerformed(iterator0.nextAtom());
        }
        shift.TE(-1.0);
        ((AtomActionTranslateBy)translator.getAction()).setTranslationVector(shift);
        
        iterator1.setList(lists[1]);
        iterator1.reset();
        while(iterator1.hasNext()) {
            translator.actionPerformed(iterator1.nextAtom());
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
        Simulation sim = new Simulation(Space3D.getInstance());
        sim.getDefaults().atomSize = 5.0;
        Space space = sim.space;
        Phase phase = new Phase(sim);
        etomica.species.SpeciesSpheresMono speciesSpheres0  = new etomica.species.SpeciesSpheresMono(sim);
        etomica.species.SpeciesSpheresMono speciesSpheres1  = new etomica.species.SpeciesSpheresMono(sim);
        speciesSpheres0.setNMolecules(32);
        speciesSpheres1.setNMolecules(32);
        etomica.graphics.ColorSchemeByType.setColor(speciesSpheres0,new java.awt.Color(0,255,0));
        etomica.graphics.ColorSchemeByType.setColor(speciesSpheres1, java.awt.Color.red);
        new ConfigurationZincblende(space).initializeCoordinates(phase);

        etomica.graphics.SimulationGraphic simGraphic = new etomica.graphics.SimulationGraphic(sim);
        simGraphic.makeAndDisplayFrame();
    }//end of main
    
}//end of ConfigurationZincBlende