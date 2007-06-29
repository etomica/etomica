package etomica.config;

import etomica.action.AtomActionTranslateBy;
import etomica.action.AtomGroupAction;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayPhase;
import etomica.lattice.LatticeCubicFcc;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;

/**
 * Sets the configuration to the zincblende structure, which consists
 * of two fcc lattices, with one shifted in each direction by one-quarter
 * of the lattice constant.
 */
public class ConfigurationZincblende extends ConfigurationLattice {
    
    private static final long serialVersionUID = 2L;
    private AtomGroupAction translator0, translator1;
    
    public ConfigurationZincblende(double latticeConstant) {
        super(new LatticeCubicFcc(latticeConstant));
    }
    
    
    /**
     * Initializes positions of atoms to the zincblende structure.  The given
     * array should hold exactly two AtomLists, each with the same number of atoms.
     */
    public void initializeCoordinates(Phase phase) {
        translator0 = new AtomGroupAction(new AtomActionTranslateBy(phase.getSpace()));
        translator1 = new AtomGroupAction(new AtomActionTranslateBy(phase.getSpace()));
        AtomSet[] lists = getMoleculeLists(phase);
        if(lists == null || lists.length != 2) {//need an exception for this
            throw new IllegalArgumentException("inappropriate argument to ConfigurationZincBlende");
        }
        if(lists[0].getAtomCount() != lists[1].getAtomCount()) {
            System.err.println("Warning: different numbers of molecules for two species in ConfigurationZincBlende");
        }
        
        int nCells = (int) Math.ceil(lists[0].getAtomCount() / 4.0);

        // determine scaled shape of simulation volume
        IVector shape = phase.getSpace().makeVector();
        shape.E(phase.getBoundary().getDimensions());
        IVector latticeConstantV = Space.makeVector(lattice.getLatticeConstants());
        shape.DE(latticeConstantV);

        // determine number of cells in each direction
        int[] latticeDimensions = calculateLatticeDimensions(nCells, shape);
        if (indexIterator.getD() > latticeDimensions.length) {
            int[] iteratorDimensions = new int[latticeDimensions.length+1];
            System.arraycopy(latticeDimensions, 0, iteratorDimensions, 0,
                    latticeDimensions.length);
            iteratorDimensions[latticeDimensions.length] = 4;
            indexIterator.setSize(iteratorDimensions);
        }
        else {
            indexIterator.setSize(latticeDimensions);
        }

        //shift lattice in all three directions by one-quarter the lattice constant
        Vector3D shift = new Vector3D();
        shift.Ea1Tv1(-0.5,phase.getBoundary().getDimensions());
        shift.PE(0.125*((LatticeCubicFcc)lattice).getLatticeConstant());
        ((AtomActionTranslateBy)translator0.getAction()).setTranslationVector(shift);

        shift.PE(0.25*((LatticeCubicFcc)lattice).getLatticeConstant());
        ((AtomActionTranslateBy)translator1.getAction()).setTranslationVector(shift);

        // Place molecules
        indexIterator.reset();
        int i = 0;
        while (indexIterator.hasNext()) {
            int[] ii = indexIterator.next();
            IVector site = (IVector) lattice.site(ii);
            atomActionTranslateTo.setDestination(site);

            IAtom a0 = lists[0].getAtom(i);
            IAtom a1 = lists[1].getAtom(i);
            i++;
            atomActionTranslateTo.actionPerformed(a0);
            atomActionTranslateTo.actionPerformed(a1);

            translator0.actionPerformed(a0);
            translator1.actionPerformed(a1);
        }
    }        
    
    /**
     * Displays configuration without setting up full simulation.
     */
    public static void main(String[] args) {
    	final String APP_NAME = "Configuration Zinc Blende";

        Simulation sim = new Simulation(Space3D.getInstance());
        final Phase phase = new Phase(sim);
        phase.getBoundary().setDimensions(new etomica.space3d.Vector3D(30.0, 30.0, 30.0));
        sim.addPhase(phase);
        etomica.species.SpeciesSpheresMono speciesSpheres0  = new etomica.species.SpeciesSpheresMono(sim);
        etomica.species.SpeciesSpheresMono speciesSpheres1  = new etomica.species.SpeciesSpheresMono(sim);
        sim.getSpeciesManager().addSpecies(speciesSpheres0);
        sim.getSpeciesManager().addSpecies(speciesSpheres1);
        ((AtomTypeSphere)speciesSpheres0.getMoleculeType()).setDiameter(5.0);
        ((AtomTypeSphere)speciesSpheres1.getMoleculeType()).setDiameter(5.0);
        phase.getAgent(speciesSpheres0).setNMolecules(32);
        phase.getAgent(speciesSpheres1).setNMolecules(32);
        ConfigurationZincblende config = new ConfigurationZincblende(15);
        config.initializeCoordinates(phase);

        final etomica.graphics.SimulationGraphic simGraphic = new etomica.graphics.SimulationGraphic(sim, APP_NAME);
        simGraphic.add(new DisplayPhase(phase));
        ColorSchemeByType colorScheme = (ColorSchemeByType)simGraphic.getDisplayPhase(phase).getColorScheme();
        colorScheme.setColor(speciesSpheres0.getMoleculeType(),new java.awt.Color(0,255,0));
        colorScheme.setColor(speciesSpheres1.getMoleculeType(), java.awt.Color.red);

        simGraphic.getController().getSimRestart().setConfiguration(config);
        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getDisplayPhasePaintAction(phase));

        simGraphic.makeAndDisplayFrame(APP_NAME);
        simGraphic.getDisplayPhase(phase).graphic().repaint();
    }
    
}
