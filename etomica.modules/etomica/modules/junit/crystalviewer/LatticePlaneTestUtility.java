package etomica.modules.junit.crystalviewer;

import etomica.config.ConfigurationLattice;
import etomica.lattice.BravaisLattice;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.LatticeCubicBcc;
import etomica.lattice.LatticeCubicDiamond;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeCubicSimple;
import etomica.lattice.LatticeHcp;
import etomica.lattice.LatticePlane;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.PrimitiveTetragonal;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.lattice.crystal.PrimitiveOrthorhombic;
import etomica.lattice.crystal.PrimitiveMonoclinic;
import etomica.lattice.crystal.PrimitiveTriclinic;
import etomica.box.Box;
import etomica.simulation.ISimulation;
import etomica.simulation.Simulation;
import etomica.space.IVector;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

public class LatticePlaneTestUtility {

	public final int SIMPLE_CUBIC  = 0;
	public final int TETRAGONAL    = 1;
	public final int HEXAGONAL     = 2;
	public final int ORTHORHOMBIC  = 3;
	public final int MONOCLINIC    = 4;
	public final int TRICLINIC     = 5;
	public final int FCC           = 6;
	public final int BCC           = 7;
	public final int HCP           = 8;
	public final int CUBIC_DIAMOND = 9;

	private ISimulation sim = null;
	private Species species = null;
	private BravaisLattice lattice = null;
	private Box box = null;
	private LatticePlane latticePlane = null;

	public LatticePlaneTestUtility() {
    	// Create a simulation
    	sim = new Simulation(Space3D.getInstance());
    	sim.getSpace().makeVector();

    	// Create a species with one atom / molecule
    	species = new SpeciesSpheresMono(sim);
    	sim.getSpeciesManager().addSpecies(species);

	}

	public void createLatticeAndBox(int index, int[] millerIndices, int[] boxSize) {

		BasisMonatomic basisMonatomic = null;

		// Create the lattice
    	switch(index) {
    	
    	case SIMPLE_CUBIC:
    		lattice = new LatticeCubicSimple();
    		break;
    	case TETRAGONAL:
            basisMonatomic = new BasisMonatomic(sim.getSpace());
            lattice = new BravaisLatticeCrystal(new PrimitiveTetragonal(sim.getSpace()), basisMonatomic);
    		break;
    	case HEXAGONAL:
            basisMonatomic = new BasisMonatomic(sim.getSpace());
            lattice = new BravaisLatticeCrystal(new PrimitiveHexagonal(sim.getSpace()), basisMonatomic);
    		break;
    	case ORTHORHOMBIC:
            basisMonatomic = new BasisMonatomic(sim.getSpace());
            lattice = new BravaisLatticeCrystal(new PrimitiveOrthorhombic(sim.getSpace()), basisMonatomic);
    		break;
    	case MONOCLINIC:
            basisMonatomic = new BasisMonatomic(sim.getSpace());
            lattice = new BravaisLatticeCrystal(new PrimitiveMonoclinic(sim.getSpace()), basisMonatomic);
    		break;
    	case TRICLINIC:
            basisMonatomic = new BasisMonatomic(sim.getSpace());
            lattice = new BravaisLatticeCrystal(new PrimitiveTriclinic(sim.getSpace()), basisMonatomic);
    		break;
    	case FCC:
    		lattice = new LatticeCubicFcc();
    		break;
    	case BCC:
    		lattice = new LatticeCubicBcc();
    		break;
    	case HCP:
    		lattice = new LatticeHcp();
    		break;
    	case CUBIC_DIAMOND:
    		lattice = new LatticeCubicDiamond();
    		break;
    	default:
    		break;
    	}
    		
        // Create the lattice plane
        latticePlane =
        	new LatticePlane(lattice.getPrimitive(), millerIndices);
    	setLatticePlanePosition(0.0);

	    // Create a box
    	if(box != null) {
    	 sim.removeBox(box);
    	}
	    box = new Box(
	    		new etomica.space.BoundaryDeformableLattice(
	                  lattice.getPrimitive(),
	              	  (etomica.util.IRandom)null, boxSize));
	    sim.addBox(box);

	}

	public void setDimensions(int size) {
        // Set the number of atoms
        int numAtoms = size*size*size;

	    // Set the dimensions for the box
	    IVector dimensions = box.getBoundary().getDimensions();
	    dimensions.E(lattice.getPrimitive().getSize());
	    dimensions.TE(size);
	    box.setDimensions(dimensions);

	    // Set the number of molecules for the box and
	    // initialze the positions.
	    box.getAgent(species).setNMolecules(numAtoms);
	    ConfigurationLattice config = new ConfigurationLattice(lattice);
	    config.initializeCoordinates(box);
		
	}

	public void setLatticePlanePosition(double pos) {
		latticePlane.setPosition(pos);
	}

	public double getLatticePlaneSpacePosition() {
		return latticePlane.getSpacePosition();
	}

	public ISimulation getSimulation() {
		return sim;
	}
	
	public Species getSpecies() {
		return species;
	}
	
	public BravaisLattice getLattice() {
		return lattice;
	}
	
	public Box getBox() {
		return box;
	}
	
	public LatticePlane getLatticePlane() {
		return latticePlane;
	}
}
