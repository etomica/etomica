package etomica.config;

import java.util.ArrayList;
import java.util.HashMap;

import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.integrator.IntegratorHard;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.IndexIteratorSizable;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.SpaceLattice;
import etomica.math.geometry.Plane;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

/**
 * Constructs configuration that has the molecules placed on the sites of a
 * lattice. Molecules are placed sequentially on the sites of the lattice in an
 * order specified by an index iterator set at construction.
 * <p>
 * Iteration over the indexes yields integer arrays, and with each iteration the
 * array is passed to the site method of the lattice which returns the position
 * for placement of the next molecule. Each array index is iterated to a maximum
 * value determined by the number of molecules to be placed, the dimensions of
 * the phase in which they are placed, and the lattice constants of the lattice.
 * <p>
 * An instance of this class may be configured to place atoms such that they
 * uniformly fill the volume of the phase. It will attempt this by scaling the
 * lattice constants of the configuration in an appropriate way. Success in
 * getting a good spatial distribution may vary.
 * <p>
 * An instance can also be configured to remember the indices used to get
 * lattice position for each molecule that is placed. This can be useful if it
 * is desired to associate each molecule with a lattice site.
 */
public class ConfigurationLatticeWithPlane extends ConfigurationLattice {

    private static final long serialVersionUID = 2L;
    private Plane plane;
    private ArrayList species;
    private HashMap allocation;
    
    private final int LEFT  = 0;
    private final int RIGHT = 1;

    /**
     * Constructs class using instance of IndexIteratorRectangular as the default
     * index iterator.
     */
    public ConfigurationLatticeWithPlane(SpaceLattice lattice, Plane plane) {
        this(lattice, plane, new IndexIteratorRectangular(lattice.D()));
    }

    /**
     * Construct class that will place atoms on sites of the given lattice,
     * proceeding in the order resulting from iteration through the given index
     * iterator.
     */
    private ConfigurationLatticeWithPlane(SpaceLattice lattice,
            Plane plane, IndexIteratorSizable indexIterator) {
    	super(lattice, indexIterator);

    	
        if(indexIterator.getD() != lattice.D()) {
            throw new IllegalArgumentException("Dimension of index iterator and lattice are incompatible");
        }
        this.plane = plane;
        species = new ArrayList();
        allocation = new HashMap();
    }

    /**
     * 
     * @param newSpecies : species that will need to be configured
     */
    public void addSpecies(Species newSpecies) {
    	if(species.indexOf(newSpecies) == -1) {
    	    species.add(newSpecies);
    	    allocation.put(newSpecies, new Float(0.5));
    	}
    }

    /**
     * 
     * @param newSpecies : species that should be removed from configuration
     */
    public void removeSpecies(Species remSpecies) {
    	if(species.indexOf(remSpecies) != -1) {
    		species.remove(remSpecies);
    	    allocation.remove(remSpecies);
    	}
    }

    /**
     * Sets the initial allocation of the species on the left side
     * of the plane.  The remaining molecules in the species will
     * be initially allocated to the right side of the plane.
     * 
     * @param sp : Species to set allocation for
     * @param pct : Percentage of species to initially allocate
     * on the left side of the plane (0.0 <= pct <= 1.0)
     */
    public void setSpeciesAllocation(Species sp, float pct) {
    	if(species.indexOf(sp) != -1 && pct <= 1.0) {
    		allocation.remove(sp);
    		allocation.put(sp, new Float(pct));
    	}
    }

    /**
     * Sets the initial allocation of the species on the left side
     * of the plane.  The remaining molecules in the species will
     * be initially allocated to the right side of the plane.
     * 
     * @param sp : Species to set allocation for
     * @param pct : Percentage of species to initially allocate
     * on the left side of the plane (0.0 <= pct <= 1.0)
     */
    public float getSpeciesAllocation(Species sp) {
    	float fValue = 0.0f;
    	if(allocation.containsKey(sp)) {
    		Float value = ((Float)allocation.get(sp));
    		fValue = value.floatValue();
    	}
		return fValue;
    }

    /**
     * Places the molecules in the given phase on the positions of the
     * lattice.  The molecules must be of a species added with addSpecies()
     * or they will not be displayed.
     */
// Method is assuming plane is in the middle.  Calculations making this
// assumption are noted.
// Method also currently hardcoded for 3D
    public void initializeCoordinates(Phase phase) {
//System.out.println("---------------------------------------");
//System.out.println("mc count = " + phase.getSpeciesMaster().moleculeCount());

        int numSpecies = species.size();
        int[] speciesCount = new int[numSpecies];
//System.out.println("# species : " + numSpecies);

        int sumOfMolecules = 0;
        for (int i = 0; i < numSpecies; i++) {
//System.out.println(" sp = " + ((Species)species.get(i)).getAgent(phase).getNMolecules());
            speciesCount[i] = ((Species)species.get(i)).getAgent(phase).getNMolecules();
            sumOfMolecules = sumOfMolecules + speciesCount[i];
        }

        if (sumOfMolecules == 0) {
//System.out.println("Should not be showing any molecules");
            return;
        }
//System.out.println("sumOfMolecules : " + sumOfMolecules);

        // Allocate species to left or right side of plane
// THE FOLLOWING WILL NEED TO BE CHANGED :
// for now, going to allocate 50/50 for each species.
        int [][] molecules = new int[RIGHT+1][numSpecies];
        int moleculesOnLeft = 0;
        int moleculesOnRight = 0;

        for (int i = 0; i < numSpecies; i++) {
        	Species sp = ((Species)species.get(i));
// NEXT LINE NEEDS TO DETERMINE ACTUAL PCT ON LEFT AND NOT JUST DIVIDE BY 2
            molecules[LEFT][i] =  (int)(sp.getAgent(phase).getNMolecules() * getSpeciesAllocation(sp));
            molecules[RIGHT][i] = sp.getAgent(phase).getNMolecules() - molecules[LEFT][i];
            moleculesOnLeft += molecules[LEFT][i];
            moleculesOnRight += molecules[RIGHT][i];
//System.out.println("  on left : " + molecules[LEFT][i]);
//System.out.println("  on right : " + molecules[RIGHT][i]);
        }

        // determine scaled shape of simulation volume
        IVector shape = phase.getSpace().makeVector();
        shape.E(phase.getBoundary().getDimensions());

//System.out.println("shapeX(pre-plane adjustment) = " + shape.x(0));
//System.out.println("shapeY(pre-plane adjustment) = " + shape.x(1));
//System.out.println("shapeZ(pre-plane adjustment) = " + shape.x(2));

//System.out.println("plane a = " + plane.getA());
//System.out.println("plane b = " + plane.getB());
//System.out.println("plane c = " + plane.getC());

    int planeDimIdx = 0;
    if(plane.getA() > plane.epsilon) {
        planeDimIdx = 0;
    }
    else if(plane.getB() > plane.epsilon) {
    	planeDimIdx = 1;
    }
    else if(phase.getSpace().D() == 3 &&
    		plane.getC() > plane.epsilon) {
    	planeDimIdx = 2;
    }

    // Need to adjust shape on either side of plane.
    //  NOTE, JUST DIVIDING BY 2 ASSUMES PLANE DOWN CENTER OF PHASE
    shape.setX(0, shape.x(planeDimIdx) / 2);

//System.out.println("shapeX(post-plane adjustment) = " + shape.x(0));
//System.out.println("shapeY(post-plane adjustment) = " + shape.x(1));
//System.out.println("shapeZ(post-plane adjustment) = " + shape.x(2));

        IVector latticeConstantV = Space.makeVector(lattice.getLatticeConstants());
        shape.DE(latticeConstantV);

        int maxMolecules = Math.max(moleculesOnLeft, moleculesOnRight);
        // determine number of cells in each direction
        int[] latticeDimensions = calculateLatticeDimensions((maxMolecules*2), shape);

//System.out.println("indexIterator.getD() = " + indexIterator.getD());
//System.out.println("latticeDimensions.length = " + latticeDimensions.length);
//System.out.println("  latticeDimensions[0] = " + latticeDimensions[0]);
//System.out.println("  latticeDimensions[1] = " + latticeDimensions[1]);
//System.out.println("  latticeDimensions[2] = " + latticeDimensions[2]);
// NOTE : Division by 2.0 assumes plane is in center of phase
        if(Math.abs(Math.IEEEremainder(latticeDimensions[planeDimIdx], 2.0)) > 0.95) {
	        latticeDimensions[planeDimIdx]++;
        }
//System.out.println("  latticeDimensions[0] = " + latticeDimensions[0]);
//System.out.println("  latticeDimensions[1] = " + latticeDimensions[1]);
//System.out.println("  latticeDimensions[2] = " + latticeDimensions[2]);

        if (indexIterator.getD() > latticeDimensions.length) {
            int[] iteratorDimensions = new int[latticeDimensions.length+1];
            System.arraycopy(latticeDimensions, 0, iteratorDimensions, 0,
                    latticeDimensions.length);
            iteratorDimensions[latticeDimensions.length] = 1;
            indexIterator.setSize(iteratorDimensions);
        }
        else {
            indexIterator.setSize(latticeDimensions);
        }

        // determine lattice constant
        IVector latticeScaling = phase.getSpace().makeVector();
        if (rescalingToFitVolume) {
            // in favorable situations, this should be approximately equal
            // to 1.0
            latticeScaling.E(phase.getBoundary().getDimensions());
            latticeScaling.DE(latticeConstantV);
            latticeScaling.DE(Space.makeVector(latticeDimensions));
        } else {
            latticeScaling.E(1.0);
        }

        // determine amount to shift lattice so it is centered in volume
        IVector offset = phase.getSpace().makeVector();
        offset.E(phase.getBoundary().getDimensions());
        IVector vectorOfMax = phase.getSpace().makeVector();
        IVector vectorOfMin = phase.getSpace().makeVector();
        IVector site = phase.getSpace().makeVector();
        vectorOfMax.E(Double.NEGATIVE_INFINITY);
        vectorOfMin.E(Double.POSITIVE_INFINITY);

        // XXX this can do strange things. it's probably not needed for 
        // periodic boundaries, but gets the atoms off the boundaries for 
        // non-periodic boundaries
        indexIterator.reset();

        while (indexIterator.hasNext()) {
            site.E((IVector) lattice.site(indexIterator.next()));
            site.TE(latticeScaling);
            for (int i=0; i<site.getD(); i++) {
                vectorOfMax.setX(i, Math.max(site.x(i),vectorOfMax.x(i)));
                vectorOfMin.setX(i, Math.min(site.x(i),vectorOfMin.x(i)));
            }
        }
        offset.Ev1Mv2(vectorOfMax, vectorOfMin);
        offset.TE(-0.5);
        offset.ME(vectorOfMin);
//latticeScaling.setX(0, latticeScaling.x(0) / (plane.getA() + 1.0));
//latticeScaling.setX(1, latticeScaling.x(1) / (plane.getB() + 1.0));
//latticeScaling.setX(2, latticeScaling.x(2) / (plane.getC() + 1.0));
//System.out.println("latticeScaling = " + latticeScaling.x(0));
//System.out.println("latticeScaling = " + latticeScaling.x(1));
//System.out.println("latticeScaling = " + latticeScaling.x(2));
//System.out.println("offset = " + offset.x(0));
//System.out.println("offset = " + offset.x(1));
//System.out.println("offset = " + offset.x(2));

        myLat = new MyLattice(lattice, latticeScaling, offset);

        indexIterator.reset();

        AtomIteratorArrayListSimple[] atomIterator = new AtomIteratorArrayListSimple[numSpecies];

        int maxOnSide = latticeDimensions[0] * latticeDimensions[1] * latticeDimensions[2] / 2;

        // Place molecules
        for (int i = 0; i < numSpecies; i++) {
            atomIterator[i] = new AtomIteratorArrayListSimple(((Species)species.get(i)).getAgent(phase).getChildList());
            atomIterator[i].reset();
        }
//System.out.println("maxOnSide = " + maxOnSide);
//System.out.println("moleculesOnLeft = " + moleculesOnLeft);
//System.out.println("moleculesOnRight = " + moleculesOnRight);

        for(int x = LEFT; x <= RIGHT; x++) {
	        for (int i = 0; i < numSpecies; i++) {

	        	for (int y = 0; y < molecules[x][i]; y++) {
	                IAtom a = atomIterator[i].nextAtom();

			        int[] idx = indexIterator.next();

			        atomActionTranslateTo.setDestination((IVector)myLat.site(idx));
			        atomActionTranslateTo.actionPerformed(a);
	
		        }
	        }
	        if(x == LEFT) {
	        	for(int i = 0; i < (maxOnSide - moleculesOnLeft); i++) {
	        		indexIterator.next();
	        	}
	        }
        }
    }

    public static void main(String[] args) {
        Simulation sim = new Simulation(Space3D.getInstance());
        PotentialMaster potentialMaster = new PotentialMaster(sim.getSpace());
        Phase phase = new Phase(sim);
        sim.addPhase(phase);
        SpeciesSpheresMono species = new SpeciesSpheresMono(sim);
        sim.getSpeciesManager().addSpecies(species);
        ((AtomTypeSphere)species.getMoleculeType()).setDiameter(5.0);
        int k = 4;
        phase.getAgent(species).setNMolecules(4 * k * k * k);
        IntegratorHard integrator = new IntegratorHard(sim, potentialMaster);
        integrator.setPhase(phase);
//	        ColorSchemeByType colorScheme = new ColorSchemeByType();
        // CubicLattice lattice = new LatticeCubicBcc();
        BravaisLatticeCrystal lattice = new LatticeCubicFcc();
        // CubicLattice lattice = new LatticeCubicSimple();
        ConfigurationLattice configuration = new ConfigurationLattice(lattice);
        // phase.boundary().setDimensions(new Space3D.Vector(15.,30.,60.5));
        configuration.initializeCoordinates(phase);
        // etomica.graphics.DisplayPhase display = new
        // etomica.graphics.DisplayPhase(phase);

        etomica.graphics.SimulationGraphic simGraphic = new etomica.graphics.SimulationGraphic(
                sim);
//	        ((ColorSchemeByType) ((DisplayPhase) simGraphic.displayList()
//	                .getFirst()).getColorScheme()).setColor(species
//	                .getMoleculeType(), java.awt.Color.red);
        simGraphic.makeAndDisplayFrame();
    }

}
