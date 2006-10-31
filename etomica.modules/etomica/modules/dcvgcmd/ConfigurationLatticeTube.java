package etomica.modules.dcvgcmd;

import etomica.action.AtomActionTranslateTo;
import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.iterator.AtomIteratorArrayListCompound;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.config.Configuration;
import etomica.config.Conformation;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayPhase;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.IndexIteratorSequential;
import etomica.lattice.IndexIteratorSizable;
import etomica.lattice.LatticeCubicFcc;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

/**
 * Creates a configuration using a CubicLattice to specify positions.  Has
 * capability to assign lattice site to atoms when specifying their coordinates.
 * See setAssigningSitesToAtoms method.
 */
public class ConfigurationLatticeTube extends Configuration {

    public ConfigurationLatticeTube(BravaisLatticeCrystal lattice, double length, SpeciesTube tubeSpecies) {
        this(lattice, length, tubeSpecies, new IndexIteratorSequential(lattice.D()));//need a default iterator
    }
	/**
	 * Constructor for ConfigurationLatticeTube.
	 * @param space
	 */
	public ConfigurationLatticeTube(BravaisLatticeCrystal lattice, double length, SpeciesTube tubeSpecies, IndexIteratorSizable indexIterator) {
	    super(lattice.getSpace());
        this.lattice = lattice;
        this.indexIterator = indexIterator;
        this.length = length;
        this.tubeSpecies = tubeSpecies;
        atomActionTranslateTo = new AtomActionTranslateTo(lattice.getSpace());
        atomIterator = new AtomIteratorArrayListCompound();
        work = space.makeVector();
	}
	
    public void initializePositions(AtomArrayList[] lists) {
        if(lists.length == 0) return;
        atomIterator.setLists(lists);
        int sumOfMolecules = atomIterator.size();
        if(sumOfMolecules == 0) {return;}
        if(sumOfMolecules == 1) {
            atomIterator.reset();
            work.E(0.0);
            atomActionTranslateTo.setDestination(work);
            atomActionTranslateTo.actionPerformed(atomIterator.nextAtom());
            return;
        }
        
        atomIterator.reset();        
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.nextAtom();
            if (a.type.getSpecies()==tubeSpecies){
            	sumOfMolecules--;
            }
        }
        int basisSize = lattice.getBasis().getScaledCoordinates().length;
        int nCells = (int)Math.ceil((double)sumOfMolecules/(double)basisSize);
        
        //determine scaled shape of simulation volume
        double[] latticeConstant = lattice.getPrimitive().getSize();
        double[] shape = (double[])dimensions.clone();
        shape[2] = dimensions[2]*length;
        for(int i=0; i<shape.length; i++) {
            shape[i] /= latticeConstant[i];
        }

                
        //determine number of cells in each direction
        int[] latticeDimensions = calculateLatticeDimensions((nCells+1)/2, shape);
        int[] iteratorDimensions = new int[latticeDimensions.length + 1];
        System.arraycopy(latticeDimensions, 0, iteratorDimensions, 0, latticeDimensions.length);
        iteratorDimensions[latticeDimensions.length] = basisSize;
        indexIterator.setSize(iteratorDimensions);
    
        //determine lattice constant
        Vector latticeScaling = space.makeVector();
        for(int i=0; i<latticeDimensions.length-1; i++) {
            latticeScaling.setX(i,dimensions[i]/(latticeDimensions[i]*latticeConstant[i]));
        }
        latticeScaling.setX(2,dimensions[2]*length/(latticeDimensions[2]*latticeConstant[2]));

        //determine amount to shift lattice so it is centered in volume
        double[] offset = (double[])dimensions.clone();
        double[] vectorOfMax = new double[lattice.getSpace().D()]; 
        double[] vectorOfMin = new double[lattice.getSpace().D()]; 
        for(int i=0; i<lattice.getSpace().D(); i++) {
            vectorOfMax[i] = Double.NEGATIVE_INFINITY;
            vectorOfMin[i] = Double.POSITIVE_INFINITY;
        }

        //XXX this looks scary and asking for trouble
        indexIterator.reset();
        while(indexIterator.hasNext()) {
            Vector site = (Vector)lattice.site(indexIterator.next());
            site.TE(latticeScaling);
            for(int i=0; i<site.D(); i++) {
                vectorOfMax[i] = Math.max(vectorOfMax[i],site.x(i));
                vectorOfMin[i] = Math.min(vectorOfMin[i],site.x(i));
            }
        }
        for(int i=0; i<lattice.getSpace().D(); i++) {
            offset[i] = -0.5*(vectorOfMax[i]-vectorOfMin[i]) - vectorOfMin[i];
        }
        offset[2] -= 0.5*dimensions[2]*(1-length);
                
        Vector offsetVector = Space.makeVector(offset);
        
        // Place molecules  
        atomIterator.reset();
        indexIterator.reset();
        int counterNAtoms = 0;
        int halfwaypoint = atomIterator.size()/2; 
        
        AtomArrayList list = new AtomArrayList();
                
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.nextAtom();
            if (a.type.getSpecies()==tubeSpecies){
            	list.add(a);
            	continue;
            }
            if (!a.node.isLeaf()) {
                //initialize coordinates of child atoms
                Conformation config = a.type.creator().getConformation();
                config.initializePositions(((AtomTreeNodeGroup)a.node).childList);
            }
            
            if (counterNAtoms == halfwaypoint){
            	double z = offsetVector.x(2);
                offsetVector.setX(2,z+dimensions[2]*(1-length));
                indexIterator.reset();
            }
            
            counterNAtoms++;
            Vector site = (Vector)lattice.site(indexIterator.next());
            site.TE(latticeScaling);
            site.PE(offsetVector);
            atomActionTranslateTo.setDestination(site);
            atomActionTranslateTo.actionPerformed(a);
        }
        //loop for multiple tubes.
        AtomIteratorArrayListSimple tubeiterator = new AtomIteratorArrayListSimple(list);
        tubeiterator.reset();
        atomActionTranslateTo.setAtomPositionDefinition(new AtomPositionGeometricCenter(space));
        while (tubeiterator.hasNext()){
        	Atom a = tubeiterator.nextAtom();
        	Conformation config = a.type.creator().getConformation();
            config.initializePositions(((AtomTreeNodeGroup)a.node).childList);
            Vector site = space.makeVector();
            atomActionTranslateTo.setDestination(site);
            atomActionTranslateTo.actionPerformed(a);
        }
        
    }
    
    private int[] calculateLatticeDimensions(int nCells, double[] shape) {
        int dimLeft = shape.length;
        int nCellsLeft = nCells;
        int[] latticeDimensions = new int[shape.length];
        while (dimLeft > 0) {
            double smin = Double.POSITIVE_INFINITY;
            int dmin = 0;
            double product = 1.0;
            for (int idim=0; idim<shape.length; idim++) {
                if (latticeDimensions[idim] > 0) continue;
                if (shape[idim] < smin) {
                    smin = shape[idim];
                    dmin = idim;
                }
                product *= shape[idim];
            }
            // round off except for last dimension (then round up)
            if (dimLeft > 1) {
                latticeDimensions[dmin] = (int)Math.round(shape[dmin]*Math.pow((nCellsLeft/product),1.0/dimLeft));
            }
            else {
                latticeDimensions[dmin] = (int)Math.ceil(shape[dmin]*nCellsLeft/product);
            }
            nCellsLeft = (nCellsLeft + latticeDimensions[dmin] - 1) / latticeDimensions[dmin];
            dimLeft--;
        }
        return latticeDimensions;
    }
	
	private final BravaisLatticeCrystal lattice;
    private final IndexIteratorSizable indexIterator;
    private final Vector work;
    private final AtomActionTranslateTo atomActionTranslateTo;
    private final AtomIteratorArrayListCompound atomIterator;
    private final double length;
    private final SpeciesTube tubeSpecies;

	public static void main(String[] args) {
        Simulation sim = new Simulation(Space3D.getInstance());
		sim.getDefaults().atomSize = 5.0;
		Phase phase = new Phase(sim);
		SpeciesSpheresMono species = new SpeciesSpheresMono(sim);
		int k = 4;
		phase.getAgent(species).setNMolecules(4*k*k*k);
//        CubicLattice lattice = new LatticeCubicBcc();
        BravaisLatticeCrystal lattice = new LatticeCubicFcc();
//        CubicLattice lattice = new LatticeCubicSimple();
		ConfigurationLatticeTube configuration = new ConfigurationLatticeTube(lattice, .25, new SpeciesTube(sim, 10,10));
//        phase.boundary().setDimensions(new Space3D.Vector(15.,30.,60.5));
        configuration.initializeCoordinates(phase);
//		etomica.graphics.DisplayPhase display = new etomica.graphics.DisplayPhase(phase);
		
        etomica.graphics.SimulationGraphic simGraphic = new etomica.graphics.SimulationGraphic(sim);
        ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayPhase)simGraphic.displayList().getFirst()).getColorScheme());
        colorScheme.setColor(species.getMoleculeType(), java.awt.Color.red);
		simGraphic.makeAndDisplayFrame();
	}

}
