package etomica.modules.dcvgcmd;

import etomica.action.AtomActionTranslateTo;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeGroup;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.config.ConfigurationLattice;
import etomica.config.Conformation;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayPhase;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.IndexIteratorSizable;
import etomica.lattice.LatticeCubicFcc;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

/**
 * Creates a configuration using a CubicLattice to specify positions.  Has
 * capability to assign lattice site to atoms when specifying their coordinates.
 * See setAssigningSitesToAtoms method.
 */
public class ConfigurationLatticeTube extends ConfigurationLattice {

    public ConfigurationLatticeTube(BravaisLatticeCrystal lattice, double length) {
        this(lattice, length, new IndexIteratorRectangular(lattice.D()));//need a default iterator
    }
	/**
	 * Constructor for ConfigurationLatticeTube.
	 * @param space
	 */
	public ConfigurationLatticeTube(BravaisLatticeCrystal lattice, double length, IndexIteratorSizable indexIterator) {
	    super(lattice);
        this.indexIterator = indexIterator;
        this.length = length;
        atomActionTranslateTo = new AtomActionTranslateTo(lattice.getSpace());
	}
	
    public void initializeCoordinates(Phase phase) {
        AtomSet[] lists = getMoleculeLists(phase);
        if(lists.length == 0 || lists[0].getAtomCount() == 0) return;
        
        int basisSize = 1;
        if (lattice instanceof BravaisLatticeCrystal) {
            basisSize = ((BravaisLatticeCrystal)lattice).getBasis().getScaledCoordinates().length;
        }
        int nCells = (int)Math.ceil((double)lists[0].getAtomCount()/(double)basisSize);
        
        //determine scaled shape of simulation volume
        IVector shape = phase.getSpace().makeVector();
        shape.E(phase.getBoundary().getDimensions());
        shape.setX(2,shape.x(2)*length);
        IVector latticeConstantV = Space.makeVector(lattice.getLatticeConstants());
        shape.DE(latticeConstantV);

        // determine number of cells in each direction
        int[] latticeDimensions = calculateLatticeDimensions(nCells, shape);
        if (indexIterator.getD() > latticeDimensions.length) {
            int[] iteratorDimensions = new int[latticeDimensions.length+1];
            System.arraycopy(latticeDimensions, 0, iteratorDimensions, 0,
                    latticeDimensions.length);
            iteratorDimensions[latticeDimensions.length] = basisSize;
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
            latticeScaling.E(shape);
            latticeScaling.DE(Space.makeVector(latticeDimensions));
        } else {
            latticeScaling.E(1.0);
        }

        // determine amount to shift lattice so it is centered in volume
        IVector offset = phase.getSpace().makeVector();
        offset.E(phase.getBoundary().getDimensions());
        IVector vectorOfMax = phase.getSpace().makeVector();
        IVector vectorOfMin = phase.getSpace().makeVector();
        vectorOfMax.E(Double.NEGATIVE_INFINITY);
        vectorOfMin.E(Double.POSITIVE_INFINITY);

        // XXX this can do strange things. it's probably not needed for 
        // periodic boundaries, but gets the atoms off the boundaries for 
        // non-periodic boundaries
        indexIterator.reset();
        while (indexIterator.hasNext()) {
            IVector site = (IVector) lattice.site(indexIterator.next());
            site.TE(latticeScaling);
            for (int i=0; i<site.getD(); i++) {
                vectorOfMax.setX(i, Math.max(site.x(i),vectorOfMax.x(i)));
                vectorOfMin.setX(i, Math.min(site.x(i),vectorOfMin.x(i)));
            }
        }
        offset.Ev1Mv2(vectorOfMax, vectorOfMin);
        offset.TE(-0.5);
        offset.ME(vectorOfMin);
        offset.setX(2, offset.x(2) - 0.5*phase.getBoundary().getDimensions().x(2)*(1-length));

        myLat = new MyLattice(lattice, latticeScaling, offset);

        // Place molecules  
        indexIterator.reset();
        
        // first species (mono spheres)
        int nSpheres = lists[0].getAtomCount();
        for (int i=0; i<nSpheres; i++) {
            IAtom a = lists[0].getAtom(i);
            
            int[] ii = indexIterator.next();
            IVector site = (IVector) myLat.site(ii);
            atomActionTranslateTo.setDestination(site);
            atomActionTranslateTo.actionPerformed(a);
        }
        
        double z = offset.x(2);
        offset.setX(2,z+phase.getBoundary().getDimensions().x(2)*(1-length));
        myLat = new MyLattice(lattice, latticeScaling, offset);
        indexIterator.reset();
        
        nSpheres = lists[1].getAtomCount();
        // second species (mono spheres)
        for (int i=0; i<nSpheres; i++) {
            IAtom a = lists[1].getAtom(i);
            
            int[] ii = indexIterator.next();
            IVector site = (IVector) myLat.site(ii);
            atomActionTranslateTo.setDestination(site);
            atomActionTranslateTo.actionPerformed(a);
        }
        
        //loop for multiple tubes.
        int nTubes = lists[2].getAtomCount();
        atomActionTranslateTo.setAtomPositionDefinition(new AtomPositionGeometricCenter(phase.getSpace()));
        // put them all at 0.  oops
        atomActionTranslateTo.setDestination(phase.getSpace().makeVector());
        for (int i=0; i<nTubes; i++) {
            IAtomGroup a = (IAtomGroup)lists[2].getAtom(i);
        	Conformation config = ((AtomTypeGroup)a.getType()).getConformation();
            config.initializePositions(a.getChildList());
            atomActionTranslateTo.actionPerformed(a);
        }
        
    }
    
    private static final long serialVersionUID = 1L;
    private final IndexIteratorSizable indexIterator;
    private final AtomActionTranslateTo atomActionTranslateTo;
    private final double length;

	public static void main(String[] args) {
        Simulation sim = new Simulation(Space3D.getInstance());
		Phase phase = new Phase(sim);
        sim.addPhase(phase);
        SpeciesSpheresMono species1 = new SpeciesSpheresMono(sim);
		SpeciesSpheresMono species2 = new SpeciesSpheresMono(sim);
        sim.getSpeciesManager().addSpecies(species1);
        sim.getSpeciesManager().addSpecies(species2);
        ((AtomTypeSphere)species1.getMoleculeType()).setDiameter(3.0);
        ((AtomTypeSphere)species2.getMoleculeType()).setDiameter(3.0);
		int k = 4;
		phase.getAgent(species1).setNMolecules(2*k*k*k);
        phase.getAgent(species2).setNMolecules(2*k*k*k);
        SpeciesTube speciesTube = new SpeciesTube(sim, 10, 10);
        sim.getSpeciesManager().addSpecies(speciesTube);
        ((AtomTypeSphere)((AtomTypeGroup)speciesTube.getMoleculeType()).getChildTypes()[0]).setDiameter(3.0);
        
        phase.getAgent(speciesTube).setNMolecules(1);
//        CubicLattice lattice = new LatticeCubicBcc();
        BravaisLatticeCrystal lattice = new LatticeCubicFcc();
//        CubicLattice lattice = new LatticeCubicSimple();
		ConfigurationLatticeTube configuration = new ConfigurationLatticeTube(lattice, .25);
//        phase.boundary().setDimensions(new Space3D.Vector(15.,30.,60.5));
        configuration.initializeCoordinates(phase);
//		etomica.graphics.DisplayPhase display = new etomica.graphics.DisplayPhase(phase);
		
        etomica.graphics.SimulationGraphic simGraphic = new etomica.graphics.SimulationGraphic(sim);
        simGraphic.add(new DisplayPhase(phase));
        ColorSchemeByType colorScheme = (ColorSchemeByType)simGraphic.getDisplayPhase(phase).getColorScheme();
        colorScheme.setColor(species1.getMoleculeType(), java.awt.Color.blue);
        colorScheme.setColor(species2.getMoleculeType(), java.awt.Color.white);
		simGraphic.makeAndDisplayFrame();
	}

}
