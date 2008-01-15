package etomica.modules.dcvgcmd;

import etomica.action.AtomActionTranslateTo;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.IMolecule;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.config.Conformation;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.IndexIteratorSizable;
import etomica.lattice.LatticeCubicFcc;
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
	
	public void setSpeciesSpheres(SpeciesSpheresMono[] speciesSpheres) {
	    this.speciesSpheres = speciesSpheres;
	}
	
	public void setSpeciesTube(SpeciesTube speciesTube) {
	    this.speciesTube = speciesTube;
	}
	
    public void initializeCoordinates(Box box) {
        AtomSet[] spheresLists = new AtomSet[]{box.getMoleculeList(speciesSpheres[0]), box.getMoleculeList(speciesSpheres[1])};
        
        int basisSize = 1;
        if (lattice instanceof BravaisLatticeCrystal) {
            basisSize = ((BravaisLatticeCrystal)lattice).getBasis().getScaledCoordinates().length;
        }
        int nCells = (int)Math.ceil((double)spheresLists[0].getAtomCount()/(double)basisSize);
        
        //determine scaled shape of simulation volume
        IVector shape = box.getSpace().makeVector();
        shape.E(box.getBoundary().getDimensions());
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
        IVector latticeScaling = box.getSpace().makeVector();
        if (rescalingToFitVolume) {
            // in favorable situations, this should be approximately equal
            // to 1.0
            latticeScaling.E(shape);
            latticeScaling.DE(Space.makeVector(latticeDimensions));
        } else {
            latticeScaling.E(1.0);
        }

        // determine amount to shift lattice so it is centered in volume
        IVector offset = box.getSpace().makeVector();
        offset.E(box.getBoundary().getDimensions());
        IVector vectorOfMax = box.getSpace().makeVector();
        IVector vectorOfMin = box.getSpace().makeVector();
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
        offset.setX(2, offset.x(2) - 0.5*box.getBoundary().getDimensions().x(2)*(1-length));

        myLat = new MyLattice(lattice, latticeScaling, offset);

        // Place molecules  
        indexIterator.reset();
        
        // first species (mono spheres)
        int nSpheres = spheresLists[0].getAtomCount();
        for (int i=0; i<nSpheres; i++) {
            IAtom a = spheresLists[0].getAtom(i);
            
            int[] ii = indexIterator.next();
            IVector site = (IVector) myLat.site(ii);
            atomActionTranslateTo.setDestination(site);
            atomActionTranslateTo.actionPerformed(a);
        }
        
        double z = offset.x(2);
        offset.setX(2,z+box.getBoundary().getDimensions().x(2)*(1-length));
        myLat = new MyLattice(lattice, latticeScaling, offset);
        indexIterator.reset();
        
        nSpheres = spheresLists[1].getAtomCount();
        // second species (mono spheres)
        for (int i=0; i<nSpheres; i++) {
            IAtom a = spheresLists[1].getAtom(i);
            
            int[] ii = indexIterator.next();
            IVector site = (IVector) myLat.site(ii);
            atomActionTranslateTo.setDestination(site);
            atomActionTranslateTo.actionPerformed(a);
        }
        
        //loop for multiple tubes.
        AtomSet tubeList = box.getMoleculeList(speciesTube);
        int nTubes = tubeList.getAtomCount();
        atomActionTranslateTo.setAtomPositionDefinition(new AtomPositionGeometricCenter(box.getSpace()));
        // put them all at 0.  oops
        atomActionTranslateTo.setDestination(box.getSpace().makeVector());
        for (int i=0; i<nTubes; i++) {
            IMolecule a = (IMolecule)tubeList.getAtom(i);
        	Conformation config = ((AtomTypeMolecule)a.getType()).getConformation();
            config.initializePositions(a.getChildList());
            atomActionTranslateTo.actionPerformed(a);
        }
        
    }
    
    private static final long serialVersionUID = 1L;
    private final IndexIteratorSizable indexIterator;
    private final AtomActionTranslateTo atomActionTranslateTo;
    protected SpeciesSpheresMono[] speciesSpheres;
    protected SpeciesTube speciesTube;
    private final double length;

	public static void main(String[] args) {
        Simulation sim = new Simulation(Space3D.getInstance());
		Box box = new Box(sim);
        sim.addBox(box);
        SpeciesSpheresMono species1 = new SpeciesSpheresMono(sim);
		SpeciesSpheresMono species2 = new SpeciesSpheresMono(sim);
        sim.getSpeciesManager().addSpecies(species1);
        sim.getSpeciesManager().addSpecies(species2);
        ((AtomTypeSphere)species1.getLeafType()).setDiameter(3.0);
        ((AtomTypeSphere)species2.getLeafType()).setDiameter(3.0);
		int k = 4;
		box.setNMolecules(species1, 2*k*k*k);
        box.setNMolecules(species2, 2*k*k*k);
        SpeciesTube speciesTube = new SpeciesTube(sim, 10, 10);
        sim.getSpeciesManager().addSpecies(speciesTube);
        ((AtomTypeSphere)speciesTube.getMoleculeType().getChildTypes()[0]).setDiameter(3.0);
        
        box.setNMolecules(speciesTube, 1);
//        CubicLattice lattice = new LatticeCubicBcc();
        BravaisLatticeCrystal lattice = new LatticeCubicFcc();
//        CubicLattice lattice = new LatticeCubicSimple();
		ConfigurationLatticeTube configuration = new ConfigurationLatticeTube(lattice, .25);
//        box.boundary().setDimensions(new Space3D.Vector(15.,30.,60.5));
        configuration.initializeCoordinates(box);
//		etomica.graphics.DisplayBox display = new etomica.graphics.DisplayBox(box);
		
        etomica.graphics.SimulationGraphic simGraphic = new etomica.graphics.SimulationGraphic(sim);
        simGraphic.add(new DisplayBox(box));
        ColorSchemeByType colorScheme = (ColorSchemeByType)simGraphic.getDisplayBox(box).getColorScheme();
        colorScheme.setColor(species1.getMoleculeType(), java.awt.Color.blue);
        colorScheme.setColor(species2.getMoleculeType(), java.awt.Color.white);
		simGraphic.makeAndDisplayFrame();
	}

}
