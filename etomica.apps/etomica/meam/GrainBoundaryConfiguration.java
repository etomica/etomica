/*
 * Created on Aug 2, 2006
 *
 * Adapted from ConfigurationLattice for the specific case of a two-lattice system,
 * separated by a 2-D grain boundary. latticeA is above the grain boundary, and
 * latticeB is below.  Because the lattice types may be differect, PBC exist only
 * in x- and y-directions.  Along the boundaries perpendicular to the z-direction 
 * are two levels of unit cells filled with fixed atoms (atoms with infinite mass).
 */
package etomica.meam;

import etomica.action.AtomActionTranslateTo;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPositionDefinitionSimple;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.AtomIteratorArrayList;
import etomica.atom.iterator.IteratorDirective;
import etomica.config.Configuration;
import etomica.config.Conformation;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.IndexIteratorSizable;
import etomica.lattice.SpaceLattice;
import etomica.phase.Phase;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space3d.Vector3D;

/**
 * @author K.R. Schadel with help from A. Schultz
 */
public class GrainBoundaryConfiguration extends Configuration {

	/**
     * Constructs class using instances of IndexIteratorRectangular as the default
     * index iterators.
     */
    public GrainBoundaryConfiguration(BravaisLatticeCrystal latticeA, BravaisLatticeCrystal latticeB) {
        this(latticeA, latticeB, new IndexIteratorRectangular(latticeA.D()), 
        		new IndexIteratorRectangular(latticeB.D()));
    }

    /**
     * Construct class that will place atoms on sites of the given lattices,
     * proceeding in the order resulting from iteration through the given index
     * iterator.
     */
    public GrainBoundaryConfiguration(BravaisLatticeCrystal latticeA, BravaisLatticeCrystal latticeB,
            IndexIteratorSizable indexIteratorA, 
			IndexIteratorSizable indexIteratorB) {
    	/** Lattices A + B share same space.  Only need to getSpace() for one 
    	 *  of the two.
    	 */
        super(); 
        this.latticeA = latticeA; 
        this.latticeB = latticeB;
        this.indexIteratorA = indexIteratorA; 
        this.indexIteratorB = indexIteratorB;
        if(indexIteratorA.getD() != latticeA.D()) {
            throw new IllegalArgumentException(
            		"Dimension of index iterator and lattice are incompatible");
        }
        if(indexIteratorB.getD() != latticeB.D()) {
            throw new IllegalArgumentException(
            		"Dimension of index iterator and lattice are incompatible");
        }
        atomActionTranslateTo = new AtomActionTranslateTo(latticeA.getSpace());
        atomActionTranslateTo.setAtomPositionDefinition(new AtomPositionDefinitionSimple());
        //There may be more than one species in the lattice with mobile atoms.
        atomIteratorMobileA = new AtomIteratorArrayList(IteratorDirective.Direction.UP);
        atomIteratorMobileB = new AtomIteratorArrayList(IteratorDirective.Direction.UP);
        //Probably only one type of atom will exit in the lattice with fixed atoms.
        atomIteratorFixedA = new AtomIteratorArrayList(
        		IteratorDirective.Direction.UP);
        atomIteratorFixedB = new AtomIteratorArrayList(
        		IteratorDirective.Direction.UP);
    }

    public void setDimensions (int nCellsAx, int nCellsAy, int nCellsAz, 
    		int nCellsBx, int nCellsBy, int nCellsBz, double aA, double bA,
			double cA, double aB, double bB, double cB) {
    	
    	this.nCellsAx = nCellsAx; 
    	this.nCellsAy = nCellsAy; 
    	this.nCellsAz = nCellsAz; 
    	this.nCellsBx = nCellsBx;
    	this.nCellsBy = nCellsBy;
    	this.nCellsBz = nCellsBz;
    	
    	latticeDimensionsA[0] = nCellsAx * aA;
    	latticeDimensionsA[1] = nCellsAy * bA;
    	latticeDimensionsA[2] = nCellsAz * cA;
    	latticeDimensionsB[0] = nCellsBx * aB;
    	latticeDimensionsB[1] = nCellsBy * bB;
    	latticeDimensionsB[2] = nCellsBz * cB;
    	
    	//iteratorDimensionsA or B
    	iteratorDimensionsA[0] = nCellsAx;
    	iteratorDimensionsA[1] = nCellsAy;
    	iteratorDimensionsA[2] = nCellsAz;
    	iteratorDimensionsA[3] = latticeA.getBasis().getScaledCoordinates().length;
    	iteratorDimensionsB[0] = nCellsBx;
    	iteratorDimensionsB[1] = nCellsBy;
    	iteratorDimensionsB[2] = nCellsBz;
    	iteratorDimensionsB[3] = latticeB.getBasis().getScaledCoordinates().length;
    }	

    /**
     * Places the molecules in the given phase on the positions of the
     * lattice.  
     */
    public void initializeCoordinates(Phase phase) {
    	AtomSet[] lists = getMoleculeLists(phase);
    	IVector firstAtomPosition = ((IAtomPositioned)lists[0].getAtom(0)).getPosition();
    	
    	System.out.println("At beginning of initializePositions  "+ firstAtomPosition);
    	
    	AtomSet[] listsMobileA = new AtomArrayList[lists.length/2 - 1];
    	AtomSet[] listsMobileB = new AtomArrayList[lists.length/2 - 1];
    	
    	for (int i = 1; i < lists.length/2; i++) {
    		listsMobileA[i - 1] = lists[i];
    		listsMobileB[i - 1] = lists[i + lists.length/2];
    	}
    	
        atomIteratorMobileA.setList(lists[1]);
        atomIteratorMobileB.setList(lists[lists.length-1]);
        atomIteratorFixedA.setList(lists[0]);
        atomIteratorFixedB.setList(lists[lists.length/2]);

        indexIteratorA.setSize(iteratorDimensionsA);
        indexIteratorB.setSize(iteratorDimensionsB);

        /**  The offset vectors are used to shift the lattices such that
         *  lattice A is above lattice B, and both are centered on the z axis.
         */
        
        Vector3D offsetA = (Vector3D)phase.getSpace().makeVector();
        Vector3D offsetB = (Vector3D)phase.getSpace().makeVector();
        offsetA.setX(0, -0.5 * latticeDimensionsA[0]);
        offsetA.setX(1, -0.5 * latticeDimensionsA[1]);
        offsetA.setX(2, (latticeDimensionsB[2] - latticeDimensionsA[2])/2.0 );
        offsetB.setX(0, -0.5 * latticeDimensionsB[0]);
        offsetB.setX(1, -0.5 * latticeDimensionsB[1]);
        offsetB.setX(2, -(latticeDimensionsB[2] + latticeDimensionsA[2])/2.0 );
        
        myLatA = new MyLattice(latticeA, offsetA);
        myLatB = new MyLattice(latticeB, offsetB);

        // Place molecules
        atomIteratorMobileA.reset();
        atomIteratorMobileB.reset();
        atomIteratorFixedA.reset();
        atomIteratorFixedB.reset();
        
        indexIteratorA.reset();
        while (indexIteratorA.hasNext()) {
        	//System.out.println("At start of while loop over indexIteratorA  " + firstAtomPosition);
        	IAtom a;
        	int[] ii = indexIteratorA.next();
        	//ii[2] goes from 0 to nCellsAz-1, not 1 to nCellsAz, because of 
        	//how setSize method works
            if (ii[2] > ((nCellsAz - 2) - 1)) {
            	a = atomIteratorFixedA.nextAtom();
            	
            	
            }
            else {
            	a = atomIteratorMobileA.nextAtom();	
            	//System.out.println(ii[2] + "  |  " + a);
            }
            if (a instanceof IAtomGroup) {
                // initialize coordinates of child atoms
                Conformation config = a.getType().creator().getConformation();
                config.initializePositions(((IAtomGroup)a).getChildList());
            }
            IVector site = (IVector) myLatA.site(ii);
            if (((IAtomPositioned)a).getPosition() == firstAtomPosition) {
            	System.out.println();
            }
            atomActionTranslateTo.setDestination(site);
            atomActionTranslateTo.actionPerformed(a);
            System.out.println("A  |  " +a + "  |  " + site.x(2) + "  |  " + ii[2]);
            System.out.println("At end of while loop over indexIteratorA  " + firstAtomPosition);
        }
        
        indexIteratorB.reset();
        while (indexIteratorB.hasNext()) {
        	IAtom a;
        	int[] ii = indexIteratorB.next();
        	//ii[2] goes from 0 to nCellsAz-1, not 1 to nCellsAz
            if (ii[2] < 2) {
            	a = atomIteratorFixedB.nextAtom();
            	
            }
            else {
            	a = atomIteratorMobileB.nextAtom();
            	
            }
            if (a instanceof IAtomGroup) {
                // initialize coordinates of child atoms
                Conformation config = a.getType().creator().getConformation();
                config.initializePositions(((IAtomGroup)a).getChildList());
            }
            IVector site = (IVector) myLatB.site(ii);
            atomActionTranslateTo.setDestination(site);
            atomActionTranslateTo.actionPerformed(a);
            //System.out.println("B  |  " +a + "  |  " + site.x(2) + "  |  " + ii[2]);
        }
        System.out.println("At end of initializePositions()  " + firstAtomPosition);
    }
    
    
    
    /**
     * Returns a lattice with positions the same as those used in the 
     * most recent use of initializeCoordinates.  Includes any scaling
     * or translation applied to fill the space, and thus will not necessarily
     * be the same positions as specified by the lattice given at construction.
     */
    public SpaceLattice getLatticeAMemento() {
        return myLatA;
    }
    
    public SpaceLattice getLatticeBMemento() {
        return myLatB;
    }
    
    /**
     * Used to store the state of a lattice.
     * 
     * @author nancycribbin, Andrew Schultz, Dr. Kofke
     * 
     */
    private static class MyLattice implements SpaceLattice {

        MyLattice(SpaceLattice l, IVector offset) {
            lattice = l;
            this.offset = offset;
            this.site = l.getSpace().makeVector();
        }

        public Space getSpace() {
            return lattice.getSpace();
        }

        public int D() {
            return lattice.D();
        }

        /**
         * Returns the same IVector instance with each call.
         */
        public Object site(int[] index) {
            site.E((IVector) lattice.site(index));
            site.PE(offset);
            return site;
        }

        public double[] getLatticeConstants() {
            double[] lat = lattice.getLatticeConstants();
            return lat;
        }

        SpaceLattice lattice;
        IVector offset;
        private final IVector site;
    }

    private final BravaisLatticeCrystal latticeA, latticeB;
    private final IndexIteratorSizable indexIteratorA, indexIteratorB;
    private final AtomActionTranslateTo atomActionTranslateTo;
    private final AtomIteratorArrayList atomIteratorFixedA, atomIteratorFixedB;
    private final AtomIteratorArrayList atomIteratorMobileA, atomIteratorMobileB;
    int[] iteratorDimensionsA = new int[4];
    int[] iteratorDimensionsB = new int[4];
    //int[] indexIteratorA = new int[4];
    //int[] indexIteratorB = new int[4];
    double[] latticeDimensionsA = new double[3];
    double[] latticeDimensionsB = new double[3];
    private int nCellsAx, nCellsAy, nCellsAz, nCellsBx, nCellsBy, nCellsBz;
    private MyLattice myLatA, myLatB;
    private static final long serialVersionUID = 2L;
}