package etomica;

import etomica.action.AtomActionAccelerateTo;
import etomica.action.AtomActionRandomizeVelocity;
import etomica.action.AtomGroupAction;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.iterator.AtomIteratorList;
import etomica.space.Vector;
import etomica.space1d.Space1D;

/**
 * General class for assignment of coordinates to a group of atoms.
 * 
 * @author David Kofke
 */
 
 /* History of changes
  * 09/01/02 (DAK) added field to flag whether total momentum should be set to
  * zero when initalizing (made in conjunction with change to Space.
  * CoordinateGroup classes, which no longer do randomizeMomentum to zero total
  * momentum). 
  * 09/04/02 (DAK) added Null inner class and NULL field
  * 01/21/04 (DAK) added initializeCoordinate(Phase) method
  */
public abstract class Configuration implements java.io.Serializable {

    protected boolean zeroTotalMomentum = false;
    protected double[] dimensions;
    protected AtomGroupAction atomActionRandomizeVelocity;
    protected AtomActionAccelerateTo atomActionAccelerateTo;
    protected Space space;
    protected Vector work;
    
    public Configuration(Space space) {
        this.space = space;
        atomActionRandomizeVelocity = new AtomGroupAction(new AtomActionRandomizeVelocity());
        atomActionAccelerateTo = new AtomActionAccelerateTo(space);
        work = space.makeVector();
    }

    public abstract void initializePositions(AtomIterator[] iterator);
    
    public void initializeCoordinates(Phase phase) {
    	setDimensions(phase.boundary().dimensions().toArray());
		initializeCoordinates(phase.speciesMaster.node.childAtomArray());
    }
    public void initializePositions(AtomIterator iterator) {
        initializePositions(new AtomIterator[] {iterator});
    }
    
    public void setDimensions(double[] dim) {
        dimensions = (double[])dim.clone();
    }

    public double[] getDimensions() {return dimensions;}
    
    /** 
     * Sets flag indicating if total momentum should be zeroed after initializing
     * sub-atom momenta.  Default is true.
     */
    public void setZeroTotalMomentum(boolean b) {zeroTotalMomentum = b;}
    /** 
     * Flag indicating if total momentum should be zeroed after initializing
     * sub-atom momenta.  Default is true.
     */
    public boolean isZeroTotalMomentum() {return zeroTotalMomentum;}
 /**   
  * All atom velocities are set such that all have the same total momentum (corresponding to
  * the default value of temperature), with the direction at random
  */
    public void initializeMomenta(Atom atom) {
        initializeMomenta(atom, Default.TEMPERATURE);
    }
    public void initializeMomenta(Atom atom, double temperature) {
        ((AtomActionRandomizeVelocity)atomActionRandomizeVelocity.getAction()).setTemperature(temperature);
        atomActionRandomizeVelocity.actionPerformed(atom);
        if(zeroTotalMomentum) {
            work.E(0.0);
            atomActionAccelerateTo.setTargetVelocity(work);
            atomActionAccelerateTo.actionPerformed(atom);
        }
    }//end of initializeMomenta
    
    /**
     * Initializes positions and momenta of the atoms in the given atom group.
     */
    public void initializeCoordinates(Atom group) {
        if(group.node.isLeaf()) throw new IllegalArgumentException("Error in Configuration.initializeCoordinates:  Attempt to initialize child coordinates of leaf atom");

        AtomIteratorList iterator = 
                AtomIteratorList.makeCopylistIterator(((AtomTreeNodeGroup)group.node).childList);
        initializePositions(iterator);
        if (space.isKinetic()) {
            initializeMomenta(group);
        }
    }
    
    /**
     * Initializes positions of the atoms in the given atom group.
     */
    public void initializePositions(Atom group) {
        if(group.node.isLeaf()) throw new IllegalArgumentException("Error in Configuration.initializeCoordinates:  Attempt to initialize child coordinates of leaf atom");
        initializePositions(new AtomIteratorList(((AtomTreeNodeGroup)group.node).childList));
    }
    /**
     * Initializes positions and momenta of the atoms in the given atom groups.
     */
    public void initializeCoordinates(Atom[] group) {
        AtomIterator[] iterators = new AtomIterator[group.length];
        for(int i=0; i<group.length; i++) {
            iterators[i] = AtomIteratorList.makeCopylistIterator(((AtomTreeNodeGroup)group[i].node).childList);
            if (space.isKinetic()) {
                initializeMomenta(group[i]);
            }
        }
        initializePositions(iterators);
    }

    public static etomica.space1d.Vector1D[] lineLattice(int n, double Lx) {
        etomica.space1d.Vector1D[] r = new etomica.space1d.Vector1D[n];
        double delta = Lx/n;
        for(int i=0; i<n; i++) {
            r[i] = new etomica.space1d.Vector1D();
            r[i].setX(0, (i+0.5)*delta);
        }
        return r;
    }
        
    /**
     * Returns a set of n coordinates filling a square lattice of sides Lx and Ly
     * If n is not suitable for square lattice, then last sites are left unfilled
     * Lattice is centered between (0,0) and (Lx, Ly).
     * The final argument should be passed one of the class variables VERTICAL or HORIZONTAL, indicating
     *   whether successive points fill the lattice across or down.
     */
    public final static etomica.space2d.Vector2D[] squareLattice(int n, double Lx, double Ly, boolean fillVertical) {
        etomica.space2d.Vector2D[] r = new etomica.space2d.Vector2D[n];
        for(int i=0; i<n; i++) {r[i] = new etomica.space2d.Vector2D();}

        int moleculeColumns, moleculeRows;
        double moleculeInitialSpacingX, moleculeInitialSpacingY;

    //  Number of molecules per row (moleculeColumns) and number of rows (moleculeRows)
    //  in initial configuration
        moleculeColumns = (int)Math.sqrt(Lx/Ly*n);
        moleculeRows = n/moleculeColumns;

        if(moleculeRows*moleculeColumns < n) moleculeRows++;

    //moleculeColumns may be greater than the actual number of columns drawn
    //Need to center columns in the initial position.
        int columnsDrawn = (int)((double)n/(double)moleculeRows - 1.0E-10) + 1;
        
    //moleculeColumnsShift used to center initial coordinates
        double moleculeColumnsShift = Lx/columnsDrawn/2;
        double moleculeRowsShift = Ly/moleculeRows/2;
  
    //assign distance between molecule centers
        moleculeInitialSpacingX = Lx/columnsDrawn;
        moleculeInitialSpacingY = Ly/moleculeRows;
        int i = 0;
        int ix = 0;
        int iy = 0;
        while(i < n) {
            r[i].setX(0, ix*moleculeInitialSpacingX + moleculeColumnsShift);
	        r[i].setX(1, iy*moleculeInitialSpacingY + moleculeRowsShift);
	        i++;
	        if(fillVertical) {
	            iy++;
	            if(iy >= moleculeRows) {
	                iy = 0;
	                ix++;
	            }}
	        else {
	            ix++;
	            if(ix >= columnsDrawn) {
	                ix = 0;
	                iy++;
	            }
	        }
	    }
	    return r;
    }//end of squareLattice
 
	public final static etomica.space2d.Vector2D[] hexagonalLattice(int n, double Lx, double Ly, boolean fillVertical) {
		etomica.space2d.Vector2D[] r = new etomica.space2d.Vector2D[n];
		if(n == 0) return r;
                        
		for(int i=0; i<n; i++) {r[i] = new etomica.space2d.Vector2D();}

		int[] moleculeL = new int[2];
        double[] L = new double[2];
        L[0] = Lx;
        L[1] = Ly;
 
        int fillDim = fillVertical ? 1 : 0;
        // ratio of moleculeRows/moleculeColumns (or vica versa)
        double ratio = L[fillDim]/L[1-fillDim]*Math.sqrt(3.0)*0.5;
        // round off lower of the two and round up the other, ensuring both are even
        if (ratio>1.0) {
            moleculeL[1-fillDim] = 2*Math.max(1,(int)Math.round(Math.sqrt(n/ratio)/2.0));
            moleculeL[fillDim] = ((n + 2*moleculeL[1-fillDim] - 1) / (2*moleculeL[1-fillDim])) * 2;
        }
        else {
            moleculeL[fillDim] = 2*Math.max(1,(int)Math.round(Math.sqrt(n*ratio)/2.0));
            moleculeL[1-fillDim] = ((n + 2*moleculeL[fillDim]- 1) / (2*moleculeL[fillDim])) * 2;
        }
        // pixel spacing between molecules
        double[] moleculeInitialSpacing = new double[2];
        moleculeInitialSpacing[0] = L[0] / moleculeL[0];
        moleculeInitialSpacing[1] = L[1] / moleculeL[1];
        double half = 0.0;
        // counters for each direction
        int[] j = new int[2];
        for (int i=0; i<n; i++) {
            r[i].setX(1-fillDim, (0.5 + j[1-fillDim]) * moleculeInitialSpacing[1-fillDim]);
            r[i].setX(fillDim,   (0.25 + j[fillDim] + half) * moleculeInitialSpacing[fillDim]);
            j[fillDim]++;
            if (j[fillDim] == moleculeL[fillDim]) {
                j[fillDim] = 0;
                if (half == 0.0) {
                    half = 0.5;
                }
                else {
                    half = 0.0;
                }
                j[1-fillDim]++;
            }
        }
           
		return r;
    }   
    
       
    /**
     * Configuration that does nothing to atom positions or momenta.
     */
    private static class Null extends Configuration {
        Null() {
            super(new Space1D());
        }
        public void initializePositions(AtomIterator[] iterators) {}
        public void initializeMomenta(Atom atom, double temperature) {}
    }
    public static final Null NULL = new Null();
}//end of Configuration
