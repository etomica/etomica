package etomica;

import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.iterator.AtomIteratorList;

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

    protected boolean zeroTotalMomentum = true;
    protected double[] dimensions;
    
    public Configuration() {}

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
        atom.coord.randomizeMomentum(temperature);
        if(zeroTotalMomentum) {
            Space.Vector zero = (Space.Vector)atom.coord.position().clone();
            zero.E(0.0);
            atom.coord.accelerateTo(zero);
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
        initializeMomenta(group);
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
            //new AtomIteratorSequential(group[i]);
            initializeMomenta(group[i]);
        }
        initializePositions(iterators);
    }

    public static Space1D.Vector[] lineLattice(int n, double Lx) {
        Space1D.Vector[] r = new Space1D.Vector[n];
        double delta = Lx/n;
        for(int i=0; i<n; i++) {
            r[i] = new Space1D.Vector();
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
    public final static Space2D.Vector[] squareLattice(int n, double Lx, double Ly, boolean fillVertical) {
        Space2D.Vector[] r = new Space2D.Vector[n];
        for(int i=0; i<n; i++) {r[i] = new Space2D.Vector();}

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
 
	public final static Space2D.Vector[] hexagonalLattice(int n, double Lx, double Ly, boolean fillVertical) {
		Space2D.Vector[] r = new Space2D.Vector[n];
		if(n == 0) return r;
		Space2D.Vector com = new Space2D.Vector();
		com.E(0.0); // later becomes present Center of Mass
		Space2D.Vector dcom = new Space2D.Vector();
		dcom.E(0.0); // difference in Present COM and Original COM
		Space2D.Vector ocom = new Space2D.Vector();
		ocom.setX(0, Lx);ocom.setX(1, Ly); ocom.TE(0.5); // Original Center Of Mass
                        
		for(int i=0; i<n; i++) {r[i] = new Space2D.Vector();}

		int moleculeColumns, moleculeRows;
		double moleculeInitialSpacingX, moleculeInitialSpacingY;
 
		int i = 0;
		int ix = 0;
		int iy = 0;
		boolean on = true;
		double half = 0.0;

		if(fillVertical){
			moleculeRows = (int)Math.sqrt(Ly/Lx*n*Math.sqrt(3.0)*0.5);
			if(moleculeRows == 0) moleculeRows = 1;
			moleculeColumns = n/moleculeRows;
			if(moleculeColumns == 0) moleculeColumns = 1;
			if(moleculeRows*moleculeColumns < n) moleculeColumns++;
			int rowsDrawn = (int)((double)n/(double)moleculeColumns - 1.0E-10) + 1;
			moleculeInitialSpacingX = Lx/moleculeColumns;
			moleculeInitialSpacingY = Ly / rowsDrawn;
            while (i < n) {
                r[i].setX(0, ix * moleculeInitialSpacingX);
                r[i].setX(1, (iy + half) * moleculeInitialSpacingY);
                i++;
                iy++;
                if (iy >= moleculeRows) {
                    iy = 0;
                    if (on) {
                        half = 0.5;
                        on = false;
                    }
                    else {
                        half = 0.0;
                        on = true;
                    }
                    ix++;
                }
            } //end of while
        }
        else {
            moleculeColumns = (int)Math.sqrt(Lx/Ly*n*Math.sqrt(3.0)*0.5);
            if (moleculeColumns == 0) moleculeColumns = 1;
            moleculeRows = n / moleculeColumns;
            if (moleculeRows == 0) moleculeRows = 1;
            if (moleculeRows * moleculeColumns < n) moleculeRows++;
            int columnsDrawn = (int)((double)n/(double)moleculeRows - 1.0E-10) + 1;
            moleculeInitialSpacingX = Lx/columnsDrawn;
            moleculeInitialSpacingY = Ly/moleculeRows;
            while (i < n) {
                r[i].setX(0, (half + ix) * moleculeInitialSpacingX);
                r[i].setX(1, iy * moleculeInitialSpacingY);
                i++;
                ix++;
                if (ix >= columnsDrawn) {
                    ix = 0;
                    if (on) {
                        half = 0.5;
                        on = false;
                    }
                    else {
                        half = 0.0;
                        on = true;
                    }
                    iy++;
                }
            } //end of while
        }

        for (int j = 0; j < n; j++) {
            com.PE(r[j]);
        }
        com.DE((double) n);
        dcom.E(ocom);
        dcom.ME(com);

		for(int j=0;j<n;j++){
			r[j].PE(dcom);
		}
            
		return r;
                
	}//end of HexagonalLattice   
    
       
    /**
     * Configuration that does nothing to atom positions or momenta.
     */
    private static class Null extends Configuration {
        public void initializePositions(AtomIterator[] iterators) {}
        public void initializeMomenta(Atom atom, double temperature) {}
    }
    public static final Null NULL = new Null();
}//end of Configuration
