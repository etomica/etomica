package etomica;
import java.util.*;
import etomica.units.*;

/**
 * General class for assignment of coordinates to a group of atoms.
 */
public abstract class Configuration implements java.io.Serializable {

//    protected Vector species = new Vector();
    protected double temperature = Default.TEMPERATURE;
    Random rand = new Random();
    protected final Space.Vector work;
    protected final Space space;
    
    public Configuration(Space space) {
        this.space = space;
        work = space.makeVector();
    }

    public final void setTemperature(double t) {temperature = t;}
    public final double getTemperature() {return temperature;}
    public final Dimension getTemperatureDimension() {return Dimension.TEMPERATURE;}
 /**   
  * All atom velocities are set such that all have the same total momentum (corresponding to
  * the current value of temperature), with the direction at random
  */
    public void initializeMomenta(Atom atom) {
        
        atom.coord.randomizeMomentum(temperature);
        
        /**   zero com momentum here or in coordinate?*/
  /*      work.E(0.0);
        int sum = 0;
        AtomIterator iterator = ((AtomGroup)atom).childIterator;
        iterator.reset();
        while(iterator.hasNext()) {
            Atom a = iterator.next();
            a.coord.randomizeMomentum(temperature);
            sum++;
            work.PE(a.coord.momentum());
	    }
	    if(sum < 2) return;  //don't zero COM momentum if only one atom
    //    Zero center-of-mass momentum
        work.DE((double)sum);
        iterator.reset();
        while(iterator.hasNext()) {
            iterator.next().coord.momentum().ME(work);
        }*/
    }//end of initializeMomenta
        
    
    public abstract void initializeCoordinates(Atom group);
    
    public final static boolean HORIZONTAL = false;
    public final static boolean VERTICAL = true;

    /* commented for redesign
    public static Space1D.Vector[] lineLattice(int n, double Lx) {
        Space1D.Vector[] r = new Space1D.Vector[n];
        double delta = Lx/(double)n;
        for(int i=0; i<n; i++) {
            r[i] = new Space1D.Vector();
            r[i].x = (i+0.5)*delta;
        }
        return r;
    }*/
        
    /**
     * Returns a set of n coordinates filling a square lattice of sides Lx and Ly
     * If n is not suitable for square lattice, then last sites are left unfilled
     * Lattice is centered between (0,0) and (Lx, Ly).
     * The final argument should be passed one of the class variables VERTICAL or HORIZONTAL, indicating
     *   whether successive points fill the lattice across or down.
     */
    public static Space2D.Vector[] squareLattice(int n, double Lx, double Ly, boolean fillVertical) {
        Space2D.Vector[] r = new Space2D.Vector[n];
        for(int i=0; i<n; i++) {r[i] = new Space2D.Vector();}

        int moleculeColumns, moleculeRows;
        double moleculeInitialSpacingX, moleculeInitialSpacingY;

    //  Number of molecules per row (moleculeColumns) and number of rows (moleculeRows)
    //  in initial configuration
        moleculeColumns = (int)Math.sqrt(Lx/Ly*(double)n);
        moleculeRows = (int)(n/moleculeColumns);

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
            r[i].x = ix*moleculeInitialSpacingX + moleculeColumnsShift;
	        r[i].y = iy*moleculeInitialSpacingY + moleculeRowsShift;
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
    }
}
