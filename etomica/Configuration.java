package simulate;
import java.util.*;
import java.awt.*;

/**
 * General class for assignment of coordinates to all molecules/atoms in phase
 */

public abstract class Configuration extends Component{

    protected Vector species = new Vector();
    public double temperature;
    Phase parentPhase;
    Random rand = new Random();

    public Configuration(){
        setTemperature(300.);
    }
    
    public Configuration(Species.Agent s){
        species.addElement(s);
        initializeCoordinates();
    }
    
    public void add(Species.Agent s){
 //       if(s.firstAtom().type instanceof AtomType.Wall) {return;}
        species.addElement(s);
        if(s.getNMolecules() > 0) initializeCoordinates();
    }

    public double getTemperature(){
        return temperature;
    }
    public void setTemperature(double t){
        temperature = t;
    }
 /**   
  * All atom velocities are set such that all have the same total momentum (corresponding to
  * the current value of temperature), with the direction at random
  */
    public void initializeMomenta() {
        Space.Vector momentumSum = parentPhase.parentSimulation.space.makeVector();
        int sum = 0;
        for(int j=0; j<species.size(); j++) {
            Species.Agent s = (Species.Agent)species.elementAt(j);
            for(Atom a=s.firstAtom(); a!=s.terminationAtom(); a=a.nextAtom()) {
                a.coordinate.randomizeMomentum(temperature);
                sum++;
                momentumSum.PE(a.coordinate.momentum());
	        }
	    }
    //    Zero center-of-mass momentum
        momentumSum.DE((double)sum);
        for(int j=0; j<species.size(); j++) {
            Species.Agent s = (Species.Agent)species.elementAt(j);
            for(Atom a=s.firstAtom(); a!=s.terminationAtom(); a=a.nextAtom()) {
                a.coordinate.momentum().ME(momentumSum);
            }
        }
    }
    
    public void initializeMomentum(Molecule m) {
        m.randomizeMomentum(temperature);
	}
        
    
    public abstract void initializeCoordinates();
    
    public final static boolean HORIZONTAL = false;
    public final static boolean VERTICAL = true;
    
    /**
     * Returns a set of n coordinates filling a square lattice of sides Lx and Ly
     * If n is not suitable for square lattice, then last sites are left unfilled
     * Lattice is centered between (0,0) and (Lx, Ly).
     * The final argument should be passed one of the class variables VERTICAL or HORIZONTAL, indicating
     *   whether successive points fill the lattice across or down.
     */
    public static Space2DCell.Vector[] squareLattice(int n, double Lx, double Ly, boolean fillVertical) {
        Space2DCell.Vector[] r = new Space2DCell.Vector[n];
        for(int i=0; i<n; i++) {r[i] = new Space2DCell.Vector();}

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
	//            if(ix >= moleculeColumns) {  delete if works ok
	            if(ix >= columnsDrawn) {
	                ix = 0;
	                iy++;
	            }
	        }
	    }
	    return r;
    }
}
