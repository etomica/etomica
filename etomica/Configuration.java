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
    
    public Configuration(Species s){
        species.addElement(s);
        initializeCoordinates();
    }
    
    public void add(Species s){
        if(s instanceof SpeciesWalls) {return;}
        species.addElement(s);
        initializeCoordinates();
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
        double[] momentumSum = new double[Space.D];
        int sum = 0;
        for(int j=0; j<species.size(); j++) {
            Species s = (Species)species.elementAt(j);
            double momentum = Math.sqrt(s.mass*temperature/Constants.KE2T*(double)Space.D);  //need to divide by sqrt(m) to get velocity
            for(Atom a=s.firstAtom(); a!=s.terminationAtom(); a=a.getNextAtom()) {
	            a.p[1] = Math.cos(2*Math.PI*rand.nextDouble());
    	        a.p[0] = Math.sqrt(1.0 - a.p[1]*a.p[1]);
                sum++;
	            double momentumNorm = Math.sqrt(Space.v1S(a.p));
	            Space.uTEa1(a.p, momentum/momentumNorm);
	            Space.uPEv1(momentumSum,a.p);
	        }
	    }
    //    Zero center-of-mass momentum
        Space.uDEa1(momentumSum,(double)sum);
        for(int j=0; j<species.size(); j++) {
            Species s = (Species)species.elementAt(j);
            for(Atom a=s.firstAtom(); a!=s.terminationAtom(); a=a.getNextAtom()) {
                Space.uMEv1(a.p, momentumSum);
            }
        }
    }
    
    public void initializeMomentum(Molecule m) {
        double momentum = Math.sqrt(m.getMass()*temperature/Constants.KE2T*(double)Space.D);
        for(Atom a=m.firstAtom(); a!=m.terminationAtom(); a=a.getNextAtom()) {
	        a.p[1] = Math.cos(2*Math.PI*rand.nextDouble());
            a.p[0] = Math.sqrt(1.0 - a.p[1]*a.p[1]);
	        double momentumNorm = Math.sqrt(Space.v1S(a.p));
	        Space.uTEa1(a.p, momentum/momentumNorm);
	    }
	}
        
    
    public abstract void initializeCoordinates();
    
    public final static boolean HORIZONTAL = false;
    public final static boolean VERTICAL = true;
    
    /**
     * Returns a set of n coordinates filling a square lattice of sides Lx and Ly
     * If n is not suitable for square lattice, that last sites are left unfilled
     * Lattice is centered between (0,0) and (Lx, Ly).
     * The first element of the returned array indicates the point (0...(n-1)) and the 
     *   second indicates the coordinate (x or y)
     * The final argument should be passed one of the class variables VERTICAL or HORIZONTAL, indicating
     *   whether successive points fill the lattice across or down.
     */
    public static double[][] squareLattice(int n, double Lx, double Ly, boolean fillVertical) {
        double[][] r = new double[n][Space.D];

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
            r[i][0] = ix*moleculeInitialSpacingX + moleculeColumnsShift;
	        r[i][1] = iy*moleculeInitialSpacingY + moleculeRowsShift;
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
