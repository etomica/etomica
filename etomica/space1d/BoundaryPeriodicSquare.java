package etomica.space1d;

import etomica.Default;
import etomica.Phase;
import etomica.Simulation;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
/**
 * Class for implementing rectangular periodic boundary conditions
 */
public class BoundaryPeriodicSquare extends Boundary implements Boundary.Periodic{
    public BoundaryPeriodicSquare() {this(Default.BOX_SIZE);}
    public BoundaryPeriodicSquare(Phase p) {this(p,Default.BOX_SIZE);}
    public BoundaryPeriodicSquare(Phase p, double lx) {super(p); dimensions.x = lx; updateDimensions();}
    public BoundaryPeriodicSquare(double lx) {dimensions.x = lx; updateDimensions();}
    public etomica.space.Boundary.Type type() {return Boundary.PERIODIC_SQUARE;}
    private final Vector temp = new Vector();
    private final Vector dimensions = new Vector();
    private final Vector dimensionsHalf = new Vector();
    private final Vector dimensionsCopy = new Vector();
    private final void updateDimensions() {
        dimensionsHalf.Ea1Tv1(0.5,dimensions);
        dimensionsCopy.E(dimensions);
    }
    public final etomica.space.Vector dimensions() {return dimensionsCopy;}
    public etomica.space.Vector randomPosition() {
        temp.x = dimensions.x*Simulation.random.nextDouble(); 
        return temp;
    }
    public void nearestImage(etomica.space.Vector dr) {nearestImage((Vector)dr);}
    public void nearestImage(Vector dr) {
        while(dr.x > +dimensionsHalf.x) dr.x -= dimensions.x;
        while(dr.x < -dimensionsHalf.x) dr.x += dimensions.x;
    }
    public boolean centralImage(Coordinate c) {return centralImage(c.position());}
    public boolean centralImage(etomica.space.Vector r) {return centralImage((Vector)r);}
    public boolean centralImage(Vector r) {
        boolean changed = false;
        while(r.x > dimensions.x) {r.x -= dimensions.x; changed = true;}
        while(r.x < 0.0)          {r.x += dimensions.x; changed = true;}
        return changed;
    }
    public void inflate(double scale) {
        dimensions.TE(scale); 
        updateDimensions();
        phase().boundaryEventManager.fireEvent(inflateEvent.setScale(scale));
    }
    public void inflate(etomica.space.Vector scale) {
        dimensions.TE(scale); 
        updateDimensions();
        phase().boundaryEventManager.fireEvent(inflateEvent.setScale(scale));
    }
    public void setDimensions(etomica.space.Vector v) {dimensions.E(v); updateDimensions();}
    public double volume() {return dimensions.x;}
    /** Computes origins for periodic images
     */
    public double[][] imageOrigins(int nShells) {
        int nImages = 2*nShells;
        double[][] origins = new double[nImages][1];
        int k = 0;
        for(int i=-nShells; i<=nShells; i++) {
            if(i==0) {continue;}
            origins[k][0] = i*dimensions.x;
            k++;
        }
        return origins;
    }

    /** Returns coordinate shifts needed to draw all images that overflow into central image
     * 0, or 1 shifts may be returned
     */
    int shiftX;
    Vector r;
    public float[][] getOverflowShifts(etomica.space.Vector rr, double distance) {
        r = (Vector)rr;
        shiftX = 0;
        if(r.x-distance < 0.0) {shiftX = +1;}
        else if(r.x+distance > dimensions.x) {shiftX = -1;}
        
        if(shiftX == 0)
            return shift0;
        shift = new float[1][2];
        shift[0][0] = (float)(shiftX*dimensions.x);
        return shift;
    } //end of getOverflowShifts
}
