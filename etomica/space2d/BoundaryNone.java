package etomica.space2d;

import etomica.Default;
import etomica.Phase;
import etomica.Simulation;
import etomica.Space2D;
import etomica.space.Boundary;
import etomica.space.Boundary.Type;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
/**
 * Class for implementing no periodic boundary conditions
 */
class BoundaryNone extends Boundary {
    private final Vector temp = new Vector();
    private final Vector dimensions = new Vector(Default.BOX_SIZE, Default.BOX_SIZE);
    private final Vector dimensionsCopy = new Vector();
    public Vector dimensions() {dimensionsCopy.E(dimensions); return dimensionsCopy;}
    public BoundaryNone() {super();}
    public BoundaryNone(Phase p) {super(p);}
    public Boundary.Type type() {return Boundary.NONE;}
    public final void nearestImage(Vector dr) {}
    public final boolean centralImage(Vector r) {return false;}
    public final void nearestImage(Vector dr) {}
    public final boolean centralImage(Vector r) {return false;}
    public final boolean centralImage(Coordinate c) {return false;}
    public double volume() {return dimensions.x*dimensions.y;}
    public void inflate(double s) {dimensions.TE(s);}
    public void inflate(Vector s) {dimensions.TE(s);}
    public void setDimensions(Vector v) {dimensions.E(v);}
    public double[][] imageOrigins(int nShells) {return new double[0][Space2D.D];}
    public float[][] getOverflowShifts(Vector rr, double distance) {return shift0;}
    public Vector randomPosition() {  //arbitrary choice for this method in this boundary
        temp.x = dimensions.x*Simulation.random.nextDouble(); 
        temp.y = dimensions.y*Simulation.random.nextDouble(); 
        return temp;
    }
}