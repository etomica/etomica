package etomica.space1d;

import etomica.Default;
import etomica.Phase;
import etomica.Simulation;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
/**
 * Class for implementing no periodic boundary conditions
 */
public final class BoundaryNone extends Boundary{
    private final Vector temp = new Vector();
    private final Vector dimensions = new Vector(Default.BOX_SIZE);
    private final Vector dimensionsCopy = new Vector();
    public final etomica.space.Vector dimensions() {dimensionsCopy.E(dimensions); return dimensionsCopy;}
    public BoundaryNone() {super();}
    public BoundaryNone(Phase p) {super(p);}
    public etomica.space.Boundary.Type type() {return Boundary.NONE;}
    public void nearestImage(Vector dr) {}
    public void nearestImage(etomica.space.Vector dr) {}
    public boolean centralImage(Vector r) {return false;}
    public boolean centralImage(etomica.space.Vector r) {return false;}
    public boolean centralImage(Coordinate c) {return false;}
    public double volume() {return dimensions.x;}
    public void inflate(double s) {dimensions.TE(s);}
    public void inflate(etomica.space.Vector s) {dimensions.TE(s);}
    public void setDimensions(etomica.space.Vector v) {dimensions.E(v);}
    public double[][] imageOrigins(int nShells) {return new double[0][1];}
    public float[][] getOverflowShifts(etomica.space.Vector rr, double distance) {return shift0;}
    public etomica.space.Vector randomPosition() {  //arbitrary choice for this method in this boundary
        temp.x = dimensions.x*Simulation.random.nextDouble(); 
        return temp;
    }
}
