package etomica.space3d;

import etomica.Default;
import etomica.Phase;
import etomica.Simulation;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public final class BoundaryNone extends Boundary {
    private final Vector temp = new Vector();
    public final Vector dimensions = new Vector(Default.BOX_SIZE, Default.BOX_SIZE, Default.BOX_SIZE);
    public final etomica.space.Vector dimensions() {return dimensions;}
    public BoundaryNone() {super();}
    public BoundaryNone(Phase p) {super(p);}
    public etomica.space.Boundary.Type type() {return Boundary.NONE;}
    public void nearestImage(Vector dr) {}
    public void nearestImage(etomica.space.Vector dr) {}
    public boolean centralImage(Vector r) {return false;}
    public boolean centralImage(etomica.space.Vector r) {return false;}
    public void nearestImage(Vector dr, Vector shift) {shift.E(0.0);}
    public boolean centralImage(Coordinate c) {return false;}
    public double volume() {return dimensions.x*dimensions.y*dimensions.z;}
    public void inflate(double s) {dimensions.TE(s);}
    public void inflate(etomica.space.Vector s) {dimensions.TE(s);}
    public void setDimensions(etomica.space.Vector v) {dimensions.E(v);}
    public double[][] imageOrigins(int nShells) {return new double[0][3];}
    public float[][] getOverflowShifts(etomica.space.Vector rr, double distance) {return shift0;}
    public etomica.space.Vector randomPosition() {
        temp.x = dimensions.x*Simulation.random.nextDouble();
        temp.y = dimensions.y*Simulation.random.nextDouble();
        temp.z = dimensions.z*Simulation.random.nextDouble();
        return temp;
    }
}
