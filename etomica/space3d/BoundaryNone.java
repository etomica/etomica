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
    public final Vector dimensions() {return dimensions;}
    public BoundaryNone() {super();}
    public BoundaryNone(Phase p) {super(p);}
    public Boundary.Type type() {return Boundary.NONE;}
    public void nearestImage(Vector dr) {}
    public boolean centralImage(Vector r) {return false;}
    public void nearestImage(Vector dr, Vector shift) {shift.E(0.0);}
    public boolean centralImage(Coordinate c) {return false;}
    public double volume() {return dimensions.x*dimensions.y*dimensions.z;}
    public void inflate(double s) {dimensions.TE(s);}
    public void inflate(Vector s) {dimensions.TE(s);}
    public void setDimensions(Vector v) {dimensions.E(v);}
    public double[][] imageOrigins(int nShells) {return new double[0][Space3D.D];}
    public float[][] getOverflowShifts(Vector rr, double distance) {return shift0;}
    public Vector randomPosition() {
        temp.x = dimensions.x*Simulation.random.nextDouble();
        temp.y = dimensions.y*Simulation.random.nextDouble();
        temp.z = dimensions.z*Simulation.random.nextDouble();
        return temp;
    }
}
