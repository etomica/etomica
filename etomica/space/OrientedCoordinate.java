package etomica.space;
import etomica.*;

public class OrientedCoordinate extends Coordinate implements Space.Coordinate.Angular {
    private double L = 0.0; //magnitude of angular momentum
    private final Space3D.Vector vector = new Space3D.Vector();//used to return vector quantities (be sure to keep x and y components zero)
    private final double[] I;
    private final Space.Orientation orientation;
    public OrientedCoordinate(Space space, Atom a) {
        super(space, a);
        orientation = space.makeOrientation();
        I = ((AtomType.SphericalTop)a.type).momentOfInertia();
    }
    public Space3D.Vector angularMomentum() {vector.setX(2,L); return vector;}
    public Space3D.Vector angularVelocity() {vector.setX(2,L/I[0]); return vector;}
    public void angularAccelerateBy(Space3D.Vector t) {L += t.x(2);}
    public Space.Orientation orientation() {return orientation;}
    public double kineticEnergy() {return super.kineticEnergy() + 0.5*L*L/I[0];}
    public void freeFlight(double t) {
        super.freeFlight(t);
        throw new RuntimeException("freeFlight method not implemented in OrientedCoordinate");
//fix this        orientation.rotateBy(t*L/I[0]);//all elements of I equal for spherical top
    }
}
