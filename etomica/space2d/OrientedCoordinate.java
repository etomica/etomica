package etomica.space2d;

import etomica.Atom;
import etomica.Space3D;
import etomica.atom.AtomType;
import etomica.atom.AtomType.SphericalTop;
import etomica.space.Coordinate;
import etomica.space.Orientation;
import etomica.space.Coordinate.Angular;
import etomica.space3d.Vector;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public class OrientedCoordinate extends Coordinate implements Coordinate.Angular {
    private double L = 0.0; //magnitude of angular momentum
    private final Vector vector = new Vector();//used to return vector quantities (be sure to keep x and y components zero)
    private final double[] I;
    private final Orientation orientation = new Orientation();
    public OrientedCoordinate(Atom a) {
        super(a);
        I = ((AtomType.SphericalTop)((Atom)a).type).momentOfInertia();
    }
    public Vector angularMomentum() {vector.setX(2,L); return vector;}
    public Vector angularVelocity() {vector.setX(2,L/I[0]); return vector;}
    public void angularAccelerateBy(Vector t) {L += t.x(2);}
    public Orientation orientation() {return orientation;}
    public double kineticEnergy() {return super.kineticEnergy() + 0.5*L*L/I[0];}
    public void freeFlight(double t) {
        super.freeFlight(t);
        orientation.rotateBy(t*L/I[0]);//all elements of I equal for spherical top
    }
}