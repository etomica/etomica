package etomica.space1d;

import etomica.util.random.IRandom;
import etomica.space.Vector;
import etomica.space.IOrientation;

/**
 *
 * Orientation in a 1-dimensional space. Only two values are possible, indicating
 * one of two directions. Value of orientation is given as +1 or -1.
 *
 * @author David Kofke.
 */
public class Orientation1D implements IOrientation {

    private int direction;

    /**
     * Construct with default orientation of +1.
     */
    public Orientation1D() {
        direction = +1;
    }

    public void E(IOrientation o) {
        this.direction = ((Orientation1D)o).direction;
    }


    public Vector getDirection() {
        return new Vector1D(direction);
    }


    public void setDirection(Vector newDirection) {
        double x = ((Vector1D)newDirection).x;
        if(x == 0) throw new IllegalArgumentException("Can't set direction to zero.");
        direction = x > 0 ? +1 : -1;
    }

    /**
     * Sets orientation randomly to be +1 or -1.
     * @param random used to generate random value to set orientation
     * @param theta is ignored
     */
    public void randomRotation(IRandom random, double theta) {
        direction = (random.nextDouble() < 0.5) ? -1 : +1;
    }
}
