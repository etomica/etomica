package etomica.space;

import etomica.Constants;
import etomica.NearestImageTransformer;
import etomica.Space;

/*
 * History Created on Jan 24, 2005 by kofke
 */
public abstract class Boundary implements NearestImageTransformer,
        java.io.Serializable {

    protected final float[][] shift0 = new float[0][0];//cannot be static
                                                       // because several phases
                                                       // may be using at once
    protected float[][] shift;
    protected Space space;

    public Boundary(Space space) {
        this.space = space;
    }

    public abstract Vector centralImage(Vector r);//returns true if r is
                                                  // changed by applying central
                                                  // image

    public abstract void nearestImage(Vector dr);

    public abstract double volume();

    /**
     * Returns a copy of the dimensions, as a Vector. Manipulation of this copy
     * will not cause any change to the boundary's dimensions.
     */
    public abstract Vector dimensions();

    public abstract void setDimensions(Vector v);

    public abstract Vector randomPosition();

    public abstract float[][] getOverflowShifts(Vector r, double distance);

    /**
     * Set of vectors describing the displacements needed to translate the
     * central image to all of the periodic images. Returns a two dimensional
     * array of doubles. The first index specifies each perioidic image, while
     * the second index indicates the x and y components of the translation
     * vector.
     * 
     * @param nShells
     *            the number of shells of images to be computed
     */
    public abstract double[][] imageOrigins(int nShells);

    public static abstract class Type extends Constants.TypedConstant {

        protected Type(String label) {
            super(label);
        }
    }

    /**
     * Placeholder boundary that performs no actions and returns null or zero
     * from every method.
     */
    public static final Boundary NULL = new Boundary(null) {

        public Vector dimensions() {
            return null;
        }

        public Boundary.Type type() {
            return null;
        }

        public final void nearestImage(Vector dr) {
        }

        public final Vector centralImage(Vector r) {
            return null;
        }

        public double volume() {
            return 0.0;
        }

        public void setDimensions(Vector v) {
        }

        public double[][] imageOrigins(int nShells) {
            return null;
        }

        public float[][] getOverflowShifts(Vector r, double distance) {
            return null;
        }

        public Vector randomPosition() {
            return null;
        }
    };//end of NULL

}