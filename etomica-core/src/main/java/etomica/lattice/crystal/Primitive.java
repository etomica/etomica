/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice.crystal;

import etomica.space.Vector;
import etomica.math.geometry.LineSegment;
import etomica.math.geometry.Parallelepiped;
import etomica.math.geometry.Parallelogram;
import etomica.math.geometry.Polytope;
import etomica.space.Space;

/**
 * Collection of primitive elements that specify or are determined
 * by the structure of a Bravais lattice.
 */
public abstract class Primitive implements java.io.Serializable {
    
    private static final long serialVersionUID = 1L;
    protected final Vector[] latticeVectors;
    protected final int[] idx;//used to return coordinate index
    protected final int D;
    protected final double[] size;
    protected final double[] angle;
    protected final Space space;
    protected static final double rightAngle = 0.5*Math.PI;
    
    /**
     * This constructor is called directly when a Primitive is constructing
     * its reciprocal primitive.  For construction of the direct-lattice
     * primitive, this constructor is called via the Primitive(Simulation) constructor.
     */
    public Primitive(Space space) {
        this.space = space;
        D = space.D();
        latticeVectors = new Vector[D];
        idx = new int[D];
        size = new double[D];
//        sizeCopy = new double[D];
        angle = new double[D];
        for(int i=0; i<D; i++) {
            latticeVectors[i] = space.makeVector();
            angle[i] = rightAngle;
        }
    }
    
    /**
     * Method defining and constructing reciprocal primitive;
     */
    public abstract Primitive makeReciprocal();
    
    /**
     * @return the space
     */
    public final Space getSpace() {
        return space;
    }

    /**
     * Returns a new array with values equal to the lengths of the primitive vectors.
     */
    public double[] getSize() {
        return size.clone();
    }
    
    /**
     * Sets the length of each primitive vector to the corresponding
     * value in the given array.
     */
    protected final void setSize(double[] newSize) {
        for (int i=0; i<newSize.length; i++) {
            if (newSize[i] <= 0.0) { 
                throw new IllegalArgumentException("sizes must be positive");
            }
        }
        for (int i=0; i<newSize.length; i++) {
            size[i] = newSize[i];
        }
        update();
    }
    
    /**
     * Sets the angles between the primitive vector to the corresponding
     * values in the given array.
     */
    protected final void setAngles(double[] newAngle) {
        for (int i=0; i<newAngle.length; i++) {
            if (newAngle[i] < 0 || newAngle[i] > Math.PI) {
                throw new IllegalArgumentException("Angles must be between 0 and pi");
            }
        }
        for (int i=0; i<newAngle.length; i++) {
            angle[i] = newAngle[i];
        }
        update();
    }

    protected abstract void update();
        
    /**
     * Scales (multiplies) the size of each primitive vector by the given value.
     */
    public abstract void scaleSize(double scale);
        
    /**
     * Returns the primitive vectors.  Does not return the original
     * vectors used by the class, but instead returns copies.  Thus
     * changing the vectors returned by this method does not modify
     * the primitive vectors used by this class.  Subclasses should
     * provide mutator methods that permit changes to the vectors while
     * adhering to a particular structure (e.g., cubic, fcc, etc.).
     */
    public Vector[] vectors() {
        return latticeVectors;
    }
    
    /**
     * Returns a new, identical instance of this primitive.
     */
    public abstract Primitive copy();
    
    /**
     * Returns the index which would give the unit cell containing the given
     * point if the index were passed to a the site method of a sufficiently
     * large lattice that uses this primitive.
     */
    public abstract int[] latticeIndex(Vector r);
    
    /**
     * Same as latticeIndex(Space.Vector), but gives index for periodic system
     * with number of unit cells in each direction as given by the dimensions array.
     * If lattice index corresponds to a cell outside the range of dimensions,
     * index of image in central cells is returned.
     */
    public abstract int[] latticeIndex(Vector r, int[] dimensions);
    
    /**
     * Returns the Wigner-Seitz cell specified by this primitive.
     * The returned cell does not remain tied to the primitive, and
     * will not be updated with changes to the primitive.
     */
    public abstract Polytope wignerSeitzCell();
    
    /**
     * Returns a the unit cell specified by this primitive.
     * The returned cell does not remain tied to the primitive, and
     * will not be updated with changes to the primitive.
     */
    public Polytope unitCell() {
        if (space.D() == 1) {
            LineSegment line = new LineSegment(space);
            line.setLength(latticeVectors[0].getX(0));
            return line;
        }
        if (space.D() == 2) {
            return new Parallelogram(space, latticeVectors[0], latticeVectors[1]);
        }
        if (space.D() == 3) {
            return new Parallelepiped(space, latticeVectors[0], latticeVectors[1], latticeVectors[2]);
        }
        throw new RuntimeException("I'm impressed by your ability to make a "+D+"-D space, but I really don't know how to make an appropriate unit cell");
    }
}
