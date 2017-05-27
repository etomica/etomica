/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice.crystal;
import etomica.space.Vector;
import etomica.math.geometry.Polytope;
import etomica.space.Space;

/**
 * Primitive group for a body-centered-cubic system.
 */
public class PrimitiveBcc extends Primitive {
    
    private static final long serialVersionUID = 1L;
    private double cubicSize;
    private Vector[] unitVectors;
    private static final double BCC_ANGLE = Math.acos(1.0/3.0);
    
    public PrimitiveBcc(Space space) {
        this(space, 1.0);
    }
    public PrimitiveBcc(Space space, double size) {
        super(space);
        //set up orthogonal vectors of unit size
        unitVectors = new Vector[D];
        for(int i=0; i<D; i++) {
            unitVectors[i] = space.makeVector();
            unitVectors[i].E(1.0/Math.sqrt(3.0));
            unitVectors[i].setX(i,-1.0/Math.sqrt(3.0));
        }
        setCubicSize(size); //also sets reciprocal via update
        setAngles(new double[]{BCC_ANGLE, BCC_ANGLE, BCC_ANGLE});
    }
    
    public Primitive makeReciprocal() {
        return new PrimitiveFcc(space, 4.0*Math.PI/cubicSize);
    }
    
    /**
     * Returns a new PrimitiveCubic with the same size as this one.
     */
    public Primitive copy() {
        return new PrimitiveBcc(space, cubicSize);
    }
    
    /**
     * Sets the length of all primitive vectors to the given value.
     */
    public void setCubicSize(double newCubicSize) {
        if (newCubicSize == cubicSize) {
            return;
        }
        if (newCubicSize <= 0) {
            throw new IllegalArgumentException("BCC size must be positive");
        }
        double[] sizeArray = new double[D];
        for(int i=0; i<D; i++) {
            sizeArray[i] = newCubicSize;
        }
        setSize(sizeArray);
        cubicSize = newCubicSize;
    }
    
    /**
     * Returns the common length of all primitive vectors.
     */
    public double getCubicSize() {return cubicSize;}
    
    protected void update() {
        for(int i=0; i<D; i++) latticeVectors[i].Ea1Tv1(size[0],unitVectors[i]);
    }

    /**
     * Multiplies the size of the current vectors by the given value.
     */
    public void scaleSize(double scale) {
        setCubicSize(scale*cubicSize);
    }

    public int[] latticeIndex(Vector q) {
        throw new RuntimeException("PrimitiveFcc.latticeIndex not yet implemented");
/*        for(int i=0; i<D; i++) {
            double x = q.x(i)/size;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
        }
        return idx;
*/    }
    
    public int[] latticeIndex(Vector q, int[] dimensions) {
        throw new RuntimeException("PrimitiveFcc.latticeIndex not yet implemented");
 /*       for(int i=0; i<D; i++) {
            double x = q.x(i)/size;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
            while(idx[i] >= dimensions[i]) idx[i] -= dimensions[i];
            while(idx[i] < 0)              idx[i] += dimensions[i];
        }
        return idx;
 */   }
    
    public Polytope wignerSeitzCell() {
        throw new RuntimeException("method PrimitiveFcc.wignerSeitzCell not yet implemented");
    }
    
    public String toString() {return "Bcc";}
    
}
