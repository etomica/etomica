/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice.crystal;
import etomica.space.Vector;
import etomica.math.geometry.Polytope;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;

/**
 * Primitive group for a face-centered-cubic system.
 */
public class PrimitiveFcc extends Primitive {
    
    private static final long serialVersionUID = 1L;
    //primitive vectors are stored internally at unit length.  When requested
    //from the vectors() method, copies are scaled to size and returned.
    //default size is 1.0
    private double cubicSize;
    private Vector[] unitVectors;
    private static final double FCC_ANGLE = Math.acos(0.5);
    
    public PrimitiveFcc(Space space) {
        this(space, 1.0);
    }
    public PrimitiveFcc(Space space, double size) {
        super(space);
        //set up orthogonal vectors of unit size
        unitVectors = new Vector[D];
        for(int i=0; i<D; i++) {
            unitVectors[i] = space.makeVector();
            unitVectors[i].E(1.0/Math.sqrt(2.0));
            unitVectors[i].setX(i,0.0);
        }
        setCubicSize(size); //also sets reciprocal via update
        double[] newAngles = new double[D];
        for (int i=0; i<D; i++) {
            newAngles[i] = FCC_ANGLE;
        }
        setAngles(newAngles);
    }
    
    public Primitive makeReciprocal() {
        return new PrimitiveBcc(space, Math.sqrt(6)*Math.PI/cubicSize);
    }
    
    /**
     * Returns a new PrimitiveCubic with the same size as this one.
     */
    public Primitive copy() {
        return new PrimitiveFcc(space, cubicSize);
    }
    
    /**
     * Sets the length of all primitive vectors to the given value.
     */
    public void setCubicSize(double newCubicSize) {
        if (newCubicSize == cubicSize) {
            // no change
            return;
        }
        double[] sizeArray = new double[D];
        for(int i=0; i<D; i++) {
            sizeArray[i] = newCubicSize;
        }
        setSize(sizeArray);
        cubicSize = newCubicSize;
    }

    protected void update() {
        for(int i=0; i<D; i++) latticeVectors[i].Ea1Tv1(size[0],unitVectors[i]);
    }
    
    /**
     * Returns the common length of all primitive vectors.
     */
    public double getCubicSize() {return cubicSize;}
    
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
    
    public String toString() {return "Fcc";}

    
    public static void main(String[] args) {
        PrimitiveFcc primitive = new PrimitiveFcc(Space3D.getInstance());
        PrimitiveBcc reciprocal = (PrimitiveBcc)primitive.makeReciprocal();
        Vector[] latticeVectors = primitive.vectors();
        Vector a = latticeVectors[0];
        Vector b = latticeVectors[1];
        Vector c = latticeVectors[2];
        Vector ar = reciprocal.vectors()[0];
        Vector br = reciprocal.vectors()[1];
        Vector cr = reciprocal.vectors()[2];
        System.out.println("Primitive");
        System.out.println(a);
        System.out.println(b);
        System.out.println(c);
        System.out.println("Reciprocal");
        //ar.TE(0.5); br.TE(0.5); cr.TE(0.5);
        System.out.println(ar);
        System.out.println(br);
        System.out.println(cr);
        Vector3D work = new Vector3D();
        work.E(b);
        work.XE(c);
        double norm = a.dot(work);
        work.TE(2*Math.PI/norm);
        work.ME(ar);
        System.out.println("check");
        System.out.println(work);
        work.E(c);
        work.XE(a);
        work.TE(2*Math.PI/norm);
        work.ME(br);
        System.out.println(work);
        work.E(a);
        work.XE(b);
        work.TE(2*Math.PI/norm);
        work.ME(cr);
        System.out.println(work);
    }
}
