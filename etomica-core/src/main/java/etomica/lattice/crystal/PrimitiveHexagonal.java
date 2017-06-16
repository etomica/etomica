/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice.crystal;
import etomica.space.Vector;
import etomica.math.geometry.Polytope;
import etomica.space.Space;
import etomica.space3d.Space3D;

/**
 * Primitive group for a hexagonal system.  Primitive-vector angles
 * are (90,90,120) degrees and two vectors are of equal length.
 */
public class PrimitiveHexagonal extends Primitive {
    
    private static final long serialVersionUID = 1L;
    protected double ab;
    protected static final double gamma = etomica.units.Degree.UNIT.toSim(120.);
    protected static final double cosGamma = Math.cos(gamma);
    protected static final double sinGamma = Math.sin(gamma);
    
    public PrimitiveHexagonal(Space space) {
        this(space, 1.0, 1.0);
    }
    public PrimitiveHexagonal(Space space, double ab, double c) {
        super(space);
        setSize(new double[]{ab, ab, c});
        this.ab = ab;
        setAngles(new double[]{rightAngle, rightAngle, gamma});
    }
    
    //called by superclass constructor
    public Primitive makeReciprocal() {
//        return PrimitiveHexagonalReciprocal(space, 2.0*Math.PI/(ab*sinGamma),
//                2.0*Math.PI/(ab*sinGamma), 2.0*Math.PI/size[2], 2.0*Math.PI/(size[2]*Math.sin(angle[1])), angle[1]);
        return new PrimitiveHexagonalReciprocal(space, 2.0 * Math.PI / (ab * sinGamma), 2.0 * Math.PI / size[2]);
    }
    
    public void setSizeAB(double newAB) {
        if (newAB == ab) {
            return;
        }
        setSize(new double[]{newAB, newAB, size[2]});
        ab = newAB;
    }
    
    public double getSizeAB() {
        return ab;
    }
    
    public void setSizeC(double newC) {
        if (newC == size[2]) {
            return;
        }
        setSize(new double[]{ab, ab, newC});
    }
    public double getSizeC() {return size[2];}

    public double getAngleAlpha() {
        return angle[0];
    }

    public double getAngleBeta() {
        return angle[1];
    }

    public double getAngleGamma() {
        return angle[2];
    }

    protected void update() {
        latticeVectors[0].setX(0,size[0]);
        latticeVectors[1].setX(0,size[0]*cosGamma);
        latticeVectors[1].setX(1,size[0]*sinGamma);
        latticeVectors[2].setX(2,size[2]);
    }
    
    /**
     * Returns a new PrimitiveHexagonal with the same size as this one.
     */
    public Primitive copy() {
        return new PrimitiveHexagonal(space, ab, size[2]);
    }
    
    public void scaleSize(double scale) {
        setSize(new double[]{scale*ab, scale*ab, scale*size[2]});
        ab = ab * scale;
    }

    public int[] latticeIndex(Vector q) {
        throw new RuntimeException("latticeIndex method not implemented yet in primitive");
   /*     for(int i=0; i<D; i++) {
            double x = q.x(i)/size;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
        }
        return idx;
   */ }
    
    public int[] latticeIndex(Vector q, int[] dimensions) {
        throw new RuntimeException("latticeIndex method not implemented yet in primitive");
   /*     for(int i=0; i<D; i++) {
            double x = q.x(i)/size;
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
            while(idx[i] >= dimensions[i]) idx[i] -= dimensions[i];
            while(idx[i] < 0)              idx[i] += dimensions[i];
        }
        return idx;
    */}
    
    public Polytope wignerSeitzCell() {
        throw new RuntimeException("method wignerSeitzCell not yet implemented");
    }
    
    public String toString() {return "Hexagonal";}
    
    protected static class PrimitiveHexagonalReciprocal extends PrimitiveHexagonal {
        public PrimitiveHexagonalReciprocal(Space space, double ab, double c) {
            super(space, ab, c);
        }

        protected void update() {
            latticeVectors[0].setX(0, ab*sinGamma);
            latticeVectors[0].setX(1,-ab*cosGamma);
            latticeVectors[1].setX(1, ab);
            latticeVectors[2].setX(2, size[2]);
        }
        
        private static final long serialVersionUID = 1L;
    }

    public static void main(String args[]) {
        PrimitiveHexagonal primitive = new PrimitiveHexagonal(Space3D.getInstance(), 2, 2);
        Vector[] v = primitive.vectors();
        Primitive reciprocal = primitive.makeReciprocal();
        Vector[] vr = reciprocal.vectors();
        for (int j=0; j<vr.length; j++) {
            System.out.println(j+" "+v[j]+" "+vr[j]);
        }
        for (int i=0; i<v.length; i++) {
            for (int j=0; j<vr.length; j++) {
                System.out.println(i+" "+j+" "+v[i].dot(vr[j]));
            }
        }
    }
}
    
