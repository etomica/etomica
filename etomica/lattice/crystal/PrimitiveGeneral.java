package etomica.lattice.crystal;

import etomica.api.IVector;
import etomica.math.geometry.Polytope;
import etomica.space.ISpace;

/**
 * Primitive defined simply by its vectors.  Cannot be modified by setting
 * size or angle.
 *
 * @author Andrew Schultz
 */
public class PrimitiveGeneral extends Primitive {

    public PrimitiveGeneral(ISpace space, IVector[] primitiveVectors) {
        super(space);
        for (int i=0; i<latticeVectors.length; i++) {
            latticeVectors[i].E(primitiveVectors[i]);
        }
        double[] mySize = new double[primitiveVectors.length];
        for (int i=0; i<mySize.length; i++) {
            mySize[i] = Math.sqrt(primitiveVectors[i].squared());
        }
        setSize(mySize);
        double[] myAngles = new double[primitiveVectors.length == 2 ? 1 : 3];
        myAngles[0] = Math.acos(latticeVectors[0].dot(latticeVectors[1])/
                Math.sqrt(latticeVectors[0].squared()*latticeVectors[1].squared()));
        if (myAngles.length == 3) {
            myAngles[1] = Math.acos(latticeVectors[1].dot(latticeVectors[2])/
                    Math.sqrt(latticeVectors[1].squared()*latticeVectors[2].squared()));
            myAngles[2] = Math.acos(latticeVectors[2].dot(latticeVectors[0])/
                    Math.sqrt(latticeVectors[2].squared()*latticeVectors[0].squared()));
        }
        setAngles(myAngles);
    }

    public Primitive copy() {
        return new PrimitiveGeneral(space, latticeVectors);
    }

    public int[] latticeIndex(IVector r) {
        throw new RuntimeException("nope");
    }

    public int[] latticeIndex(IVector r, int[] dimensions) {
        throw new RuntimeException("nope");
    }

    public Primitive makeReciprocal() {
        if (space.D() == 3) {
            IVector aStar = space.makeVector();
            IVector bStar = space.makeVector();
            IVector cStar = space.makeVector();
            aStar.E(latticeVectors[1]);
            aStar.XE(latticeVectors[2]);
            double factor = 2.0*Math.PI/latticeVectors[0].dot(aStar); // a . (b X c)
            aStar.TE(factor);
            bStar.E(latticeVectors[2]);
            bStar.XE(latticeVectors[0]);
            bStar.TE(factor);
            cStar.E(latticeVectors[0]);
            cStar.XE(latticeVectors[1]);
            cStar.TE(factor);
            return new PrimitiveGeneral(space, new IVector[]{aStar, bStar, cStar});
        }
        if (space.D() == 2) {
            IVector aStar = space.makeVector();
            IVector bStar = space.makeVector();
            aStar.setX(0, -latticeVectors[0].x(1));
            aStar.setX(1, latticeVectors[0].x(0));
            aStar.TE(2.0*Math.PI/aStar.dot(latticeVectors[1]));
            bStar.setX(0, -latticeVectors[1].x(1));
            bStar.setX(1, latticeVectors[1].x(0));
            bStar.TE(2.0*Math.PI/aStar.dot(latticeVectors[0]));
            return new PrimitiveGeneral(space, new IVector[]{aStar, bStar});
        }
        throw new RuntimeException("can't make a "+space.D()+"D reciprocal");
    }

    public void scaleSize(double scale) {
        for (int i=0; i<latticeVectors.length; i++) {
            latticeVectors[i].TE(scale);
        }
    }

    protected void update() {
    }

    public Polytope wignerSeitzCell() {
        throw new RuntimeException("nope");
    }

    private static final long serialVersionUID = 1L;
}
