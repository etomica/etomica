package etomica.lattice.crystal;

import etomica.api.IVector;
import etomica.math.geometry.Polytope;
import etomica.space.Space;
import etomica.space3d.IVector3D;

/**
 * Primitive defined simply by its vectors.  Cannot be modified by setting
 * size or angle.
 *
 * @author Andrew Schultz
 */
public class PrimitiveGeneral extends Primitive {

    public PrimitiveGeneral(Space space, IVector[] primitiveVectors) {
        super(space);
        for (int i=0; i<latticeVectors.length; i++) {
            latticeVectors[i].E(primitiveVectors[i]);
        }
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
            IVector3D aStar = (IVector3D)space.makeVector();
            IVector3D bStar = (IVector3D)space.makeVector();
            IVector3D cStar = (IVector3D)space.makeVector();
            aStar.E(latticeVectors[1]);
            aStar.XE((IVector3D)latticeVectors[2]);
            double factor = 2.0*Math.PI/latticeVectors[0].dot(aStar); // a . (b X c)
            aStar.TE(factor);
            bStar.E(latticeVectors[2]);
            bStar.XE((IVector3D)latticeVectors[0]);
            bStar.TE(factor);
            cStar.E(latticeVectors[0]);
            cStar.XE((IVector3D)latticeVectors[1]);
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
