/*
 * Created on Dec 6, 2004
 *
 */
package etomica.models.hexane;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveTriclinic;
import etomica.math.geometry.Polytope;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Vector3D;

/**
 * Implements the triclinic crystal structure derived from Dr. Monson's data.
 * @author nancycribbin
 */

public class PrimitiveHexane extends Primitive {
    
    public PrimitiveHexane(Space space) {
        super(space, false);
        latticeVectors[0] = new Vector3D(-0.0450226624118149, 0.533363447960696, 2.05308369446092);
        latticeVectors[1] = new Vector3D(0.355734889084923, 1.281874771828840, 0.146059929673340);
        latticeVectors[2] = new Vector3D(0.763064234689232, -0.640942796742683, -0.083271977067129);
        
        setSize(new double[] {Math.sqrt(latticeVectors[0].squared()), 
                                Math.sqrt(latticeVectors[1].squared()),
                                Math.sqrt(latticeVectors[2].squared())});
        
        
        //In units of sigma
//        super(space, 2.86842545, 1.338313769, 1.290909603, 0.779541414, 1.109209538,
//                0.768992016);
        
        //In units of Angstroms
//        super(space, 28.68425449956680, 13.38313769270750, 12.90909602969220, 
//                0.779541414, 1.109209538, 0.768992016);
    }
    
    
    
    protected Primitive makeReciprocal() {
        return new PrimitiveHexane(space);
    }
    
    //called by update method of superclass
    protected void updateReciprocal() {
        
    }
    

    public void setSize(double[] newSize) {
        if (size[0] == newSize[0] && size[1] == newSize[1] && size[2] == newSize[2]) {
            // no change
            return;
        }
        super.setSize(newSize);
    }
    
    protected void update() {
        super.update();
    }
    
    
    /**
     * Returns a new, identical instance of this primitive.
     */
    public Primitive copy() {
        return new PrimitiveHexane(space);
    }
        
    public void scaleSize(double scale) {
        setSize(new double[]{size[0]*scale, size[1]*scale, size[2]*scale});
    }        
    
    public int[] latticeIndex(Vector q) {
        for(int i=0; i<D; i++) {
            double x = q.x(i)/size[i];
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
        }
        return idx;
    }

    public int[] latticeIndex(Vector q, int[] dimensions) {
        for(int i=0; i<D; i++) {
            double x = q.x(i)/size[i];
            idx[i] = (x < 0) ? (int)x - 1 : (int)x; //we want idx to be the floor of x
            while(idx[i] >= dimensions[i]) idx[i] -= dimensions[i];
            while(idx[i] < 0)              idx[i] += dimensions[i];
        }
        return idx;
    }
    
    public Polytope wignerSeitzCell() {
        throw new RuntimeException("method PrimitiveOrthorhombic.wignerSeitzCell not yet implemented");
    }
    
    public Polytope unitCell() {
        throw new RuntimeException("method PrimitiveOrthorhombic.unitCell not yet implemented");
    }
    
    public String toString() {return "Hexane";}

}

