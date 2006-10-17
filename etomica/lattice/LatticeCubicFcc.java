package etomica.lattice;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.lattice.crystal.PrimitiveFcc;
import etomica.space3d.Space3D;

/**
 * Cubic primitive with a 4-site fcc basis.
 * 
 */


public class LatticeCubicFcc extends LatticeCrystal implements CubicLattice {

    /**
     * Cubic fcc crystal with a lattice constant that gives a
     * maximum-density structure for spheres of unit size. 
     * <p>
     * Use scaleBy method if desired to make lattice constant give
     * maximum density for another sphere size.
     */
    public LatticeCubicFcc() {
        this(Math.sqrt(2.0));
    }
    
    public LatticeCubicFcc(double latticeConstant) {
        this(new PrimitiveCubic(Space3D.getInstance()));
        primitive = (PrimitiveCubic)crystal.getLattice().getPrimitive();
        primitive.setCubicSize(latticeConstant);
    }

    /**
     * Auxiliary constructor needed to be able to pass new PrimitiveCubic and
     * new BasisCubicFcc (which needs the new primitive) to super.
     */ 
    private LatticeCubicFcc(PrimitiveCubic primitive) {
        super(new Crystal(primitive, new BasisCubicFcc(primitive)));
    }
    
    /**
     * Returns the cubic primitive on which the 4-site basis groups are placed.
     * Set the lattice constant via primitive().setSize(value).
     */
    public PrimitiveCubic primitive() {
        return primitive;
    }
    
    /**
     * Returns a new PrimitiveFcc instance corresponding to the fcc lattice. 
     * Changes to the returned instance have no effect on the lattice.
     */
    public PrimitiveFcc getPrimitiveFcc() {
        PrimitiveFcc p = new PrimitiveFcc(Space3D.getInstance());
        p.setCubicSize(getLatticeConstant()/Math.sqrt(2.0));
        return p;
    }
    
    /**
     * The lattice constant is the size of the cubic primitive vectors
     * of the lattice underlying this crystal.
     */
    public void setLatticeConstant(double latticeConstant) {
        primitive.setCubicSize(latticeConstant);
    }
    
    public double getLatticeConstant() {
        return primitive.getCubicSize();
    }
    

    /**
     * Rescales the lattice by the given factor. Multiplies the lattice constant
     * by the given value.
     */
    public void scaleBy(double scaleFactor) {
        primitive.scaleSize(scaleFactor);
   }

    /**
     * Returns "Fcc".
     */
    public String toString() {return "Fcc";}
    
    private PrimitiveCubic primitive;
}//end of CrystalFcc