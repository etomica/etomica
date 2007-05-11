package etomica.paracetamol;
import etomica.lattice.crystal.Basis;
import etomica.space3d.Vector3D;

/**
 * An 8-atom basis that for a orthorhombic paracetamol crystal.  This Form II crystal structure 
 * occurs at T = 123K. 
 *
 * @author Tai Tan
 */
public class BasisOrthorhombicParacetamol extends Basis {
    
    /**
     * Makes a orthorhombic-crystal paracetamol 8-atom basis.
     */
    public BasisOrthorhombicParacetamol() {
        super(scaledPositions); 
    }
  
    
    private static final Vector3D[] scaledPositions = new Vector3D[] {

    	new Vector3D(  0.162154,   0.936073,	0.282665),
    	new Vector3D(  0.337846, 1-0.936073,	0.782665),
    	new Vector3D(1-0.162154,   1.436073-1,	0.217335),
    	new Vector3D(  0.662154, 1-0.436073,  1-0.282665),
    	new Vector3D(1-0.162154, 1-0.936073,  1-0.282665),
    	new Vector3D(  0.662154,   0.936073,	0.217335),
    	new Vector3D(  0.162154, 1-0.436073,	0.782665),
    	new Vector3D(  0.337846,   1.436073-1,	0.282665)


    };

    private static final long serialVersionUID = 1L;
}