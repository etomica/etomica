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

    	new Vector3D(  0.172278,   0.909616,	0.275662),
    	new Vector3D(  0.327722, 1-0.909616,	0.775662),
    	new Vector3D(1-0.172278,   1.409616-1,	0.224338),
    	new Vector3D(  0.672278, 1-0.409616,  1-0.275662),
    	new Vector3D(1-0.172278, 1-0.909616,  1-0.275662),
    	new Vector3D(  0.672278,   0.909616,	0.224338),
    	new Vector3D(  0.172278, 1-0.409616,	0.775662),
    	new Vector3D(  0.327722,   1.409616-1,	0.275662)


    };

    private static final long serialVersionUID = 1L;
}