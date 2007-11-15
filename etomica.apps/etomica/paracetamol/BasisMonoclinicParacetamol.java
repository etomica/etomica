package etomica.paracetamol;
import etomica.lattice.crystal.Basis;
import etomica.space3d.Vector3D;

/**
 * An 4-atom basis that for a monoclinic paracetamol crystal.  This Form I crystal structure 
 * occurs at T = 10K. 
 *
 * @author Tai Tan
 */
public class BasisMonoclinicParacetamol extends Basis {
    
    /**
     * Makes a monoclinic-crystal paracetamol 4-atom basis.
     */
    public BasisMonoclinicParacetamol() {
        super(scaledPositions); 
    }
  
    
    private static final Vector3D[] scaledPositions = new Vector3D[] {

    	new Vector3D(  0.929048,    0.355078,   0.900198),
    	new Vector3D(1-0.429048,	0.855078, 1-0.400198),
    	new Vector3D(1-0.929048,  1-0.355078, 1-0.900198),
    	new Vector3D(  1.429048-1,	0.144922,	1.400198-1)

    };

    private static final long serialVersionUID = 1L;
}