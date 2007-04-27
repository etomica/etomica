package etomica.paracetamol;
import etomica.lattice.crystal.Basis;
import etomica.space3d.Vector3D;

/**
 * An 4-atom basis that for a monoclinic paracetamol crystal.  This Form I crystal structure 
 * occurs at T = 20K. 
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

    	new Vector3D(  0.953265,    0.337679,   0.870361),
    	new Vector3D(1-0.453265,	0.837679, 1-0.370361),
    	new Vector3D(1-0.953265,  1-0.337679, 1-0.870361),
    	new Vector3D(  1.453265-1,	0.162321,	1.370361-1)

    };

    private static final long serialVersionUID = 1L;
}