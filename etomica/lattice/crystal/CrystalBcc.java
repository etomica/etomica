package etomica.lattice.crystal;
import etomica.*;
import etomica.lattice.Crystal;

/**
 * Cubic primitive with a 2-site bcc basis.
 */

 /* History
  * 09/23/02 (DAK) new
  * 01/19/04 (DAK) modified with use of Basis class
  * 01/20/04 (DAK) added new constructors to permit specification of factory
  */
public class CrystalBcc extends Crystal {
    
    /**
	 * Cubic bcc crystal in which each bcc site is populated by a mono atom
     * @param space
     */
    public CrystalBcc(Space space) {
        this(space, new AtomFactoryMono(space, AtomSequencerFactory.SIMPLE));
    }

	/**
	 * Cubic bcc crystal in which each bcc site is populated by an atom made by
	 * the given factory
	 * @param space
	 * @param factory
	 */
	public CrystalBcc(Space space, AtomFactory factory) {
		this(new PrimitiveCubic(space), factory);
	}

	/**
	 * Auxiliary constructor needed to be able to pass new PrimitiveCubic and
	 * new BasisCubicBcc (which needs the new primitive) to super.
	 */	
	private CrystalBcc(PrimitiveCubic primitive, AtomFactory factory) {
		super(primitive, new BasisCubicBcc(primitive.space, factory, primitive));
		primitive.setSize(2.0/Math.sqrt(3.0)*Default.ATOM_SIZE);
	}
    
    public String toString() {return "Bcc";}
    
}//end of CrystalBcc