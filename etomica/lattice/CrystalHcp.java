package etomica.lattice;
import etomica.*;

/**
 * Hexagonal primitive with a 2-site hcp basis.
 */

 /* History
  * 09/27/02 (DAK) new
  * 01/20/04 (DAK) revised constructors, added one to take AtomFactory argument
  */
public class CrystalHcp extends Crystal {
    
	/**
	 * Hcp crystal in which each hcp site is populated by a mono atom
	 * @param space
	 */
	public CrystalHcp(Space space) {
		this(space, new AtomFactoryMono(space, AtomSequencerFactory.SIMPLE));
	}

	/**
	 * Hcp crystal in which each hcp site is populated by an atom made by the
	 * given factory
	 * @param space
	 * @param factory
	 */
	public CrystalHcp(Space space, AtomFactory factory) {
		this(new PrimitiveHexagonal(space), factory);
	}

	/**
	 * Auxiliary constructor needed to be able to pass new PrimitiveHexagonal
	 * and new BasisHcp (which needs the new primitive) to super.
	 */	
	private CrystalHcp(PrimitiveHexagonal primitive, AtomFactory factory) {
		super(primitive, new BasisHcp(primitive.space, factory, primitive));
		((PrimitiveHexagonal)primitive).setA(Default.ATOM_SIZE);
		((PrimitiveHexagonal)primitive).setC(Math.sqrt(8.0/3.0)*Default.ATOM_SIZE);
	}
    
    public String toString() {return "Hcp";}
    
}//end of CrystalHcp