package etomica.lattice.crystal;
import etomica.*;
import etomica.lattice.Crystal;

/**
 * Cubic primitive with a 4-site fcc basis.
 * 
 */

 /* History
  * 09/23/02 (DAK) new
  * 12/09/02 (skkwak) built another contructor to pass an instance of
  * AtomFactory class into BasisCubicFcc to change atoms to sites (concrete
  * sites under basis from BravisLattice)
  * 01/20/04 (DAK) restructured constructors
  */
public class CrystalFcc extends Crystal {

	/**
	 * Cubic fcc crystal in which each fcc site is populated by a mono atom
	 * @param space
	 */
	public CrystalFcc(Space space) {
		this(space, new AtomFactoryMono(space, AtomSequencerFactory.SIMPLE));
	}

	/**
	 * Cubic fcc crystal in which each bcc site is populated by an atom made by
	 * the given factory
	 * @param space
	 * @param factory
	 */
	public CrystalFcc(Space space, AtomFactory factory) {
		this(new PrimitiveCubic(space), factory);
	}

	/**
	 * Auxiliary constructor needed to be able to pass new PrimitiveCubic and
	 * new BasisCubicFcc (which needs the new primitive) to super.
	 */	
	private CrystalFcc(PrimitiveCubic primitive, AtomFactory factory) {
		super(primitive, new BasisCubicFcc(primitive.space, factory, primitive));
		primitive.setSize(Math.sqrt(2.0)*Default.ATOM_SIZE);
	}
   
    public String toString() {return "Fcc";}
    
}//end of CrystalFcc