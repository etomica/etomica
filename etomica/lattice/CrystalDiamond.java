package etomica.lattice;
import etomica.*;

/**
 * Cubic primitive with a 4-site fcc basis, on which each site 
 * is a 2-site diamond basis.
 */

 /* History
  * 09/26/02 (DAK) new
  * 01/20/04 (DAK) revised constructors; added one taking atomFactory argument
  */
public class CrystalDiamond extends Crystal {
    
	/**
	 * Cubic diamond crystal in which each diamond site is populated by an atom
	 * made by the given factory
	 * @param space
	 */
	public CrystalDiamond(Space space) {
		this(space, new AtomFactoryMono(space, AtomSequencerSimple.FACTORY));
	}

	/**
	 * Cubic diamond crystal in which each diamond site is populated by an atom
	 * made by the given factory
	 * @param space
	 * @param factory
	 */
	public CrystalDiamond(Space space, AtomFactory factory) {
		this(new PrimitiveCubic(space), factory);
	}

	/**
	 * Auxiliary constructor needed to be able to pass new PrimitiveCubic and
	 * new BasisCubicDiamond (which needs the new primitive) to super.
	 */	
	private CrystalDiamond(PrimitiveCubic primitive, AtomFactory factory) {
		super(primitive, new BasisCubicDiamond(primitive.space, factory, primitive));
		((PrimitiveCubic)primitive).setSize(4.0/Math.sqrt(3.0)*Default.ATOM_SIZE);
	}
   
    public String toString() {return "Diamond";}
    
}//end of CrystalDiamond