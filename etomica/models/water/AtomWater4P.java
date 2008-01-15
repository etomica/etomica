package etomica.models.water;

import etomica.atom.Molecule;
import etomica.atom.AtomType;
import etomica.atom.IAtomPositioned;

/**
 * Tree node for a 4-point water molecule.
 */
public class AtomWater4P extends Molecule {

    /**
	 * Constructor for AtomTreeNodeWater.
	 * @param atom
	 * @param parent
	 */
	public AtomWater4P(AtomType waterType) {
		super(waterType);
	}
	
	public final static int indexH1 = 0;
	public final static int indexH2 = 1;
	public final static int indexO  = 2;
	public final static int indexM  = 3;
	
    private static final long serialVersionUID = 1L;
	public IAtomPositioned H1, H2, O, M;
	
	public final static double [] Echarge = new double [4];

}
