package etomica.lattice;
import etomica.*;

/**
 * Class representing a crystal basis.  Primary purpose is to hold an
 * AtomFactory instance that produces atoms(s) or molecule that is to be placed
 * at each site of a Bravais lattice.  These atoms may either (1) represent
 * multiple lattice sites that are positioned at each Bravais site (e.g., a
 * cubic fcc lattice is obtained by placing 4 basis atoms at each cubic Bravis-
 * lattice site); or (2), they may simply be atoms/molecules that reside at the
 * Bravais site. The number of lattice sites resulting at each Bravais site is
 * given by the size() method.  In case (1) this is the number of sites in an
 * Atom produced by the AtomFactory; in case (2), this is 1.
 *
 * @author David Kofke
 */
 
 /* History
  * 09/16/02 (DAK) new
  * 01/19/04 (DAK) updated for subclassing by previously-implemented BasisXXX
  * classes.  Previously they had extended AtomFactory or a subclass of
  * AtomFactory.  Now AtomFactory is given as a field of this class.
  */
 
public class Basis {

	private AtomFactory atomFactory;
	private Space space;
	private int size; 
	    
	/**
	 * A simple 1-atom basis.
	 */
	public Basis(Space space) {
		this(space, 1, new AtomFactoryMono(space, AtomSequencerSimple.FACTORY));
	}
	
	/**
	 * A homogeneous n-atom basis arranged in the given configuration.  Each
	 * site in the basis is constructed by AtomFactoryMono.
	 * @param space 
	 * @param size Number of atoms (n) in basis
	 * @param configuration Specifies arrangement of atoms
	 */
	public Basis(Space space, int size, Configuration configuration) {
		this(space, size, configuration, new AtomFactoryMono(space,AtomSequencerSimple.FACTORY));
	}

	/**
	 * A homogeneous n-site basis arranged in the given configuration, with each
	 * site in the basis constructed using the given AtomFactory.
	 * @param space 
	 * @param size Number of atoms (n) in basis
	 * @param configuration Specifies arrangement of atoms
	 * @param atomFactory Used to constructe each site in the basis
	 */
	public Basis(Space space, int size, Configuration configuration, AtomFactory atomFactory) {
		this(space, size, new AtomFactoryHomo(space, AtomSequencerSimple.FACTORY, 
										      atomFactory, size, BondInitializer.NULL, configuration));
	}
	
	/**
	 * General basis with  the given factory making each instance of the full
	 * basis site (e.g., the factory will make the molecule occupying each site
	 * of a molecular crystal).
	 * @param space Governing Space instance
	 * @param size Number of sites in the basis.  If not unity, indicates that
	 * factory is producing the lattice sites, such as the 4 sites used to form
	 * a fcc lattice on a cubic primitive.
	 * @param factory Atom factory that makes occupant of each site of the
	 * lattice
	 */
	public Basis(Space space, int size, AtomFactory factory) {
		if(size <= 0) throw new IllegalArgumentException("Cannot have basis with non-positive size");
		this.space = space;
		this.size = size;
		this.atomFactory = factory;
	}
	
	/**
	 * Number of atoms in the basis.
	 * @return int
	 */
	public int size() {return size;}
    
    /**
     * Returns that AtomFactory that makes each basis instance.
     * @return AtomFactory
     */
    public AtomFactory atomFactory() {return atomFactory;}

}//end of Basis