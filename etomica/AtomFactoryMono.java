package etomica;

/**
 * Builder of a monoatomic atom group, which comprises just an Atom.
 *
 * @author David Kofke
 */

/* History
 * 08/26/03 (DAK) removed override of build method, instead specifying
 * AtomTreeNodeLeaf.FACTORY to superclass in constructor
 * 
 */
public class AtomFactoryMono extends AtomFactory {
    
    /**
     * Constructor with neighborSequencerFactory and AtomType.Sphere defaults.
     */
    public AtomFactoryMono(Space space, AtomSequencerFactory seqFactory) {
        super(space, seqFactory, AtomTreeNodeLeaf.FACTORY);
        init();
    }
    
    private void init() {
        setType(new AtomTypeSphere(this));//default
    }
    
    //can't pass atomtype to constructor because atomtype needs this in its constructor
/*    public AtomFactoryMono(AtomType type) {
        atomType = type;
    }
*/ 

///**
// * Overrides parent class method and builds a single atom.
// */
//	protected Atom build(AtomTreeNodeGroup parent) {
//		return new Atom(space, atomType, 
//						AtomTreeNodeLeaf.FACTORY, 
//						sequencerFactory, 
//						parent);
//	}
    
   
    public boolean isGroupFactory() {return false;}
    
    public void setType(AtomType t) {atomType = t;}
    public AtomType type() {return atomType;}
    
    
    /**
     * Simply returns the given atom.
     */
    public Atom build(Atom atom) {return atom;}
    
    
}//end of AtomFactoryMono