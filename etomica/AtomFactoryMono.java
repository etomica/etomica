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
    
 //   AtomType atomType;
//    private AtomSequencer.Factory seqFactory;
    
	//each constructor has a version that takes a Simulation, and one that takes a Space.
	//Simulation version is preferred.  It causes a handle to the simulation to be put in 
	//the AtomType, which is the way that the atom accesses the simulation

    /**
     * Constructor with neighborSequencerFactory and AtomType.Sphere defaults.
     */
    public AtomFactoryMono(Simulation sim) {
        this(sim, sim.iteratorFactory.neighborSequencerFactory());
    }
	public AtomFactoryMono(Simulation sim, AtomSequencer.Factory seqFactory) {
		super(sim, seqFactory, AtomTreeNodeLeaf.FACTORY);
		init();
	}
    public AtomFactoryMono(Space space, AtomSequencer.Factory seqFactory) {
        super(space, seqFactory, AtomTreeNodeLeaf.FACTORY);
        init();
    }
    
    private void init() {
        setType(new AtomType.Sphere(this));//default
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
    
    public boolean vetoAddition(Atom a) {return (a.type != atomType);}
    
//    public void renew(Atom a) {}
    
}//end of AtomFactoryMono