package etomica;

/**
 * Builds an atom group that comprises a set of identically formed atoms or atomgroups.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 09/23/02 (DAK) changed childFactory from private to protected so that AtomFactoryTree can implement a method to change it
  */
  
public class AtomFactoryHomo extends AtomFactory {
    
    protected AtomFactory childFactory;
    private int atomsPerGroup;
    
	//each constructor has a version that takes a Simulation, and one that takes a Space.
	//Simulation version is preferred.  It causes a handle to the simulation to be put in 
	//the AtomType, which is the way that the atom accesses the simulation
    public AtomFactoryHomo(Simulation sim, AtomFactory factory) {
        this(sim, sim.potentialMaster.sequencerFactory(), factory);
    }
	public AtomFactoryHomo(Simulation sim, AtomSequencer.Factory sequencerFactory, AtomFactory factory) {
		this(sim, sequencerFactory, factory, 1);
	}
	public AtomFactoryHomo(Simulation sim, AtomSequencer.Factory sequencerFactory, AtomFactory factory, int atoms) {
		this(sim, sequencerFactory, factory, atoms, BondInitializer.NULL, new ConfigurationLinear(sim.space));
	}
	public AtomFactoryHomo(Simulation sim, AtomSequencer.Factory sequencerFactory, AtomFactory factory, 
							int atoms, BondInitializer bondInit, Configuration config) {
		this(sim, sequencerFactory, AtomTreeNodeGroup.FACTORY, factory, atoms, bondInit, config);
	}   
	public AtomFactoryHomo(Simulation sim, AtomSequencer.Factory sequencerFactory, AtomTreeNode.Factory nodeFactory,
							AtomFactory factory, int atoms, BondInitializer bondInit, Configuration config) {	
		super(sim, sequencerFactory, nodeFactory);
		init(factory, atoms, bondInit, config);
	}
    /**
     * @param factory the factory that makes each of the identical children.
     */
    public AtomFactoryHomo(Space space, AtomSequencer.Factory sequencerFactory, AtomFactory factory) {
        this(space, sequencerFactory, factory, 1);
    }
    /**
     * @param factory the factory that makes each of the identical children.
     * @param atoms the number of identical children per group (default is 1).
     */
    public AtomFactoryHomo(Space space, AtomSequencer.Factory sequencerFactory, AtomFactory factory, int atoms) {
        this(space, sequencerFactory, factory, atoms, BondInitializer.NULL, new ConfigurationLinear(space));
    }
    /**
     * @param factory the factory that makes each of the identical children.
     * @param atoms the number of identical children per group (default is 1).
     * @param config the configuration applied to each group that is built (default is Linear).
     */
    public AtomFactoryHomo(Space space, AtomSequencer.Factory sequencerFactory, AtomFactory factory, 
                            int atoms, BondInitializer bondInit, Configuration config) {  
        this(space, sequencerFactory, AtomTreeNodeGroup.FACTORY, factory, atoms, bondInit, config);
    }
    public AtomFactoryHomo(Space space, AtomSequencer.Factory sequencerFactory, AtomTreeNode.Factory nodeFactory, 
    						AtomFactory factory, int atoms, BondInitializer bondInit, Configuration config) {
        super(space, sequencerFactory, nodeFactory);
        init(factory, atoms, bondInit, config);
    }
    
    private void init(AtomFactory factory, int atoms, 
    					BondInitializer bondInit, Configuration config) {                      
        childFactory = factory;
        atomsPerGroup = atoms;
        bondInitializer = bondInit;
        configuration = config;
        //set up fields of Group type (can't build sample atoms because factories defined by subclassing this one may not be ready to build at atom at this point)

        groupType.childSequencerClass = sequencerFactory.sequencerClass();
        groupType.childrenAreGroups = factory.isGroupFactory();
    }
    
    public boolean isGroupFactory() {return true;}
    
    /**
     * Constructs a new group using the given atom.
     */
     public Atom build(Atom group) {
        if(!group.creator().equals(this)) 
            throw new IllegalArgumentException("Error:  Cannot build atom from one created by a different factory");
        for(int i=0; i<atomsPerGroup; i++) {
            childFactory.build((AtomTreeNodeGroup)group.node);
        }
        bondInitializer.makeBonds(group);
        configuration.initializeCoordinates(group);
        return group;
     }
     
    /**
     * Returns the subfactory that produces each of the identical atoms
     * in the group made by this factory.
     */
    public AtomFactory childFactory() {return childFactory;}
        
    public boolean vetoAddition(Atom a) {return (a.creator() != childFactory);} 
        
 /*   public void renew(Atom a) {//need an exception in the case a is unrenewable
        if(a.type != groupType) return;  //throw exception
        AtomGroup group = (AtomGroup)a;
        int childCount = group.childCount();
        if(atomsPerGroup > childCount) {
            for(int i=childCount; i<atomsPerGroup; i++) group.addAtom(childFactory.makeAtom());
        }
        else if(atomsPerGroup < childCount) {
            for(int i=childCount; i>atomsPerGroup; i--) group.removeAtom(group.lastChild());
        }
        group.childIterator.reset();
        while(group.childIterator.hasNext()) {
            childFactory.renew(group.childIterator.next());
        }
        configuration.initializeCoordinates(group);
    }       */
        
    /**
     * Specifies the number of child atoms in each atom constructed by this factory.
     * 
     * @param na The new number of atoms per group
     */
    public void setAtomsPerGroup(int na) {
        atomsPerGroup = na;
    }
            
    /**
     * Accessor method for number of child atoms per group constructed.
     */
     public int getAtomsPerGroup() {return atomsPerGroup;}

}//end of AtomFactoryHomo
    
