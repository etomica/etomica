package etomica;

/**
 * Builds an atom group that comprises a set of identically formed atoms or atomgroups.
 *
 * @author David Kofke
 */
public class AtomFactoryHomo extends AtomFactory {
    
    private AtomFactory childFactory;
    private int atomsPerGroup;
    
    /**
     * @param factory the factory that makes each of the identical children.
     */
    public AtomFactoryHomo(Simulation sim, AtomFactory factory) {
        this(sim, factory, 1);
    }
    /**
     * @param factory the factory that makes each of the identical children.
     * @param atoms the number of identical children per group (default is 1).
     */
    public AtomFactoryHomo(Simulation sim, AtomFactory factory, int atoms) {
        this(sim, factory, atoms, BondInitializer.NULL, new ConfigurationLinear(sim));
    }
    /**
     * @param factory the factory that makes each of the identical children.
     * @param atoms the number of identical children per group (default is 1).
     * @param config the configuration applied to each group that is built (default is Linear).
     */
    public AtomFactoryHomo(Simulation sim, AtomFactory factory, int atoms, BondInitializer bondInit,
                            Configuration config) {    
        super(sim);
        childFactory = factory;
        atomsPerGroup = atoms;
        bondInitializer = bondInit;
        configuration = config;
    }
    
    
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
    
    public boolean producesAtomGroups() {return true;}
    
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
    
