package etomica;

/**
 * Builds an atom group that comprises a set of differently formed atoms or atomgroups.
 * Each child atoms is constructed by a different atom factory, which are set as an
 * array of atomfactories given in the constructor.
 *
 * @author David Kofke
 */
public class AtomFactoryHetero extends AtomFactory {
    
    private AtomFactory[] childFactory;
    private final AtomType.Group groupType = new AtomType.Group(this);
    
    /**
     * @param factory the factory that makes each of the identical children.
     */
    public AtomFactoryHetero(Simulation sim, AtomFactory[] factory) {
        this(sim, factory, new ConfigurationLinear(sim.space()));
    }
    /**
     * @param factory the factory that makes each of the identical children.
     * @param atoms the number of identical children per group (default is 1).
     * @param config the configuration applied to each group that is built (default is Linear).
     */
    public AtomFactoryHetero(Simulation sim, AtomFactory[] factory, 
                            Configuration config) {    
        super(sim);
        childFactory = factory;
        configuration = config;
    }
    
    /**
     * Constructs a new group.
     */
    public Atom build() {
        AtomGroup group = new AtomGroup(parentSimulation.space(), groupType);
        for(int i=0; i<childFactory.length; i++) {
            group.addAtom(childFactory[i].build());
        }
        configuration.initializeCoordinates(group);
        return group;
    }
    
    /**
     * Returns the array of subfactories that produces each of the identical atoms
     * in the group made by this factory.
     */
    public AtomFactory[] childFactory() {return childFactory;}
    
    public boolean producesAtomGroups() {return true;}
    
    public boolean vetoAddition(Atom a) {return true;} 
        
    public void renew(Atom a) {//need an exception in the case a is unrenewable
        if(a.type != groupType) return;  //throw exception
        configuration.initializeCoordinates((AtomGroup)a);
    }       
        
}//end of AtomFactoryHomo
    
