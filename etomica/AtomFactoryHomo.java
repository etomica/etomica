package etomica;

/**
 * Builds an atom group that comprises a set of identically formed atoms or atomgroups.
 *
 * @author David Kofke
 */
public class AtomFactoryHomo extends AtomFactory {
    
    private AtomFactory childFactory;
    private int atomsPerGroup;
    private final AtomType.Group groupType = new AtomType.Group(this);
    private Configuration configuration;
    
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
        this(sim, factory, atoms, new ConfigurationLinear(sim.space()));
    }
    /**
     * @param factory the factory that makes each of the identical children.
     * @param atoms the number of identical children per group (default is 1).
     * @param config the configuration applied to each group that is built (default is Linear).
     */
    public AtomFactoryHomo(Simulation sim, AtomFactory factory, int atoms,
                            Configuration config) {    
        super(sim);
        childFactory = factory;
        atomsPerGroup = atoms;
        configuration = config;
    }
    
    /**
     * Constructs a new group.
     */
    public Atom build() {
        AtomGroup group = new AtomGroup(parentSimulation.space(), groupType);
        for(int i=0; i<atomsPerGroup; i++) {
            group.addAtom(childFactory.build());
        }
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
        
    public void renew(Atom a) {//need an exception in the case a is unrenewable
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
    }       
        
    /**
     * Specifies the number of atoms in a molecule of this species.
     * Since the number of atoms in a molecule cannot be changed once the molecule 
     * is constructed, to have the atom/molecule change take effect it is necessary to 
     * create new molecules of this species in each phase.  Thus the method invokes
     * the setNMolecules method of each agent of this species, replacing all the existing
     * molecules with ones having the newly prescribed value of atomsPerMolecule.
     * 
     * @param na The new number of atoms per molecule
     */
    public void setAtomsPerGroup(int na) {
        if(na == atomsPerGroup) return;  //do nothing if value isn't changings
        atomsPerGroup = na;
//        if(parentSimulation == null) {return;}
/*        Iterator e = agents.values().iterator();
        while(e.hasNext()) {
            Agent a = (Agent)e.next();
            a.setNMolecules(a.nMolecules, true);//2nd argument indicates to make new molecules even though number of them is not changing
        }*/
    }
            
    /**
     * Accessor method for number of atoms per molecule
     * 
     * @return Present number of atoms in a molecule of this species.
     */
     public int getAtomsPerGroup() {return atomsPerGroup;}

}//end of AtomFactoryHomo
    
