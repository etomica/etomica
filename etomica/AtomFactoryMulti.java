package etomica;

/**
 * Builder of a multiatomic atom group of identical atoms.
 *
 * @author David Kofke
 */
public class AtomFactoryMulti extends AtomFactory {
    
    private AtomFactory childFactory;
    private int atomsPerGroup = 1;
    private final AtomType.Group groupType = new AtomType.Group(this);
    private Configuration configuration;
    
    public AtomFactoryMulti() {
        this(new AtomType.Sphere());
    }
    
    public AtomFactoryMulti(AtomType type) {
        childFactory = new AtomFactoryMono(type);
        //need better way to get space
        configuration = new ConfigurationLinear(Simulation.instance.space());
    }
    
    public Atom makeNewAtom(AtomGroup parent, int index) {
        AtomGroup group = new AtomGroup(parent, index, groupType, childFactory, 
                                        atomsPerGroup, true);
        configuration.initializeCoordinates(group);
    }
    
    public boolean producesAtomGroups() {return true;}
    
    public boolean vetoAddition(Atom a) {return a.type == groupType;} 
        //((GroupType)((AtomGroup)a).type).creator() == this;}
        
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
/*    public void setAtomsPerMolecule(int na) {
        if(na == atomsPerMolecule) return;  //do nothing is value isn't changings
        atomsPerMolecule = na;
//        if(parentSimulation == null) {return;}
        Iterator e = agents.values().iterator();
        while(e.hasNext()) {
            Agent a = (Agent)e.next();
            a.setNMolecules(a.nMolecules, true);//2nd argument indicates to make new molecules even though number of them is not changing
        }
    }
            
    /**
     * Accessor method for number of atoms per molecule
     * 
     * @return Present number of atoms in a molecule of this species.
     */
     public int getAtomsPerGroup() {return atomsPerGroup;}

}//end of AtomFactoryMulti
    
