package etomica;

/**
 * Builder of a multiatomatomic atom group of identical atoms.
 *
 * @author David Kofke
 */
public class AtomFactoryMulti extends AtomFactory {
    
    AtomType atomType;
    
    public AtomFactoryMulti() {
        this(new AtomType.Sphere());
    }
    
    public AtomFactoryMulti(AtomType type) {
        atomType = type;
    }
    
    public Atom makeAtom(AtomGroup parent, int index) {
        return new Atom(parent, index, atomType);
    }
    
    public boolean producesAtomGroups() {return false;}
    
    public boolean vetoAddition(Atom a) {return (a instanceof AtomGroup);}
    
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
    public void setAtomsPerMolecule(int na) {
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
    public int getAtomsPerMolecule() {return atomsPerMolecule;}

}//end of AtomFactoryMono
    
