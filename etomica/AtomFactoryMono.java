package etomica;

/**
 * Builder of a monoatomic atom group, which comprises just an Atom.
 *
 * @author David Kofke
 */
public class AtomFactoryMono extends AtomFactory {
    
    AtomType atomType;
    
    public AtomFactoryMono(Simulation sim) {
        super(sim);
        setType(new AtomType.Sphere(this));//default
    }
    
    //can't pass atomtype to constructor because atomtype needs this in its constructor
/*    public AtomFactoryMono(AtomType type) {
        atomType = type;
    }
*/    
    public void setType(AtomType t) {atomType = t;}
    public AtomType type() {return atomType;}
    
    /**
     * Overrides parent class method and builds a single atom.
     */
    protected Atom build(AtomTreeNodeGroup parent) {
        return new Atom(parentSimulation().space, atomType, 
                        AtomTreeNodeLeaf.FACTORY, 
                        parentSimulation.getIteratorFactory().neighborSequencerFactory(), 
                        parent);
    }
    
    /**
     * Simply returns the given atom.
     */
    public Atom build(Atom atom) {return atom;}
    
    public boolean vetoAddition(Atom a) {return (a.type != atomType);}
    
//    public void renew(Atom a) {}
    
}//end of AtomFactoryMono