package etomica;

/**
 * Builder of a monoatomic atom group, which comprises just an Atom.
 *
 * @author David Kofke
 */
public class AtomFactoryMono extends AtomFactory {
    
    AtomType atomType;
    private final Simulation simulation;
    
    public AtomFactoryMono(Space space) {
        super(space);
        simulation = Simulation.instance;//needs work
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
     * Builds a single atom.
     */
    protected Atom build() {
        return new Atom(space, atomType, 
                        AtomTreeNodeLeaf.FACTORY, 
                        simulation.getIteratorFactory().neighborSequencerFactory());
    }
    
    public boolean vetoAddition(Atom a) {return (a.type != atomType);}
    
//    public void renew(Atom a) {}
    
}//end of AtomFactoryMono