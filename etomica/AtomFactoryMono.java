package etomica;

/**
 * Builder of a monoatomic atom group, which comprises just an Atom.
 *
 * @author David Kofke
 */
public class AtomFactoryMono extends AtomFactory {
    
    AtomType atomType;
    
    public AtomFactoryMono(Space space) {
        super(space);
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
        return new Atom(space, atomType, new AtomTreeNodeLeaf());
    }
    
    public boolean vetoAddition(Atom a) {return (a.type != atomType);}
    
    public void renew(Atom a) {}
    
}//end of AtomFactoryMono