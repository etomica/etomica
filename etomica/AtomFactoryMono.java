package etomica;

/**
 * Builder of a monoatomic atom group, which is just an Atom.
 *
 * @author David Kofke
 */
public class AtomFactoryMono extends AtomFactory {
    
    AtomType atomType;
    
    public AtomFactoryMono() {
        this(new AtomType.Sphere());
    }
    
    public AtomFactoryMono(AtomType type) {
        atomType = type;
    }
    
    public Atom makeAtom(AtomGroup parent, int index) {
        return new Atom(parent, index, atomType);
    }
    
}//end of AtomFactoryMono