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
    
    public Atom makeNewAtom(AtomGroup parent, int index) {
        return new Atom(parent, index, atomType);
    }
    
    public boolean producesAtomGroups() {return false;}
    
    public boolean vetoAddition(Atom a) {return (a instanceof AtomGroup);}
    
}//end of AtomFactoryMono