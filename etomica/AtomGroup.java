package etomica;

/**
 * Collection of atoms (or other atom groups). 
 *
 * @author David Kofke
 */
public class AtomGroup extends Atom implements java.io.Serializable {
    
    private Atom firstAtom;
    private Atom lastAtom;
    private final AtomFactory factory;
    private Configuration configuration;
    
    public AtomGroup(AtomGroup parent, int index, AtomFactory factory, 
                        int nAtoms, Configuration configuration) {
        super(parent, new AtomType.Group(), index);
        this.factory = factory;
        this.configuration = configuration;
        firstAtom = factory.makeAtom(this,0);
        lastAtom = firstAtom;
        for(int i=1; i<atomCount; i++) {
            lastAtom.setNextAtom(factory.makeAtom(this,i));
            lastAtom = lastAtom.nextAtom();
        }
        configuration.initializeCoordinates(this);
    }
    
    public void setConfiguration(Configuration c) {}
    
    public int atomCount() {}
    public int childCount() {}
}