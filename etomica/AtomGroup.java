package etomica;

public class AtomGroup extends Atom implements java.io.Serializable {
    
    private AtomGroup first;
    private AtomGroup last;
    private AtomGroup next;
    private AtomGroup previous;
    
    public AtomGroup(AtomGroup parent, AtomType type) {
        super(parent, type, 0);
    }
    
 //   public Species[] representedSpecies  //species represented in this group of atoms
}