package etomica;

public class AtomGroup implements java.io.Serializable {
    
    private AtomGroup first;
    private AtomGroup last;
    private AtomGroup next;
    private AtomGroup previous;
    private AtomGroup parentGroup;
    
    public AtomGroup(AtomGroup parent) {
        parentGroup = parent;
    }
    
 //   public Species[] representedSpecies  //species represented in this group of atoms
}