package etomica;

//should not be used to hold information about the state of iteration.
//could be accessed by more than one iterator at a time

public abstract class AtomSequencer extends AtomLinker {
    
    public AtomSequencer(Atom a) {super(a);}
    
    /**
     * Notifies sequencer that atom has been moved to a new position.
     */
    public abstract void moveNotify();
    
    public abstract void setParentNotify(AtomTreeNodeGroup newParent);
    
    public abstract boolean preceeds(Atom a);
    
    public interface Factory {
        public AtomSequencer makeSequencer(Atom atom);
    }
    
}