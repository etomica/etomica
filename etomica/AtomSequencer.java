package etomica;

//should not be used to hold information about the state of iteration.
//could be accessed by more than one iterator at a time

public interface AtomSequencer {
    
    public Atom nextAtom();
    public Atom previousAtom();
    
    public void setNextAtom(Atom a);
    public void setPreviousAtom(Atom a);
    public void clearPreviousAtom();
    
    public boolean preceeds(Atom a);
    
}