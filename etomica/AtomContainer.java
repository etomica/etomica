package etomica;

/**
 * Interface for classes that can hold atoms.
 */
public interface AtomContainer {
    
 /** 
  * addAtom should first remove atom from original container, then
  * do what is necessary to add atom to new container
  */
    public void addAtom(Atom a);
    
 /**
  * removeAtom should be called only by the addAtom method of another container
  */
    public void removeAtom(Atom a);
    
} //end of interface AtomContainer
    