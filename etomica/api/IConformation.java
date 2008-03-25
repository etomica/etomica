package etomica.api;


public interface IConformation {

    /**
     * Defined by subclass to assign coordinates to the atoms in the given list.
     */
    public void initializePositions(IAtomSet atomList);

}