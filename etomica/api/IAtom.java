package etomica.api;


public interface IAtom {

    /**
     * @return the Atom type, holding properties held in common with other 
     * atoms made by this atom's factory.
     */
    public IAtomType getType();
}