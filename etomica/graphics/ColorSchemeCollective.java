package etomica.graphics;

public interface ColorSchemeCollective {

    /**
     * Determine color of all atoms.  Color is not re-determined when
     * getAtomColor is called.
     */
    public void colorAllAtoms();

}