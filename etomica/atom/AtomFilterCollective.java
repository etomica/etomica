package etomica.atom;

/**
 * Filter that needs to be informed whenever something has changed and the
 * filtered status of atoms should be redetermined.
 * 
 * @author Andrew Schultz
 */
public interface AtomFilterCollective extends AtomFilter {

    /**
     * Informs the filter that something in the system has changed (atom
     * positions, velocityies, etc) has changed and the filtered status of the
     * atoms should be determined.
     */
    public void resetFilter();
}
