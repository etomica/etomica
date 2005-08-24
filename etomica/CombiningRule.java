package etomica;

import etomica.atom.Atom;

/**
 * Encapulates the rule used to determine the size and/or energy parameter
 * for the interaction of a pair of atoms from the size and/or energy
 * values associated with each.
 */
public abstract class CombiningRule implements Atom.AgentSource, java.io.Serializable {
    
    public abstract double sigma(Atom[] pair);
    
    public abstract double epsilon(Atom[] pair);
    
    public abstract Object makeAgent(Atom a);
    
}
