package etomica;

/**
 * Encapulates the rule used to determine the size and/or energy parameter
 * for the interaction of a pair of atoms from the size and/or energy
 * values associated with each.
 */
public abstract class CombiningRule implements Atom.AgentSource {
    
    public abstract double sigma(AtomPair pair);
    
    public abstract double epsilon(AtomPair pair);
    
    public abstract Object makeAgent(Atom a);
    
}