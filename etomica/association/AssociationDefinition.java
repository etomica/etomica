package etomica.association;
import etomica.Atom;

/**
 * Interface for a method that defines whether two atoms are considered associated.
 */
 
public interface AssociationDefinition {
    
    public boolean isAssociated(Atom atomA, Atom atomB);
    
}
