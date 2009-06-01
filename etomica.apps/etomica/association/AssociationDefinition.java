package etomica.association;
import etomica.api.IAtom;

/**
 * Interface for a method that defines whether two atoms are considered associated.
 */
 
public interface AssociationDefinition {
    
    public boolean isAssociated(IAtom atomA, IAtom atomB);
    
}
