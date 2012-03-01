package etomica.association;
import etomica.api.IMolecule;

/**
 * Interface for a method that defines whether two molecules are considered associated.
 */
 
public interface AssociationDefinitionMolecule {
    
    public boolean isAssociated(IMolecule moleculeA, IMolecule moleculeB);
    
}
