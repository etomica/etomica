package etomica.association;
import etomica.Atom;

public interface AssociationDefinition {
    
    public boolean isAssociated(Atom atomA, Atom atomB);
    
}
