package etomica.association;
import etomica.Atom;

public abstract class BiasVolume implements AssociationDefinition, java.io.Serializable {
    
    public final Space space;
    public BiasVolume(Space s) {
        space = s;
    }
    
    public abstract double biasVolume();  
    public abstract Atom biasInsert(Atom atomA, Atom atomB);
    public abstract boolean isAssociated(Atom atomA, Atom atomB);
 
    
    
}