package etomica.association;
import etomica.api.IAtom;
import etomica.space.ISpace;

public abstract class BiasVolume implements AssociationDefinition, java.io.Serializable {
    
    public final ISpace space;
    public BiasVolume(ISpace s) {
        space = s;
    }
    
    public abstract double biasVolume();  
    public abstract void biasInsert(IAtom atomA, IAtom atomB);
    public abstract boolean isAssociated(IAtom atomA, IAtom atomB);
 
    
    
}