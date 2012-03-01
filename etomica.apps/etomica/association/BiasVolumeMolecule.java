package etomica.association;

import etomica.api.IMolecule;
import etomica.space.ISpace;

public abstract class BiasVolumeMolecule implements AssociationDefinitionMolecule, java.io.Serializable {
    
    public final ISpace space;
    public BiasVolumeMolecule(ISpace s) {
        space = s;
    }
    
    public abstract double biasVolume();  
    public abstract void biasInsert(IMolecule moleculeA, IMolecule moleculeB);
    public abstract boolean isAssociated(IMolecule moleculeA, IMolecule moleculeB);
 
    
    
}