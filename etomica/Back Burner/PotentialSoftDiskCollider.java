package simulate;
import java.io.*;
import java.beans.*;
import java.awt.*;

public class PotentialSoftDiskCollider extends PotentialSoftDisk {
        
    public PotentialSoftDiskCollider() {
        super();
    }
    
  public PairInteraction computePairInteraction(SpeciesElement element1, SpeciesElement element2) {
      PairInteraction pair = new PairInteraction();
      pair.rij = space.rij((MoleculeAtomic)element1,(MoleculeAtomic)element2);
      pair.rSquared = pair.rij[0]*pair.rij[0] + pair.rij[1]*pair.rij[1];  //change if not 2-d
      pair.force = new double[2];  //change if not 2-d
      if(pair.rSquared < potentialCutoffSquared) {
        factor = sigma*sigma/pair.rSquared;
        factor = factor*factor;  // (sig/r)^4
        factor = factor*factor*factor; //(sig/r)^12
        pair.energy = epsilon*factor;
        pair.virial = -n*pair.energy;
        factor = factor*n*epsilon/pair.rSquared;
        pair.force[0] = pair.rij[0]*factor;       //force on 2 due to 1
        pair.force[1] = pair.rij[1]*factor;
      }
      else {
        pair.energy = pair.virial = 0.0;
        pair.force[0] = pair.force[1] = 0.0;
      }
      return pair;
  }    
}