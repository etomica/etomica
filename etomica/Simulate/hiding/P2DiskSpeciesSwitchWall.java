package simulate;
import java.beans.Beans;
import java.awt.*;

public class P2DiskSpeciesSwitchWall extends Potential2 {

  int changeSpeciesIndex;
    
  public P2DiskSpeciesSwitchWall() {
    super();
    setSize(30,30);
    nAtoms1 = 1;
    nAtoms2 = 1;
    potential = new Potential[nAtoms1][nAtoms2];
    potential[0][0] = new PotentialHardDiskSpeciesSwitchWall();
    setChangeSpeciesIndex(0);
  }
  
  public final boolean isNeighbor(Molecule m1, Molecule m2) {
    return true;
  }
  
  public final Potential getPotential(Atom a1, Atom a2) {return potential[0][0];}
  
  public void setPhase(Phase p) {
    super.setPhase(p);
    setChangeSpeciesIndex(changeSpeciesIndex);
  }
  
  public final int getChangeSpeciesIndex() {
    return changeSpeciesIndex;
  }
  
  public final void setChangeSpeciesIndex(int i) {
    changeSpeciesIndex = i;
    Phase phase = potential[0][0].parentPhase;
    if (phase == null) {return;}
    for (Species s=phase.firstSpecies(); s!=null; s=s.getNextSpecies()) {
        if (s.getSpeciesIndex() == i) {
            ((PotentialHardDiskSpeciesSwitchWall)potential[0][0]).setChangeSpecies(s);
            return;
        }
    }
  }  
}


