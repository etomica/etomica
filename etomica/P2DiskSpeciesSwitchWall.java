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
  
  public void setSimulation(Simulation s) {
    super.setSimulation(s);
    setChangeSpeciesIndex(changeSpeciesIndex);
  }
  
  public final int getChangeSpeciesIndex() {
    return changeSpeciesIndex;
  }
  
  public final void setChangeSpeciesIndex(int i) {
    changeSpeciesIndex = i;
    Simulation simulation = potential[0][0].parentSimulation;
    if (simulation == null) {return;}
    for (Species s=simulation.firstSpecies(); s!=null; s=s.nextSpecies()) {
        if (s.getSpeciesIndex() == i) {
            ((PotentialHardDiskSpeciesSwitchWall)potential[0][0]).setChangeSpecies(s);
            return;
        }
    }
  }  
}


