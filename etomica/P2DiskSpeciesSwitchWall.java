package simulate;
import java.beans.Beans;
import java.awt.*;

public class P2DiskSpeciesSwitchWall extends Potential2 {

  int changeSpeciesIndex;
  private final PotentialHardDiskSpeciesSwitchWall onlyPotential;
    
  public P2DiskSpeciesSwitchWall() {
    this(Simulation.instance);
  }
  public P2DiskSpeciesSwitchWall(Simulation sim) {
    super(sim);
    onlyPotential = new PotentialHardDiskSpeciesSwitchWall();
    setChangeSpeciesIndex(0);
  }
  
  public final boolean isNeighbor(Molecule m1, Molecule m2) {
    return true;
  }
  
  public final Potential getPotential(Atom a1, Atom a2) {return onlyPotential;}
  
  public final int getChangeSpeciesIndex() {
    return changeSpeciesIndex;
  }
  
  public final void setChangeSpeciesIndex(int i) {
    changeSpeciesIndex = i;
    for (Species s=parentSimulation().firstSpecies(); s!=null; s=s.nextSpecies()) {
        if (s.getSpeciesIndex() == i) {
            onlyPotential.setChangeSpecies(s);
            return;
        }
    }
  }  
}


