package etomica;
import java.awt.*;

public class P2HardDisk extends Potential2 {

  private double collisionDiameter;
  private PotentialHardDisk onlyPotential;

  public P2HardDisk() {
    this(Simulation.instance);
  }
  public P2HardDisk(Simulation sim) {
    super(sim);
    collisionDiameter = Default.ATOM_SIZE;
    onlyPotential = new PotentialHardDisk(collisionDiameter);
    setCollisionDiameter(collisionDiameter);  //sets up neighbor distances 
  }
  
  public final Potential getPotential(Atom a1, Atom a2) {return onlyPotential;}

  public final double getCollisionDiameter() {return collisionDiameter;}
  public final void setCollisionDiameter(double d) {
    collisionDiameter = d;
    onlyPotential.setCollisionDiameter(d);
    setPotentialCutoff(d);
  }
}


