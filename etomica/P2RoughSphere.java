package simulate;
import java.awt.*;

public class P2RoughSphere extends Potential2 {

  private double collisionDiameter;
  private PotentialHardDisk onlyPotential;

  public P2RoughSphere() {
    this(Simulation.instance);
  }
  public P2RoughSphere(Simulation sim) {
    super(sim);
    collisionDiameter = Default.ATOM_SIZE;
    onlyPotential = new PotentialRoughSphere(collisionDiameter);
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


