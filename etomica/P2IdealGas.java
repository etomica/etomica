package simulate;
import java.beans.Beans;
import java.awt.*;

public class P2IdealGas extends Potential2 {

  private double collisionDiameter = 0.1;

  public P2IdealGas() {
    super();
    nAtoms1 = 1;
    nAtoms2 = 1;
    potential = new Potential[nAtoms1][nAtoms2];
    potential[0][0] = new PotentialIdealGas();
  }
  
  public final Potential getPotential(Atom a1, Atom a2) {return potential[0][0];}
  
}


