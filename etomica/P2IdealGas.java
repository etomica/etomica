package simulate;
import java.beans.Beans;
import java.awt.*;

public class P2IdealGas extends Potential2 {

  private final PotentialIdealGas onlyPotential;
  
  public P2IdealGas() {
    this(Simulation.instance);
  }
  public P2IdealGas(Simulation sim) {
    super(sim);
    onlyPotential = new PotentialIdealGas(sim);
  }
  
  public final Potential getPotential(Atom a1, Atom a2) {return onlyPotential;}
  
}


