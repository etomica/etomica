package simulate;

public class P1Null extends Potential1 {

  Potential p;
  
  public P1Null() {
    this(Simulation.instance);
  }
  public P1Null(Simulation sim) {
    super(sim);
    p = new PotentialIdealGas(sim);
  }
    
  public Potential getPotential(Atom a1, Atom a2) {return p;}
}


