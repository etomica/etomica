package simulate;

public class P1Null extends Potential1 {

  Potential p;
  
  public P1Null() {
    super();
    p = new PotentialIdealGas();
  }
    
  public Potential getPotential(Atom a1, Atom a2) {return p;}
}


