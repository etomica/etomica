package etomica;

public class P1Null extends Potential1 implements EtomicaElement {

  Potential p;
  
  public P1Null() {
    this(Simulation.instance);
  }
  public P1Null(Simulation sim) {
    super(sim);
    p = new PotentialIdealGas(sim);
  }
    
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("No intramolecular interaction");
        return info;
    }
    
  public Potential getPotential(Atom a1, Atom a2) {return p;}
}


