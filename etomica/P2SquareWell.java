package simulate;
import simulate.units.*;

public class P2SquareWell extends Potential2 {

  private double coreDiameter;
  private double epsilon;
  private PotentialSquareWell onlyPotential;
  private double lambda = 1.5;

  public P2SquareWell() {
    this(Simulation.instance);
  }
  public P2SquareWell(Simulation sim) {
    super(sim);
    setCoreDiameter(Default.ATOM_SIZE);
    setEpsilon(Default.POTENTIAL_WELL);
    setLambda(lambda);  //set potentialCutoff, etc.
    onlyPotential = new PotentialSquareWell(coreDiameter,lambda,epsilon);
  }
    
  public final Potential getPotential(Atom a1, Atom a2) {return onlyPotential;}
  
  public final void setCoreDiameter(double d) {
    coreDiameter = d;
    setPotentialCutoff(coreDiameter*lambda);
    if(onlyPotential != null) onlyPotential.setCoreDiameter(coreDiameter);
  }
  public final double getCoreDiameter() {return coreDiameter;}
  public final Dimension getCoreDiameterDimension() {return Dimension.LENGTH;}
  
  //Well-width multiplier lambda
  public final double getLambda() {return lambda;}
  public final void setLambda(double lam) {
    lambda = lam;
    setPotentialCutoff(coreDiameter*lambda);
    if(onlyPotential != null) onlyPotential.setLambda(lambda);
  }
  
  public final void setEpsilon(double eps) {
    epsilon = eps;
    if(onlyPotential != null) onlyPotential.setEpsilon(epsilon);
  }
  public final double getEpsilon() {return epsilon;}
  public Dimension getEpsilonDimension() {return Dimension.ENERGY;}
}


