package simulate;
import simulate.units.*;

public class P2LennardJones extends Potential2 {

  private final PotentialLJ onlyPotential;

  public P2LennardJones() {
    this(Simulation.instance);
  }
  public P2LennardJones(Simulation sim) {
    super(sim);
    onlyPotential = new PotentialLJ(Default.ATOM_SIZE,Default.POTENTIAL_WELL,Default.POTENTIAL_CUTOFF);
  }
  
  public final Potential getPotential(Atom a1, Atom a2) {return onlyPotential;}
  
  public final double getCutoff() {return onlyPotential.getCutoff();}
  public final void setCutoff(double c) {onlyPotential.setCutoff(c);}
  public final Dimension getCutoffDimension() {return Dimension.NULL;}
  
  public final void setSigma(double sigma) {onlyPotential.setSigma(sigma);}
  public final double getSigma() {return onlyPotential.getSigma();}
  public final Dimension getSigmaDimension() {return Dimension.LENGTH;}
    
  public final void setEpsilon(double epsilon) {onlyPotential.setEpsilon(epsilon);}
  public final double getEpsilon() {return onlyPotential.getEpsilon();}
  public final Dimension getEpsilonDimension() {return Dimension.ENERGY;}
  
}


