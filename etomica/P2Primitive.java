package simulate;
import simulate.units.*;

/**
 * "Primitive model" potential.  Hard sphere plus Coulomb interaction
 */
public class P2Primitive extends Potential2 {

  private double sigma;
  private double cutoff = 2.5;
  PotentialPrimitive onlyPotential = new PotentialPrimitive(sigma,cutoff);

  public P2Primitive() {
    this(Simulation.instance);
  }
  public P2Primitive(Simulation sim) {
    super(sim);
    setSigma(Default.ATOM_SIZE);
  }
  
  public final Potential getPotential(Atom a1, Atom a2) {return onlyPotential;}
  
  public final double getCutoff() {return cutoff;}
  public final void setCutoff(double c) {
    cutoff = c;
    onlyPotential.setCutoff(c);}
  
  public final void setSigma(double d) {
    sigma = d;
    onlyPotential.setSigma(sigma);
  }
  public final double getSigma() {return sigma;}
  public final double sigma() {return sigma;}
  public final Dimension getSigmaDimension() {return Dimension.LENGTH;}
    
}


