package simulate;
import java.beans.Beans;
import java.awt.*;

public class P2LennardJones extends Potential2 {

  private double sigma = 0.1;
  private double epsilon = 300;
  private double cutoff = 2.5;
  PotentialLJ onlyPotential = new PotentialLJ(sigma,epsilon,cutoff);

  public P2LennardJones() {
    super();
    nAtoms1 = 1;
    nAtoms2 = 1;
    potential = new Potential[1][1];
    potential[0][0] = onlyPotential;  //another copy for superclass methods
  }
  
//  public final boolean isNeighbor(Molecule m1, Molecule m2) {
//    return (space.r1Mr2_S(m1.COM(), m2.COM()) < squareNeighborRadius);
//  }
  
  public final Potential getPotential(Atom a1, Atom a2) {return onlyPotential;}
  
  public final double getCutoff() {return cutoff;}
  public final void setCutoff(double c) {
    cutoff = c;
    onlyPotential.setCutoff(c);}
  
  public final double getSigma() {return sigma;}
  public final void setSigma(double s) {
    sigma = s;
    onlyPotential.setSigma(s);
  }
    
  public final double getEpsilon() {return epsilon;}
  public final void setEpsilon(double eps) {
    epsilon = eps;
    onlyPotential.setEpsilon(eps);
  }
  public final void setEpsilon(int eps) {
    setEpsilon((double)eps);
  }
  
  // Paint a red disk at design time to show size of sigma
  
  public void paint(Graphics g) {
    int moleculePixelPositionX, moleculePixelPositionY, moleculePixelDiameter, wellPixelDiameter;
    int[] simulationPixelDimensions = {-1, -1}; // pixel width (0) and height (1) of simulation box less boundaries. 
    if(Beans.isDesignTime()) {
        if(getParent() != null) {
            Component par = getParent();
            simulationPixelDimensions[0] = par.getSize().width;
            simulationPixelDimensions[1] = par.getSize().height;
            double scale = Math.max(simulationPixelDimensions[0], simulationPixelDimensions[1]);
            moleculePixelDiameter = (int)(scale*sigma);
            g.setColor(Color.red);
            g.fillOval(0,0,moleculePixelDiameter,moleculePixelDiameter);
        }
    }
  }  
}


