package simulate;
import java.beans.Beans;
import java.awt.*;

public class P2LennardJones extends Potential2 {

  private double sigma = 0.1;
  private double epsilon = 300;

  public P2LennardJones() {
    super();
    nAtoms1 = 1;
    nAtoms2 = 1;
    potential = new Potential[nAtoms1][nAtoms2];
    potential[0][0] = new PotentialLennardJones(sigma,epsilon);
  }
  
  public final boolean isNeighbor(Molecule m1, Molecule m2) {
    return (space.r1Mr2_S(m1.COM(), m2.COM()) < squareNeighborRadius);
  }
  
  public final Potential getPotential(Atom a1, Atom a2) {return potential[0][0];}
  
  public final double getSigma() {return sigma;}
  public final void setSigma(double s) {
    sigma = s;
    ((PotentialLennardJones)potential[0][0]).setSigma(s);
  }
    
  public final double getEpsilon() {return epsilon;}
  public final void setEpsilon(double eps) {
    epsilon = eps;
    ((PotentialSquareWell)potential[0][0]).setEpsilon(eps);
  }
  public final void setEpsilonInt(int eps) {
    setEpsilon((double)eps);
  }
  
  // Paint a red disk and blue circle at design time to show size of 
  // core diameter and well
  
  // Paint a red disk at design time to show size of collision diameter
  
  public void paint(Graphics g) {
    int moleculePixelPositionX, moleculePixelPositionY, moleculePixelDiameter;
    int[] simulationPixelDimensions = {-1, -1}; // pixel width (0) and height (1) of simulation box less boundaries. 
    if(Beans.isDesignTime()) {
        if(getParent() != null) {
            Component par = getParent();
            simulationPixelDimensions[0] = par.getSize().width;
            simulationPixelDimensions[1] = par.getSize().height;
            double scale = Math.max(simulationPixelDimensions[0], simulationPixelDimensions[1]);
            moleculePixelDiameter = (int)(scale*sigma);
            g.setColor(Color.red);
            moleculePixelPositionX = getLocation().x;
            moleculePixelPositionY = getLocation().y;
            g.fillOval(0,0,moleculePixelDiameter,moleculePixelDiameter);
        }
    }
  }
}


