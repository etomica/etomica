package simulate;
import java.beans.Beans;
import java.awt.*;

public class P2SquareWell extends Potential2 {

  private double coreDiameter = 0.1;
  private double lambda = 1.5;
  private double epsilon = 300;

  public P2SquareWell() {
    super();
    nAtoms1 = 1;
    nAtoms2 = 1;
    potential = new Potential[nAtoms1][nAtoms2];
    potential[0][0] = new PotentialSquareWell(coreDiameter,lambda,epsilon);
    setLambda(lambda);  //set potentialCutoff, etc.
  }
    
  public final Potential getPotential(Atom a1, Atom a2) {return potential[0][0];}
  
  public final double getCoreDiameter() {return coreDiameter;}
  public final void setCoreDiameter(double d) {
    coreDiameter = d;
    setPotentialCutoff(coreDiameter*lambda);
    ((PotentialSquareWell)potential[0][0]).setCoreDiameter(d);
  }
  
  public final double getLambda() {return lambda;}
  public final void setLambda(double lam) {
    lambda = lam;
    setPotentialCutoff(coreDiameter*lambda);
    ((PotentialSquareWell)potential[0][0]).setLambda(lam);
  }
  
  public final double getEpsilon() {return epsilon;}
  public final void setEpsilon(double eps) {
    epsilon = eps;
    ((PotentialSquareWell)potential[0][0]).setEpsilon(eps);
  }
  public final void setEpsilon(int eps) {
    setEpsilon((double)eps);
  }
  
  // Paint a red disk and blue circle at design time to show size of 
  // core diameter and well
  
  public void paint(Graphics g) {
    int moleculePixelPositionX, moleculePixelPositionY, moleculePixelDiameter, wellPixelDiameter;
    int[] simulationPixelDimensions = {-1, -1}; // pixel width (0) and height (1) of simulation box less boundaries. 
    if(Beans.isDesignTime()) {
        if(getParent() != null) {
            Component par = getParent();
            simulationPixelDimensions[0] = par.getSize().width;
            simulationPixelDimensions[1] = par.getSize().height;
            double scale = Math.max(simulationPixelDimensions[0], simulationPixelDimensions[1]);
            moleculePixelDiameter = (int)(scale*coreDiameter);
            wellPixelDiameter = (int)(scale*lambda*coreDiameter);
            g.setColor(Color.red);
 //           moleculePixelPositionX = getLocation().x;
 //           moleculePixelPositionY = getLocation().y;
            g.fillOval((wellPixelDiameter-moleculePixelDiameter)/2,(wellPixelDiameter-moleculePixelDiameter)/2,moleculePixelDiameter,moleculePixelDiameter);
            g.setColor(Color.blue);
            g.drawOval(0,0,wellPixelDiameter,wellPixelDiameter);
        }
    }
  }  
}


