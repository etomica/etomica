package simulate;
import java.beans.Beans;
import java.awt.*;

public class P2DiskWellWall extends Potential2 {

  private double collisionDiameter = 0.1;
  double lambda = 1.5;
  double epsilon = 300;
  double elasticity = 1.0;

  public P2DiskWellWall() {
    super();
    nAtoms1 = 1;
    nAtoms2 = 1;
    setSpecies2Index(1);
    potential = new Potential[nAtoms1][nAtoms2];
    potential[0][0] = new PotentialWellWall(collisionDiameter,lambda,epsilon,elasticity);
  }
  
  public final boolean isNeighbor(Molecule m1, Molecule m2) {
    return true;
  }
  
  public final Potential getPotential(Atom a1, Atom a2) {return potential[0][0];}
  
  public final double getCollisionDiameter() {return collisionDiameter;}
  public final void setCollisionDiameter(double d) {
    collisionDiameter = d;
    ((PotentialWellWall)potential[0][0]).setCoreDiameter(d);
    setPotentialCutoff(d);
  }

  public final double getEpsilon() {return ((PotentialWellWall)potential[0][0]).getEpsilon();}
  public final void setEpsilon(double d) {
    ((PotentialWellWall)potential[0][0]).setEpsilon(d);
  }
  
  public final double getLambda() {return ((PotentialWellWall)potential[0][0]).getLambda();}
  public final void setLambda(double d) {
    ((PotentialWellWall)potential[0][0]).setLambda(d);
  }
  
  public double getElasticity() {return elasticity;}
  public void setElasticity(double e) {
      elasticity = e;
  }

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
            moleculePixelDiameter = (int)(scale*collisionDiameter);
            g.setColor(Color.red);
            moleculePixelPositionX = getLocation().x;
            moleculePixelPositionY = getLocation().y;
            g.fillOval(0,0,moleculePixelDiameter,moleculePixelDiameter);
        }
    }
  }
    
  
}


