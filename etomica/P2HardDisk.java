package simulate;
import java.beans.Beans;
import java.awt.*;

public class P2HardDisk extends Potential2 {

  private double collisionDiameter = 0.1;

  public P2HardDisk() {
    super();
    setSize(30,30);
    nAtoms1 = 1;
    nAtoms2 = 1;
    potential = new Potential[nAtoms1][nAtoms2];
    potential[0][0] = new PotentialHardDisk(collisionDiameter);
    setCollisionDiameter(collisionDiameter);  //sets up neighbor distances 
  }
  
  public final boolean isNeighbor(Molecule m1, Molecule m2) {
    return (space.r1Mr2_S(m1.COM(), m2.COM()) < squareNeighborRadius);
  }
  
  public final Potential getPotential(Atom a1, Atom a2) {return potential[0][0];}
  
  public final double getCollisionDiameter() {return collisionDiameter;}
  public final void setCollisionDiameter(double d) {
    collisionDiameter = d;
    ((PotentialHardDisk)potential[0][0]).setCollisionDiameter(d);
    setPotentialCutoff(d);
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


