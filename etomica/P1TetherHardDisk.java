package simulate;
import simulate.units.Dimension;

public class P1TetherHardDisk extends Potential1 {

  public PotentialTether potentialTether;
  public PotentialHardDisk potentialHD;
  
  public P1TetherHardDisk() {
    this(Simulation.instance);
  }
  public P1TetherHardDisk(Simulation sim) {
    super(sim);
    potentialTether = new PotentialTether();
    potentialHD = new PotentialHardDisk(Default.ATOM_SIZE);
  }
    
  public Potential getPotential(Atom a1, Atom a2) {
    if(a1.nextAtom() == a2 || a2.nextAtom() == a1) { 
        return potentialTether;
    }
    else {
        return potentialHD;
    }
  }
  
  public final double getTetherLength() {return potentialTether.getTetherLength();}
  public final void setTetherLength(double t) {potentialTether.setTetherLength(t);}
  public final Dimension getTetherLengthDimension() {return Dimension.LENGTH;}

  public final double getCollisionDiameter() {return potentialHD.getCollisionDiameter();}
  public final void setCollisionDiameter(double d) {potentialHD.setCollisionDiameter(d);}
  public final Dimension getCollisionDiameterDimension() {return Dimension.LENGTH;}

}


