package simulate;

public class P2DiskWall extends Potential2 {

  private double collisionDiameter;
  private PotentialHardDiskWall onlyPotential;

  public P2DiskWall() {
    this(Simulation.instance);
  }
  public P2DiskWall(Simulation sim) {
    super(sim);
    collisionDiameter = Default.ATOM_SIZE;
    onlyPotential = new PotentialHardDiskWall(collisionDiameter);
    setCollisionDiameter(collisionDiameter);   
  }
  
  public final Potential getPotential(Atom a1, Atom a2) {return onlyPotential;}
  
  public final double getCollisionDiameter() {return collisionDiameter;}
  public final void setCollisionDiameter(double d) {
    collisionDiameter = d;
    onlyPotential.setCollisionDiameter(d);
    setPotentialCutoff(d);
  }
    public void setIsothermal(boolean b) {onlyPotential.setIsothermal(b);}
    public boolean isIsothermal() {return onlyPotential.isIsothermal();}
    public void setTemperature(double t) {onlyPotential.setTemperature(t);}
    public double getTemperature() {return onlyPotential.getTemperature();}
}


