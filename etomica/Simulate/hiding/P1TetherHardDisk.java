package simulate;

public class P1TetherHardDisk extends Potential1 {

  public PotentialTether potentialTether;
  public PotentialHardDisk potentialHD;
  
  public P1TetherHardDisk() {
    super();
    potentialTether = new PotentialTether();
    potentialHD = new PotentialHardDisk(0.1);
  }
    
  public Potential getPotential(Atom a1, Atom a2) {
    if(a1.getNextAtom() == a2 || a2.getNextAtom() == a1) { 
        return potentialTether;
    }
    else {
        return potentialHD;
    }
  }
  
  public final void setSpace(Space s) {
    this.space = s;
//    potentialTether.space = s;
//    potentialHD.space = s;
  }

  public final double getTetherLength() {return potentialTether.getTetherLength();}
  public final void setTetherLength(double t) {potentialTether.setTetherLength(t);}

  public final double getCollisionDiameter() {return potentialHD.getCollisionDiameter();}
  public final void setCollisionDiameter(double d) {potentialHD.setCollisionDiameter(d);}

}


