package simulate;

public class P1TetherDumbbell extends Potential1 {

  public P1TetherDumbbell() {
    super();
    nAtoms = 2;
    potential = new Potential[nAtoms-1][nAtoms-1];
    potential[0][0] = new PotentialTether();
  }
    
  public Potential getPotential(Atom a1, Atom a2) {return potential[0][0];}
  
  public final double getTetherLength() {return ((PotentialTether)potential[0][0]).getTetherLength();}
  public final void setTetherLength(double t) {((PotentialTether)potential[0][0]).setTetherLength(t);}

}


