package simulate;

public class MoleculeCell extends Molecule {
  
  public double;
  public double cellRadius = 0.1;
  public double bondLength = 0.01;
  
  public MoleculeCell(Species s, int nT) {
    super(s, nT+1);
    nTethers = nT;
  }
  
  protected void makeAtoms() {
    double[] e = new double[Space.D];
    double[] r0 = new double[Space.D];
    atom = new Atom[nAtoms];
    atom[0] = new Atom(this,0,1.0,cellRadius);
    for(int i=1; i<=nTethers; i++) {
        double theta = 2*Math.PI*(i-1)/(double)nTethers;
        e[0] = Math.sin(theta);
        e[1] = Math.cos(theta);
        Space.uEa1Tv1(r0,cellRadius,e);
        atom[i] = new AtomTether(this,i,r0,bondLength,e);
    }
  }
  
  public void rotate(double theta) {
    for(AtomTether a=(AtomTether)lastAtom; a!=firstAtom; a=(AtomTether)a.getNextAtom()) {
        Space.rotate(a.r0,firstAtom.r,theta);
        Space.rotate(a.r1,firstAtom.r,theta);
    }
  }
  
}
