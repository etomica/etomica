package simulate;

public class MoleculeWall extends Molecule {
  
  public MoleculeWall(Species s) {
    super(s, 1);  //Molecule is one AtomWall atom
  }
  
  //AtomWall is an unconventional atom
  protected void makeAtoms() {
    atom = new AtomWall[nAtoms];
    for(int i=0; i<nAtoms; i++) {atom[i] = new AtomWall(this,i);}
  }
  
}
