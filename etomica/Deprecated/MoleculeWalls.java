package simulate;

public class MoleculeWalls extends Molecule {
  
  public MoleculeWalls(Species s, int n) {
    super(s, n);
  }
  
  //Wall is an unconventional atom
  protected void makeAtoms() {
    atom = new Wall[nAtoms];
    for(int i=0; i<nAtoms; i++) {atom[i] = new Wall(this,i);}
  }
  
}
