package simulate;
import java.awt.Component;

public abstract class Potential1 extends Component {

  Space space;
  int speciesIndex, nAtoms;
  Potential[][] potential;

  public Potential1() {
    speciesIndex = 0;
  }
    
  public abstract Potential getPotential(Atom a1, Atom a2);
  
  public void setPhase(Phase p) {
    for(int i=0; i<nAtoms-1; i++) {
        for(int j=0; j<nAtoms-1; j++) {
            potential[i][j].setParentPhase(p);
        }
    }
  }
  
  public void setSpace(Space s) {
    this.space = s;
/*    for(int i=0; i<nAtoms-1; i++) {
        for(int j=0; j<nAtoms-1; j++) {
            potential[i][j].space = s;
        }
    }
*/  }
        
  public final int getSpeciesIndex() {return this.speciesIndex;}  
  public final void setSpeciesIndex(int index) {this.speciesIndex = index;}
}


