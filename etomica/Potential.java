package simulate; 

public class Potential {

  PhaseSpace parentPhase;
  public Potential() {;}
    
 /**
  * Sets parentPhase
  *
  * @param p the phase to be identified as this Potential's parentPhase
  */
  public void setParentPhase(PhaseSpace p) {
        parentPhase = p;
  }
      
  public double energy(PhaseSpace.AtomPair pair) {
    System.out.println("super energy");
    return 0.0;
  }
  public final double energy(Atom atom1, Atom atom2) {     //prefer call with AtomPair argument to this one, if possible
    return energy(parentPhase.makeAtomPair(atom1, atom2));
  }
  
  public boolean overlap(PhaseSpace.AtomPair pair) {return false;}
}


