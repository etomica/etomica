package simulate; 

public class Potential {

  Phase parentPhase;
  public Potential() {;}
    
 /**
  * Sets parentPhase
  *
  * @param p the phase to be identified as this Potential's parentPhase
  */
  public void setParentPhase(Phase p) {
        parentPhase = p;
  }
      
  public double energy(Atom atom1, Atom atom2) {
    System.out.println("super energy");
    return 0.0;
  }
  
  public double energy(AtomC atom1, AtomC atom2) {
    System.out.println("super energy C");
    return 0.0;
  }

  public boolean overlap(Atom atom1, Atom atom2) {return false;}
}


