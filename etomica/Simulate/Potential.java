package simulate; 

public class Potential {

  Simulation parentSimulation;
  public Potential() {;}
    
 /**
  * Sets parentPhase
  *
  * @param p the phase to be identified as this Potential's parentPhase
  */
  public void setParentSimulation(Simulation s) {
        parentSimulation = s;
  }
      
  public double energy(AtomPair pair) {
    System.out.println("super energy");
    return 0.0;
  }
  
  public boolean overlap(AtomPair pair) {return false;}
}


