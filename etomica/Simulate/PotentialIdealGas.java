package simulate;

public class PotentialIdealGas extends Potential implements PotentialHard, PotentialSoft {
    
    public PotentialIdealGas() {
    }
    
    public void bump(AtomPair pair) {return;}
    public double collisionTime(AtomPair pair) {return Double.MAX_VALUE;}
    
//    public double[] force(Atom atom1, Atom atom2) {return zero;}
//    public double virial(Atom atom1, Atom atom2) {return 0.0;}
    
}

