package simulate;

public class PotentialIdealGas extends Potential implements PotentialHard, PotentialSoft {
    
    public PotentialIdealGas() {
    }
    
    public void bump(Atom atom1, Atom atom2) {return;}
    public double collisionTime(Atom atom1, Atom atom2) {return Double.MAX_VALUE;}
    
//    public double[] force(Atom atom1, Atom atom2) {return zero;}
//    public double virial(Atom atom1, Atom atom2) {return 0.0;}
    
}

