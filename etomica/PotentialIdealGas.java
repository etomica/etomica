package simulate;

public class PotentialIdealGas extends Potential implements PotentialHard, PotentialSoft {

    private final double[] zero = new double[Space.D];
    
    public PotentialIdealGas() {
        Space.uEa1(zero, 0.0);
    }
    
    public void bump(AtomC atom1, AtomC atom2) {return;}
    public double collisionTime(AtomC atom1, AtomC atom2) {return Double.MAX_VALUE;}
    
    public double[] force(Atom atom1, Atom atom2) {return zero;}
    public double virial(Atom atom1, Atom atom2) {return 0.0;}
    
}

