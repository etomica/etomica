package simulate;

public interface PotentialHard {
    
    public void bump(Atom atom1, Atom atom2);
    public double collisionTime(Atom atom1, Atom atom2);
    
}