package simulate;

public interface PotentialHard {
    
    public void bump(AtomHard atom1, AtomHard atom2);
    public double collisionTime(AtomHard atom1, AtomHard atom2);
    
}