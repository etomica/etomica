package simulate;

public interface PotentialHard {
    
    public void bump(AtomC atom1, AtomC atom2);
    public double collisionTime(AtomC atom1, AtomC atom2);
    
}