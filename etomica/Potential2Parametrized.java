package etomica;

public class Potential2Parametrized extends Potential2 {
    
    final Potential2 potential;
    final CombiningRule combiningRule;
    
    public Potential2SoftParametrized(PotentialGroup parent, Potential2 potential,
                                        CombiningRule rule) {
        this.potential = potential;
    }
    
    public double setParameters(AtomPair pair) {
        
}