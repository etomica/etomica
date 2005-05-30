package etomica;

/**
 * Lorenz-Berthelot combining rule, in which the pair diameter is given
 * by the arithmetic mean of the individual atom diamters, and the pair
 * energy parameter is given by the geometric mean of the atom values.
 */
public class CombiningRuleLorenzBerthelot extends CombiningRule {
    
    public double sigma(Atom[] pair) {
        double sigma1 = ((Agent)pair[0].allatomAgents[0]).sigma;
        double sigma2 = ((Agent)pair[0].allatomAgents[0]).sigma;
        return (sigma1 == sigma2) ? sigma1 : 0.5*(sigma1+sigma2);
    }
    
    public double epsilon(Atom[] pair) {
        double epsilon1 = ((Agent)pair[0].allatomAgents[0]).epsilon;
        double epsilon2 = ((Agent)pair[1].allatomAgents[0]).epsilon;
        return (epsilon1 == epsilon2) ? epsilon1 : Math.sqrt(epsilon1*epsilon2);
    }
    
    public Object makeAgent(Atom a) {
        return new Agent();
    }
    
    public static class Agent {
        public double sigma = Default.ATOM_SIZE;
        public double epsilon = Default.POTENTIAL_WELL;
    }
        
}