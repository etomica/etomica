package etomica;

/**
 * Master potential that sits that the top of the hierarchy of
 * potentials in a simulation.  Agents from master potential form
 * the base of all potential agents in each phase.
 *
 * @author David Kofke
 */
public class PotentialMaster extends PotentialGroup {
    
    public PotentialMaster(Simulation sim) {
        super(sim);
    }
    
    public PotentialAgent makeAgent(Phase p) {
        return new Agent(this, p);
    }
    
    //---------------------------------------//
    
    public class Agent extends PotentialGroup.Agent {
        
        public Agent(Potential potential, Phase phase) {
            super(potential, phase);
        }
        
        /**
         * Convenience extension of calculate method in PotentialAgent.  This method 
         * is applied if the PotentialCalculation argument implements PotentialCalculation.Sum.
         * The method returns this argument so that the sum can be accessed in-line with the call
         * to the calculate method.
         */
        public final PotentialCalculation.Sum calculate(IteratorDirective id, PotentialCalculation.Sum pa) {
            super.calculate(id, pa);
            return pa;
        }
        
    }//end of Agent
}//end of PotentialGroup
    