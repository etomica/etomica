package etomica;

/**
 * Master potential that sits that the top of the hierarchy of
 * potentials in a simulation.  Agents from master potential form
 * the base of all potential agents in each phase.
 *
 * @author David Kofke
 */
public class PotentialMaster extends PotentialGroup {
    
//    private PotentialGroup parentGroup;
    
    public PotentialMaster(Simulation sim) {
        super(sim);
    }
    
    public PotentialAgent makeAgent(Phase p) {
        return new Agent(this, p);
    }
    
    //---------------------------------------//
    
    public class Agent extends PotentialGroup.Agent {

//        private PotentialAgent.Hard.Linker firstHard;
        
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
 /*       public void addPotential(PotentialAgent potential) {
            super.addPotential(potential);
            if(potential instanceof PotentialAgent.Hard) {
                firstHard = new PotentialAgent.Hard.Linker((PotentialAgent.Hard)potential, firstHard);
            }
        }//end of addPotential
        
        //hard interface
        public void findCollisions(IteratorDirective id, IntegratorHardAbstract.CollisionHandler c) {
            for(PotentialAgent.Hard.Linker link=firstHard; link!=null; link=link.next()) {
                link.agent.findCollisions(id, c);
            }
        }
        
        //this will never be called
        public void bump(IntegratorHardAbstract.Agent agent) {
            System.out.println("Unexpected call to PotentialGroupMaster.Agent bump method");
        }
   */     
        
    }//end of Agent
}//end of PotentialGroup
    