package etomica;

public class PotentialGroupMaster extends PotentialGroup {
    
//    private PotentialGroup parentGroup;
    
    public PotentialGroupMaster(Simulation sim) {
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
    