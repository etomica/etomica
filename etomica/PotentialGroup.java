package etomica;

public class PotentialGroup extends Potential {
    
//    private PotentialGroup parentGroup;
    
    public PotentialGroup() {
        this(Simulation.instance);
    }
    public PotentialGroup(Simulation sim) {
        super(sim);
    }
    
    public PotentialAgent makeAgent(Phase p) {
        return new Agent(this, p);
    }
    
    //---------------------------------------//
    public class Agent extends PotentialAgent {

        private PotentialAgent first;
        private PotentialAgent last;
        
        public Agent(Potential potential, Phase phase) {
            super(potential, phase);
        }
        
        public void calculate(IteratorDirective id, PotentialCalculation pc) {
            for(PotentialAgent p=first; p!=null; p=p.nextAgent()) {
                if(id.excludes(p)) continue; //see if potential is ok with iterator directive
                p.calculate(id, pc);
            }
        }
        
        /**
         * No iterator is associated with this agent.  This method has no action.
         */
        public void makeDefaultIterator() {}
        
        public void addPotential(PotentialAgent potential) {
            if(first == null) first = potential;
            if(last != null) last.setNextAgent(potential);
            last = potential;
        }
    }//end of Agent
        
}//end of PotentialGroup
    