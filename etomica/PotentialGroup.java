package etomica;

public class PotentialGroup extends PotentialAbstract {
    
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
        
        public Agent(PotentialAbstract potential, Phase phase) {
            super(potential, phase);
        }
        
        public double energy(IteratorDirective id) {
            double sum = 0.0;
            for(PotentialAgent p=first; p!=null; p=p.nextAgent()) {
                sum += p.energy(id);
            }
            return sum;
        }
        
        public void makeIterator() {}
        
        public void addPotential(PotentialAgent potential) {
            if(first == null) first = potential;
            if(last != null) last.setNextAgent(potential);
            last = potential;
        }
    }
        
}//end of PotentialGroup
    