package etomica;

public abstract class PotentialAgent implements java.io.Serializable {

    protected Phase parentPhase;  //want final, but won't compile
    protected PotentialAbstract parentPotential;
    private PotentialAgent next;
    
    /**
        * @param p The phase in which this agent will be placed
        */
    public PotentialAgent(PotentialAbstract potential, Phase phase) {
        parentPotential = potential;
        parentPhase = phase;
        parentPotential.agents().put(phase, this);
        makeIterator();
    }
    public final Phase parentPhase() {return parentPhase;}
        
    public PotentialAgent next() {return next;}
//    public PotentialAgent previous() {return previous;}
    public void setNext(PotentialAgent potl) {next = potl;}

    protected abstract void makeIterator();
    
    public abstract double energy(IteratorDirective id);

    //PotentialAgent.Hard
    public interface Hard {
        public double energy(IteratorDirective id);
        public void findCollisions(IteratorDirective id, IntegratorHardAbstract.CollisionHandler c);
        public void bump(IntegratorHardAbstract.Agent agent);
    }//end of Hard
    
    /**
     * No-op implementation of the Hard interface.  Methods have no action, energy returns zero.
     */
    public static final Hard HARD_NULL = new Hard() {
        public double energy(IteratorDirective id) {return 0.0;}
        public void findCollisions(IteratorDirective id, IntegratorHardAbstract.CollisionHandler c) {}
        public void bump(IntegratorHardAbstract.Agent agent) {}
    };
        
}//end of Agent
