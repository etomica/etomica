package etomica;

public abstract class PotentialAgent implements java.io.Serializable {

    private Phase parentPhase;  //want final, but won't compile
    private PotentialAgent next;
    
    /**
        * @param p The phase in which this agent will be placed
        */
    public PotentialAgent(Phase p) {
        parentPhase = p;
    }
    public final Phase parentPhase() {return parentPhase;}
        
    public PotentialAgent next() {return next;}
//    public PotentialAgent previous() {return previous;}
    public void setNext(PotentialAgent potl) {next = potl;}

    public abstract double energy(IteratorDirective id);

    //PotentialAgent.Hard
    public interface Hard {
        public double energy(IteratorDirective id);
        public void findCollisions(IteratorDirective id, IntegratorHardAbstract.CollisionHandler c);
        public void bump(IntegratorHardAbstract.Agent agent);
    }//end of Hard
}//end of Agent
        
    
