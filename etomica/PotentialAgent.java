package etomica;

/**
 * Representative of a potential in a phase.  Each potential constructs
 * an agent specific to it (and ususally defined as an inner class within it), and
 * places that agent in every phase during initiation of the simulation.
 * One of the primary responsibilities of the agent is to define the iterator
 * that loops over the atoms affected by the potential.  Agents from the same 
 * potential but in different phases have in common the parent potential, which
 * holds the common parameters (e.g., diameter, well depth) referenced by the 
 * potential functions.  Any request for properties of the potential (e.g. a value 
 * for the energy) are directed to the phase and are routed to the appropriate
 * potential agents in the phase.
 *
 * @author David Kofke
 */
public abstract class PotentialAgent implements java.io.Serializable {

    protected Phase parentPhase;  //want final, but won't compile
    protected Potential parentPotential;
    private PotentialAgent next;
    
    /**
     * @param p The phase in which this agent will be placed
     */
    public PotentialAgent(Potential potential, Phase phase) {
        parentPotential = potential;
        parentPhase = phase;
        parentPotential.agents().put(phase, this);
        makeDefaultIterator();
    }
    public final Phase parentPhase() {return parentPhase;}
    public final Potential parentPotential() {return parentPotential;}
        
    public PotentialAgent nextAgent() {return next;}
    public void setNextAgent(PotentialAgent potl) {next = potl;}

    protected abstract void makeDefaultIterator();
    
    public abstract void calculate(IteratorDirective id, PotentialCalculation pa);
    
    //PotentialAgent.Hard
 /*   public interface Hard {
        public double energy(IteratorDirective id);
        public void findCollisions(IteratorDirective id, IntegratorHardAbstract.CollisionHandler c);
        public void bump(IntegratorHardAbstract.Agent agent);
        
        public static final class Linker implements java.io.Serializable {
            public final PotentialAgent.Hard agent;
            public Linker next = null;
            //Constructors
            public Linker(Hard a) {agent = a;}
            public Linker(Hard a, Linker l) {agent = a; next = l;}
            //Access methods
 //           public PotentialAgent.Hard agent() {return agent;}
            public Linker next() {return next;}
            public void setNext(Linker l) {next = l;}
        }//end of Linker
    }//end of Hard
    
    /**
     * No-op implementation of the Hard interface.  Methods have no action, energy returns zero.
     * /
    public static final Hard HARD_NULL = new Hard() {
        public double energy(IteratorDirective id) {return 0.0;}
        public void findCollisions(IteratorDirective id, IntegratorHardAbstract.CollisionHandler c) {}
        public void bump(IntegratorHardAbstract.Agent agent) {}
    };
   */     
}//end of Agent
