package etomica; 

/**
 * Potential acting on a single atom.
 *
 * @author David Kofke
 */
public abstract class Potential1 extends Potential {
  
    public static String VERSION = "Potential1:01.06.27/"+Potential.VERSION;
    
    public Potential1(Simulation sim) {
        super(sim);
    }
    
    /**
     * Returns the energy of the given atom.
     */
    public abstract double energy(Atom atom);
          
    public PotentialAgent makeAgent(Phase p) {return new Agent(this, p);}
    
   
    //***************** end of methods for Potential1 class *****************//
    
    //Potential1.Agent
    public class Agent extends PotentialAgent {
        
        protected AtomIterator iterator;
        protected Potential1 parentPotential1;
        
        /**
         * @param potential The parent potential making this agent
         * @param phase The phase in which this agent will be placed
         */
        public Agent(Potential potential, Phase phase) {
            super(potential, phase);
            parentPotential1 = (Potential1)potential;
        }
        
        /**
         * Default iterator gives no atoms.
         */
        protected void makeDefaultIterator() {
            iterator = AtomIterator.NULL;
        }
            
        public void setIterator(AtomIterator iterator) {
            this.iterator = iterator;
        }
        public AtomIterator iterator() {return iterator;}
    
        public void calculate(IteratorDirective id, PotentialCalculation pc) {
            if( !(pc instanceof Potential1Calculation) ) return;
            iterator.reset(id);
            ((Potential1Calculation)pc).calculate(iterator, parentPotential1); 
        }
        
    }//end of Agent    
}//end of Potential1



