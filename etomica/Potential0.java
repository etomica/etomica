package etomica; 

/**
 * Potential that does not depend on any atom positions.
 * Typically used to implement long-range corrections for potential truncation.
 * Potential thus depends on phase parameters, such as the number of molecules and the volume.
 *
 * @author David Kofke
 */
 
 //under development
public abstract class Potential0 extends Potential {
  
    public static String VERSION = "Potential0:01.07.08/"+Potential.VERSION;
    
    public Potential0(Simulation sim) {
        super(sim);
    }
              
    public PotentialAgent makeAgent(Phase p) {return new Agent(this, p);}
    
   
    //***************** end of methods for Potential0 class *****************//
    
    //Potential1.Agent
    public class Agent extends PotentialAgent {
        
        /**
         * @param potential The parent potential making this agent
         * @param phase The phase in which this agent will be placed
         */
        public Agent(Potential potential, Phase phase) {
            super(potential, phase);
            parentPotential1 = (Potential1)potential;
        }
        
        /**
         * No iterator is needed; this method has no action.
         */
        protected void makeDefaultIterator() { }
            
        public void calculate(IteratorDirective id, PotentialCalculation pc) {
   //         if( !(pc instanceof Potential0Calculation) ) return;
            ((Potential0Calculation)pc).calculate(parentPotential1); 
        }
        
    }//end of Agent    
}//end of Potential0



