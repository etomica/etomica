package etomica; 

/**
 * Potential acting on a pair of atoms.
 *
 * @author David Kofke
 */
public abstract class Potential2 extends Potential {
  
    public static String VERSION = "Potential2:01.07.03/"+Potential.VERSION;
    
    protected Space.Vector work1;
    
    public Potential2(Simulation sim) {
        super(sim);
        work1 = sim.space().makeVector();
    }
    
    /**
     * Returns the energy of the given atom pair.
     */
    public abstract double energy(AtomPair pair);
    
    public PotentialAgent makeAgent(Phase p) {return new Agent(this, p);}
    
            
    //***************** end of methods for Potential2 class *****************//
    
    //Potential2.Agent
    public class Agent extends PotentialAgent {
        
        protected AtomPairIterator iterator;
        protected Potential2 parentPotential2;
        /**
         * @param p The phase in which this agent will be placed
         */
        public Agent(Potential potential, Phase phase) {
            super(potential, phase);
            parentPotential2 = (Potential2)potential;
        }
        
        /**
         * Default iterator yields no pairs.
         */
        protected void makeDefaultIterator() {
            iterator = AtomPairIterator.NULL;
        }
            
        public void setIterator(AtomPairIterator iterator) {
            this.iterator = iterator;
        }
        public AtomPairIterator iterator() {return iterator;}
    
//        public final Potential parentPotential() {return Potential2.this;}
        
        public void calculate(IteratorDirective id, PotentialCalculation pc) {
            if( !(pc instanceof Potential2Calculation) ) return;
            //at this point we have identified a (pair)iterator and a (pair)potential
            //that are compatible, in that the iterates form correct arguments to the potential;
            //but the PotentialCalculation class must handle arbitrary iterator/potential sets 
            //it calls the methods of the potential using the iterates as arguments, but it
            //doesn't know that the iterates are atomPairs.  We want to avoid casting each
            //iterate to atomPair, so we need to define separate methods for each iterator/potential set
            //(atom, atomPair, atom3, and so on)
            iterator.reset(id);
            
            //do we give the iterator to the calculation...
            ((Potential2Calculation)pc).calculate(iterator, parentPotential2); 
                //inconvenience:  must put loop construct in every PotentialCalculation
                //problem:  must have different calculate methods for each Potentialx type

            //...or the calculation to the iterator
            //pc.setPotential(parentPotential);
            //iterator.allPairs(pc); 
                //problem:  cannot abort iteration once started (e.g., overlap detected)
                //          but this could be resolved by adding a boolean abort() method to Calculation to end iteration
                //problem: setPotential would have to be defined separately for each Potentialx type,
                //         and different copies of the potential kept in the PotentialCalculation object
        }//end of calculate
        
    }//end of Agent    

}//end of Potential2



