package etomica;

/**
 * Methods for properties obtained for a soft, differentiable pair potential.
 *
 * @author David Kofke
 */

public abstract class Potential2Soft extends Potential2 {
    
    
    public Potential2Soft(PotentialGroup parent) {
        super(parent);
    }
    public Potential2Soft(PotentialGroup parent, PotentialTruncation trunc) {
        super(parent, trunc);
    }
    
    public abstract double virial(AtomPair pair);
    
    public abstract double hyperVirial(AtomPair pair);
    
    public abstract Space.Vector gradient(AtomPair pair);
    
    /**
     * Integral used to evaluate correction to truncation of potential.
     */
    public abstract double integral(double rC);
    
    
    public final void calculate2(IteratorDirective id, Potential2Calculation pc) {
//        if( !(pc instanceof Potential2Calculation) ) return;
        //at this point we have identified a (pair)iterator and a (pair)potential
        //that are compatible, in that the iterates form correct arguments to the potential;
        //but the PotentialCalculation class must handle arbitrary iterator/potential sets 
        //it calls the methods of the potential using the iterates as arguments, but it
        //doesn't know that the iterates are atomPairs.  We want to avoid casting each
        //iterate to atomPair, so we need to define separate methods for each iterator/potential set
        //(atom, atomPair, atom3, and so on)
//        iterator.reset(id);
            
        //do we give the iterator to the calculation...
//        ((Potential2Calculation)pc).calculate(iterator, this); 
        pc.calculate(iterator,this);
            //inconvenience:  must put loop construct in every PotentialCalculation
            //problem:  must have different calculate methods for each Potentialx type

        //...or the calculation to the iterator
  ///      iterator.allPairs(((PotentialCalculationEnergySum)pc).getAtomPairCalculation(this)); 
 //       System.out.println();
            //problem:  cannot abort iteration once started (e.g., overlap detected)
            //          but this could be resolved by adding a boolean abort() method to Calculation to end iteration
            //problem: setPotential would have to be defined separately for each Potentialx type,
            //         and different copies of the potential kept in the PotentialCalculation object
    }//end of calculate

}//end of Potential2Soft