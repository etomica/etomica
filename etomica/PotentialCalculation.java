package etomica;

public interface PotentialCalculation {
    
    public interface Sum extends PotentialCalculation {
        public double sum();
    }
    
    /**
     * Evaluates the energy summed over all iterated atoms.
     */
    public static final class EnergySum implements PotentialCalculation.Sum, Potential1Calculation, Potential2Calculation {
        
        private double sum = 0.0;
        
        public EnergySum() {}
        
        public PotentialCalculation.Sum reset() {sum = 0.0; return this;}
        public double sum() {return sum;}
        
        //atom
        public void calculate(AtomIterator iterator, Potential1 potential) {
            while(iterator.hasNext()) {
                sum += potential.energy(iterator.next());
                if(sum >= Double.MAX_VALUE) return;
            }//end while
        }//end of calculate

        //pair
        public void calculate(AtomPairIterator iterator, Potential2 potential) {
            while(iterator.hasNext()) {
                sum += potential.energy(iterator.next());
                if(sum >= Double.MAX_VALUE) return;
            }//end while
        }//end of calculate
    }//end EnergySum

    /**
     * Sums the force on each iterated atom and adds it to the integrator agent
     * associated with the atom.
     */
    public static final class ForceSum implements Potential1Calculation, Potential2Calculation {
        
        private final Space.Vector f;
        public ForceSum(Space space) {
             f = space.makeVector();
        }
        
        //atom
        public void calculate(AtomIterator iterator, Potential1 potential) {
            Potential1Soft potentialSoft = (Potential1Soft)potential;
            while(iterator.hasNext()) {
                Atom atom = iterator.next();
                f.E(potentialSoft.gradient(atom));
                ((Integrator.Agent.Forcible)atom.ia).force().ME(f);
            }//end while
        }//end of calculate

        //pair
        public void calculate(AtomPairIterator iterator, Potential2 potential) {
            Potential2Soft potentialSoft = (Potential2Soft)potential;
            while(iterator.hasNext()) {
                AtomPair pair = iterator.next();
                f.E(potentialSoft.gradient(pair));
                ((Integrator.Agent.Forcible)pair.atom1().ia).force().PE(f);
                ((Integrator.Agent.Forcible)pair.atom2().ia).force().ME(f);
            }//end while
        }//end of calculate
    }//end ForceSums
    
}//end of PotentialCalculation