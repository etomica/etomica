package etomica;
import etomica.performance.*;

/**
 * Evaluates the energy summed over all iterated atoms.
 */
public class PotentialCalculationEnergySumPerformance extends PotentialCalculationEnergySum {
        
    protected double sum = 0.0;
        
    public PotentialCalculationEnergySumPerformance() {}
        
    public PotentialCalculation.Sum reset() {sum = 0.0; return this;}
    public double sum() {;return sum;}
    
    public void setSum(double s){sum=s;}
        
    public AtomPairAction getAtomPairCalculation(PotentialBase potential) {
        atomPairAction.potential = (PotentialBase)potential;
        return atomPairAction;
    }//end of getAtomPairCalculation

    private final MyAtomPairAction atomPairAction = new MyAtomPairAction();
    public final class MyAtomPairAction extends AtomPairAction {
        PotentialBase potential;
        public void action(AtomPair pair) {}
        public void action(int i, int j) {
            sum += potential.energy(i,j);
        }
    }
}//end EnergySum
