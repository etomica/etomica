package etomica;
import etomica.units.Dimension;

/**
 * Meter for evaluation of the potential energy in a phase
 * Includes several related methods for computing the potential energy of a single
 * atom or molecule with all neighboring atoms
 *
 * @author David Kofke
 */
 
public class MeterPotentialEnergy extends Meter implements EtomicaElement
{
    private IteratorDirective all;  
    public MeterPotentialEnergy() {
        this(Simulation.instance);
    }
    public MeterPotentialEnergy(Simulation sim) {
        super(sim);
        setLabel("Potential Energy");
        all = new IteratorDirective();
        all.set();
    }
      
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Total intermolecular potential energy in a phase");
        return info;
    }

    public Dimension getDimension() {return Dimension.ENERGY;}
      
 /**
  * Computes total potential energy for all atom pairs in phase
  * Returns infinity (MAX_VALUE) as soon as overlap is detected
  * Currently, does not include long-range correction to truncation of energy
  */
    public final double currentValue() {
        phase.potential().calculate(all, energy);
        return energy.sum;
    }
    
    private final CalculationEnergy energy = new CalculationEnergy();
    
    private static final class CalculationEnergy implements Potential1Calculation, Potential2Calculation {
        
        double sum = 0.0;
        
        public void reset() {sum = 0.0;}
        
        //atom pair
        public void calculate(AtomPairIterator iterator, Potential2 potential) {
            while(iterator.hasNext()) {
                sum += potential.energy(iterator.next());
            }//end while
        }//end of calculate(AtomPair...
        
        //single atom
        public void calculate(AtomIterator iterator, Potential1 potential) {
            while(iterator.hasNext()) {
                sum += potential.energy(iterator.next());
            }//end while
        }//end of calculate(Atom...
    } //end of collisionHandlerUp

    
}
