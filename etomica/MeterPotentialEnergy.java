package etomica;
import etomica.units.Dimension;

/**
 * Meter for evaluation of the potential energy in a phase
 * Includes several related methods for computing the potential energy of a single
 * atom or molecule with all neighboring atoms
 *
 * @author David Kofke
 */
 
public class MeterPotentialEnergy extends Meter implements EtomicaElement {
    
    private IteratorDirective iteratorDirective;
    private final PotentialCalculation.EnergySum energy = new PotentialCalculation.EnergySum();
    private final PotentialMaster potential;
    
    public MeterPotentialEnergy() {
        this(Simulation.instance);
    }
    public MeterPotentialEnergy(Simulation sim) {
        super(sim);
        setLabel("Potential Energy");
        iteratorDirective = new IteratorDirective();
        potential = sim.hamiltonian.potential;
    }
      
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Total intermolecular potential energy in a phase");
        return info;
    }

    public Dimension getDimension() {return Dimension.ENERGY;}
    
    /**
     * Iterator directive specifies the atoms for which the potential energy is measured.
     * Default value measures all atoms in phase.
     */
    public void setIteratorDirective(IteratorDirective directive) {iteratorDirective = directive;}
    
    /**
     * Accessor method for iterator directive.
     */
    public IteratorDirective getIteratorDirective() {return iteratorDirective;}
      
 /**
  * Computes total potential energy for all atom pairs in phase.
  * Currently, does not include long-range correction to truncation of energy
  */
    public final double currentValue() {
        return potential.set(phase).calculate(iteratorDirective, energy.reset()).sum();
    }
    
}//end of MeterPotentialEnergy
