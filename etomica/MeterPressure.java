package etomica;
import etomica.units.Dimension;

/**
 * Meter for evaluation of the soft-potential pressure in a phase.
 *
 * @author David Kofke
 */
 
public class MeterPressure extends MeterScalar implements EtomicaElement {
    
    private IteratorDirective iteratorDirective;
    private final PotentialCalculationVirialSum virial;
    private final PotentialMaster potential;
    private final Integrator integrator;
    private final double rD;
    
/*    public MeterPressure() {
        this(Simulation.instance);
    }*/
    //requires Integrator for temperature
    public MeterPressure(Integrator integrator) {
        super(integrator.parentSimulation());
        this.integrator = integrator;
        rD = 1.0/(double)parentSimulation().space.D();
        setLabel("Pressure");
        iteratorDirective = new IteratorDirective();
        iteratorDirective.includeLrc = true;
        potential = parentSimulation().hamiltonian.potential;
        virial = new PotentialCalculationVirialSum();
    }
      
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Total pressure in a phase (requires soft-potential model)");
        return info;
    }

    public Dimension getDimension() {return Dimension.PRESSURE;}
    
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
  * Computes total pressure in phase by summing virial over all pairs, and adding
  * ideal-gas contribution.
  * Currently, does not include long-range correction to truncation of energy.
  */
    public final double currentValue() {
        double dbv = potential.calculate(phase, iteratorDirective.set(), virial.reset()).sum();
        double p = phase.getDensity()*integrator.temperature() - dbv*rD/phase.boundary().volume();
        return p;
    }
    
}//end of MeterPressure
