package etomica.data.meter;
import etomica.EtomicaInfo;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorPhase;
import etomica.phase.Phase;
import etomica.potential.PotentialCalculationVirialSum;
import etomica.space.Space;
import etomica.units.Pressure;

/**
 * Meter for evaluation of the soft-potential pressure in a phase.
 * Requires that temperature be set in order to calculation ideal-gas
 * contribution to pressure; default is to use zero temperature, which
 * causes this contribution to be omitted.
 *
 * @author David Kofke
 */
 
public class MeterPressure extends DataSourceScalar {
    
    public MeterPressure(IntegratorPhase integrator, Space space) {
    	super("Pressure",Pressure.dimension(space.D()));
        rD = 1.0/space.D();
        iteratorDirective = new IteratorDirective();
        iteratorDirective.includeLrc = true;
        this.integrator = integrator;
        virial = new PotentialCalculationVirialSum();
    }
      
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Total pressure in a phase (requires soft-potential model)");
        return info;
    }

    /**
     * Sets flag indicating whether calculated energy should include
     * long-range correction for potential truncation (true) or not (false).
     */
    public void setIncludeLrc(boolean b) {
    	iteratorDirective.includeLrc = b;
    }
    
    /**
     * Indicates whether calculated energy should include
     * long-range correction for potential truncation (true) or not (false).
     */
    public boolean isIncludeLrc() {
    	return iteratorDirective.includeLrc;
    }

	 /**
	  * Computes total pressure in phase by summing virial over all pairs, and adding
	  * ideal-gas contribution.
	  */
    public double getDataAsScalar() {
    	virial.zeroSum();
        Phase phase = integrator.getPhase();
        integrator.getPotential().calculate(phase, iteratorDirective, virial);
        return phase.getDensity()*integrator.getTemperature() - virial.getSum()*rD/phase.getBoundary().volume();
    }

    private IntegratorPhase integrator;
    private IteratorDirective iteratorDirective;
    private final PotentialCalculationVirialSum virial;
    private final double rD;
}
