package etomica.data.meter;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorBox;
import etomica.potential.PotentialCalculationHyperVirialSum;
import etomica.units.Null;

/**
 * Meter for measurement of the hypervirial, r^2 d2u/dr2 + r du/dr
 */
public class MeterHyperVirial extends DataSourceScalar {

    private IntegratorBox integrator;
    private IteratorDirective iteratorDirective;
    protected final PotentialCalculationHyperVirialSum hyperVirial;
    protected IPotentialMaster potentialMaster;
    protected IBox box;

    public MeterHyperVirial() {
        super("hypervirial", Null.DIMENSION);
        iteratorDirective = new IteratorDirective();
        iteratorDirective.includeLrc = true;
        hyperVirial = new PotentialCalculationHyperVirialSum();
    }

    public void setPotentialMaster(IPotentialMaster newPotentialMaster) {
        potentialMaster = newPotentialMaster;
    }
    
    public void setBox(IBox newBox) {
        box = newBox;
    }

    /**
     * Returns the integrator associated with this instance.  The pressure is 
     * calculated for the box the integrator acts on and integrator's 
     * temperature is used for the ideal gas contribution.
     */
    public IntegratorBox getIntegrator() {
        return integrator;
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
	  * Computes total pressure in box by summing virial over all pairs, and adding
	  * ideal-gas contribution.
	  */
    public double getDataAsScalar() {
        if (potentialMaster == null || box == null) {
            throw new IllegalStateException("You must call setIntegrator before using this class");
        }
        hyperVirial.zeroSum();
        potentialMaster.calculate(box, iteratorDirective, hyperVirial);
        return hyperVirial.getSum();
    }
}
