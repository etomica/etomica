package etomica.data.meter;
import etomica.EtomicaInfo;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataTensor;
import etomica.phase.Phase;
import etomica.potential.PotentialCalculationPressureTensor;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.units.Pressure;

/**
 * Meter for evaluation of the soft-potential pressure tensor in a phase.  This
 * should only be used when using an MC Integrator.  For MD, use a
 * MeterPressureTensorFromIntegrator.
 *
 * @author Andrew Schultz
 */
public class MeterPressureTensor implements DataSource, java.io.Serializable {
    
    public MeterPressureTensor(PotentialMaster potentialMaster) {
    	super();
        this.potentialMaster = potentialMaster;
        Space space = potentialMaster.getSpace();
        data = new DataTensor(space);
        tag = new DataTag();
        dataInfo = new DataTensor.DataInfoTensor("Pressure",Pressure.dimension(space.D()), space);
        dataInfo.addTag(tag);
        rD = 1.0/space.D();
        iteratorDirective = new IteratorDirective();
        iteratorDirective.includeLrc = true;
        pc = new PotentialCalculationPressureTensor(space);
    }
      
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Total pressure in a phase (requires soft-potential model)");
        return info;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    /**
     * Sets the integrator associated with this instance.  The pressure is 
     * calculated for the phase the integrator acts on and integrator's 
     * temperature is used for the ideal gas contribution.
     */
    public void setPhase(Phase newPhase) {
        pc.setPhase(newPhase);
        phase = newPhase;
    }
    
    /**
     * Returns the integrator associated with this instance.  The pressure is 
     * calculated for the phase the integrator acts on and integrator's 
     * temperature is used for the ideal gas contribution.
     */
    public Phase getPhase() {
        return phase;
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
    public Data getData() {
        if (phase == null) {
            throw new IllegalStateException("You must call setPhase before using this class");
        }
    	pc.zeroSum();
        potentialMaster.calculate(phase, iteratorDirective, pc);
        data.x.E(pc.getPressureTensor());
        data.x.TE(1/phase.getBoundary().volume());
        return data;
    }

    private static final long serialVersionUID = 1L;
    protected final DataTag tag;
    protected final DataTensor data;
    protected final DataInfo dataInfo;
    protected final PotentialMaster potentialMaster;
    protected Phase phase;
    protected IteratorDirective iteratorDirective;
    protected final PotentialCalculationPressureTensor pc;
    protected final double rD;
}
