package etomica.virial;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveEvent;
import etomica.integrator.mcmove.MCMoveListener;
import etomica.integrator.mcmove.MCMovePhase;
import etomica.units.Null;

/**
 * Measures value of clusters in a phase and returns the values
 * divided by the sampling bias of the integrator.
 */

public class MeterVirial implements DataSource, MCMoveListener, java.io.Serializable {

	/**
	 * Constructor for MeterVirial.
	 */
	public MeterVirial(ClusterAbstract[] aClusters) {
		clusters = aClusters;
        data = new DataDoubleArray(clusters.length);
        dataInfo = new DataInfoDoubleArray("Cluster Value",Null.DIMENSION, new int[]{clusters.length});
        oldValues = new double[clusters.length];
        trialValues = new double[clusters.length];
	}

	public DataInfo getDataInfo() {
        return dataInfo;
    }
    
    public void setIntegrator(IntegratorMC newIntegrator) {
        if (integrator != null) {
            integrator.getMoveEventManager().removeListener(this);
        }
        integrator = newIntegrator;
        if (integrator != null) {
            integrator.getMoveEventManager().addListener(this);
        }
    }
    
    public IntegratorMC getIntegrator() {
        return integrator;
    }
    
    public void actionPerformed(MCMoveEvent event) {
        if (!event.isTrialNotify) {
            if (event.wasAccepted) {
                // if a move is accepted make the "trial" values into the "old" values
                for (int i=0; i<clusters.length; i++) {
                    oldValues[i] = trialValues[i];
                }
                oldPi = trialPi;
            }
            // do nothing on notification of rejection
            return;
        }
        
        if (nTrials == 0) {
            //we just returned data, so re-zero
            data.E(0);
        }
        
        // this is a trial
        PhaseCluster phase = (PhaseCluster)((MCMovePhase)event.mcMove).getPhase();
        CoordinatePairSet cPairSet = phase.getCPairSet();
        AtomPairSet aPairSet = phase.getAPairSet();
        trialPi = phase.getSampleCluster().value(cPairSet,aPairSet);
        double x[] = data.getData();
        for (int i=0; i<clusters.length; i++) {
            trialValues[i] = clusters[i].value(cPairSet,aPairSet);
            // actually add to the sum on notification of the trial
            x[i] += (oldValues[i]+trialValues[i]) / (oldPi+trialPi);
        }
        nTrials++;
    }
    
    //returns the average gamma/pi since getData was last called
	public Data getData() {
        // if you divide by 0 here, you've called the method twice without taking 
        // any data (via actionPerformed)!
        data.TE(1/(double)nTrials);
        // flag data to get re-zero'd next call to actionPerformed
        nTrials = 0;
        return data;
	}

    protected final ClusterAbstract clusters[];
    private final double oldValues[], trialValues[];
    private double oldPi, trialPi;
    private int nTrials;
    private IntegratorMC integrator;
	private final DataDoubleArray data;
	private final DataInfo dataInfo;
}
