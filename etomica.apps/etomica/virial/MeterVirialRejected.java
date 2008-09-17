package etomica.virial;

import etomica.api.IEvent;
import etomica.api.IListener;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.DataTag;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.integrator.mcmove.MCMoveEvent;
import etomica.integrator.mcmove.MCMoveTrialCompletedEvent;
import etomica.units.Null;

/**
 * Measures value of clusters in a Box and returns the value of the system's
 * cluster divided by the sampling bias in addition to the Bennett overlap
 * integrands using appropriate values of alpha.
 * 
 * This meter uses information from rejected moves instead of using cluster
 * values after an MC move.  Data returned is equal to
 * 
 * gamma = (gamma1 + gamma2) / (pi1 + pi2)
 * 
 * for the system's cluster.  An equivalent formula is used for the overlap
 * clusters.
 */
public class MeterVirialRejected implements DataSource, IListener, java.io.Serializable {

	/**
	 * Constructor for MeterVirialRejected.
	 */
	public MeterVirialRejected(ClusterAbstract[] aClusters, int aNBennetPoints, boolean aIsReference) {
        nBennetPoints = aNBennetPoints;
        isReference = aIsReference;
        if (nBennetPoints%2 == 0) {
            throw new IllegalArgumentException("try again with an odd aNPoints");
        }
		clusters = aClusters;
        data = new DataDoubleArray(nBennetPoints+1);
        dataInfo = new DataInfoDoubleArray("Cluster Value",Null.DIMENSION, new int[]{nBennetPoints+1});
        oldValues = new double[nBennetPoints+1];
        trialValues = new double[nBennetPoints+1];
        tag = new DataTag();
        dataInfo.addTag(tag);
        alpha = new double[nBennetPoints];
        setBennetParam(1.0,5);
	}
	
	public ClusterAbstract[] getClusters() {
	    return clusters;
	}

    /**
     * sets the range of parameter values used for Bennets method.
     * Default is a span of 5 centered about 1 (exp(-5) to (exp(5)).
     * @param aCenter geometric mean of all values
     * @param aSpan natural log of ratio of max value to aCenter
     */
    public void setBennetParam(double aCenter, double aSpan) {
        if (aSpan < 0.0 || (aSpan == 0 && nBennetPoints > 1) || aCenter <= 0.0 ) throw new IllegalArgumentException("span and center must be positive");
        if (nBennetPoints==1) {
            alpha[0] = aCenter;
            return;
        }
        for (int i=0; i<nBennetPoints; i++) {
            alpha[i] = Math.exp(2.0*aSpan*(i-(nBennetPoints-1)/2)/(nBennetPoints-1))*aCenter;
        }
        // our oldValues are bogus.
        oldPi = 0;
    }

    public int getNBennetPoints() {
        return nBennetPoints;
    }
    
    /**
     * Returns iParam'th factor used in the Bennet sum.  Higher values
     * indicate the overlap function is more like the 0th value.
     */
    public double getBennetBias(int iParam) {
        return alpha[iParam];
    }
    
	public DataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
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
    
    public void actionPerformed(IEvent event) {
        if (event instanceof MCMoveTrialCompletedEvent) {
            if (((MCMoveTrialCompletedEvent)event).isAccepted()) {
                // if a move is accepted make the "trial" values into the "old" values
                for (int i=0; i<trialValues.length; i++) {
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
        BoxCluster box = (BoxCluster)((MCMoveBox)((MCMoveEvent)event).getMCMove()).getBox();
        trialPi = box.getSampleCluster().value(box);
        trialValues[0] = clusters[0].value(box);
        double x[] = data.getData();
        x[0] = (oldValues[0]+trialValues[0]) / (oldPi+trialPi);
        double value1 = clusters[1].value(box);
        for (int i=0; i<nBennetPoints; i++) {
            trialValues[i+1] = trialPi * value1;
            if (isReference) {
                trialValues[i+1] /= (value1 + alpha[i]*trialPi);
            }
            else {
                trialValues[i+1] /= (trialPi + alpha[i]*value1);
            }
        }
        for (int i=1; i<x.length; i++) {
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
    private final DataTag tag;
    private final int nBennetPoints;
    private final double[] alpha;
    private final boolean isReference;
}
