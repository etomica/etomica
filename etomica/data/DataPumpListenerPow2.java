package etomica.data;

import etomica.api.IIntegratorEvent;
import etomica.api.IIntegratorListener;

/**
 * DataPump class that takes data at steps 2^i, i=0,1,2,3...
 * This can be used to take (and retain) data over many orders of magnitude.
 * 
 * @author Andrew Schultz
 */
public class DataPumpListenerPow2 extends DataPump implements IIntegratorListener {

    private static final long serialVersionUID = 1L;
    private long interval;
    private long intervalCount;
    
    public DataPumpListenerPow2(IEtomicaDataSource dataSource, IDataSink dataSink) {
        super(dataSource, dataSink);
        reset();
    }

    public void integratorInitialized(IIntegratorEvent e) {}
    
    public void integratorStepStarted(IIntegratorEvent e) {}
    
    public void integratorStepFinished(IIntegratorEvent e) {
        if(++intervalCount < interval) return;
        if (interval == 0) {
            interval = 1;
        }
        else {
            interval *= 2;
        }
        intervalCount = 0;
        actionPerformed();
    }
    
    public void reset() {
        interval = 0;
    }
    
    public long getInterval() {
        return interval;
    }
}
