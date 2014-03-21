package etomica.data;

import etomica.api.IIntegratorEvent;
import etomica.api.IIntegratorListener;

public class DataPumpListener extends DataPump implements IIntegratorListener {

    protected long interval;
    protected long intervalCount;
    
    public DataPumpListener(IEtomicaDataSource dataSource, IDataSink dataSink) {
        this(dataSource, dataSink, 1);
    }
    
    public DataPumpListener(IEtomicaDataSource dataSource, IDataSink dataSink, int interval) {
        super(dataSource, dataSink);
        setInterval(interval);
    }

    public void integratorInitialized(IIntegratorEvent e) {}
    
    public void integratorStepStarted(IIntegratorEvent e) {}
    
    public void integratorStepFinished(IIntegratorEvent e) {
        if(++intervalCount < interval) return;
        intervalCount = 0;
        actionPerformed();
    }
    
    public void setInterval(long i) {
        interval = i;
    }
    
    public long getInterval() {
        return interval;
    }
}
