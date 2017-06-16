/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorEvent;

public class DataPumpListener extends DataPump implements IntegratorListener {

    protected long interval;
    protected long intervalCount;
    
    public DataPumpListener(IEtomicaDataSource dataSource, IDataSink dataSink) {
        this(dataSource, dataSink, 1);
    }
    
    public DataPumpListener(IEtomicaDataSource dataSource, IDataSink dataSink, int interval) {
        super(dataSource, dataSink);
        setInterval(interval);
    }

    public void integratorInitialized(IntegratorEvent e) {}
    
    public void integratorStepStarted(IntegratorEvent e) {}
    
    public void integratorStepFinished(IntegratorEvent e) {
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
