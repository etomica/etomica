/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;

/**
 * DataPump class that takes data at steps pow^i, i=0,1,2,3...
 * This can be used to take (and retain) data over many orders of magnitude.
 * Pow is 2 by default.
 * 
 * @author Andrew Schultz
 */
public class DataPumpListenerPowStrict extends DataPump implements IntegratorListener {

    protected long interval;
    protected long intervalCount;
    protected int pow;

    public DataPumpListenerPowStrict(IDataSource dataSource, IDataSink dataSink) {
        this(dataSource, dataSink, 2);
    }

    public DataPumpListenerPowStrict(IDataSource dataSource, IDataSink dataSink, int pow) {
        super(dataSource, dataSink);
        this.pow = pow;
        reset();
    }
    
    public void integratorInitialized(IntegratorEvent e) {}
    
    public void integratorStepStarted(IntegratorEvent e) {}
    
    public void integratorStepFinished(IntegratorEvent e) {
        if(++intervalCount < interval) {
            return;
        }
        if (interval == 0) {
            interval = 1;
        }
        else {
            interval *= pow;
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
