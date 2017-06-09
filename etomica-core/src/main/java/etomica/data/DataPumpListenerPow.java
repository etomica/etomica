/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorEvent;

/**
 * DataPump class that takes data at steps 2^i, i=0,1,2,3...
 * This can be used to take (and retain) data over many orders of magnitude.
 * 
 * With pow set, the pump behaves differently.  For instance, pow=10
 * 
 * 1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400
 * 
 * @author Andrew Schultz
 */
public class DataPumpListenerPow extends DataPump implements IntegratorListener {

    protected long interval;
    protected long stepCount;
    protected long intervals;
    protected long totalSteps;
    protected int maxRatio;
    protected int pow;
    protected IData lastData;
    
    public DataPumpListenerPow(IEtomicaDataSource dataSource, IDataSink dataSink) {
        this(dataSource, dataSink, 0, 10, 100);
    }

    public DataPumpListenerPow(IEtomicaDataSource dataSource, IDataSink dataSink, long interval, int pow, int maxRatio) {
        super(dataSource, dataSink);
        this.pow = pow;
        this.interval = interval;
        this.maxRatio = maxRatio;
    }
    
    public void integratorInitialized(IntegratorEvent e) {}
    
    public void integratorStepStarted(IntegratorEvent e) {}
    
    public void integratorStepFinished(IntegratorEvent e) {
        totalSteps++;
        if(++stepCount < interval) {
            return;
        }
        if (interval == 0) {
            interval = 1;
            intervals = 1;
        }
        else if (totalSteps/interval > maxRatio) {
            interval *= pow;
            intervals = 1;
        }
        else {
            intervals++;
        }
        stepCount = 0;
        actionPerformed();
    }
    
    public void actionPerformed() {
        lastData = dataSource.getData();
        if (dataSourceInfo != dataSource.getDataInfo()) {
            dataSourceInfo = dataSource.getDataInfo();
            if (dataSink != null) {
                dataSink.putDataInfo(dataSourceInfo);
            }
        }
        putData(lastData);
    }
    
    public void reset() {
        totalSteps = 0;
        interval = 0;
    }
    
    public long getInterval() {
        return interval;
    }
}
