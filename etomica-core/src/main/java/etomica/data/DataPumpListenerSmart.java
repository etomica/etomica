/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.data.history.HistoryCollapsingDiscard;

/**
 * Works with an instance of HistoryCollapsingDiscard to only retrieve new data
 * from the datasource when the history will actually use the data.
 * 
 * @author Andrew Schultz
 */
public class DataPumpListenerSmart extends DataPumpListener {

    protected final HistoryCollapsingDiscard history;
    protected IData lastData;
    
    public DataPumpListenerSmart(IEtomicaDataSource dataSource, IDataSink dataSink, HistoryCollapsingDiscard historySink) {
        this(dataSource, dataSink, 1, historySink);
    }
    
    public DataPumpListenerSmart(IEtomicaDataSource dataSource, IDataSink dataSink, int interval, HistoryCollapsingDiscard historySink) {
        super(dataSource, dataSink, interval);
        this.history = historySink;
    }

    public void actionPerformed() {
        if (lastData == null || !history.willDiscardNextData()) {
            lastData = dataSource.getData();
            if (dataSourceInfo != dataSource.getDataInfo()) {
                dataSourceInfo = dataSource.getDataInfo();
                if (dataSink != null) {
                    dataSink.putDataInfo(dataSourceInfo);
                }
            }
        }
        putData(lastData);
    }
}
