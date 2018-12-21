/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data;

import etomica.action.IAction;

/**
 * A DataProcessor whose action is to actively take Data from a DataSource and send it to
 * DataSinks.  
 */
public class DataPump extends DataProcessor implements IAction {

    /**
	 * Constructs DataPump with the given DataSource and
	 * DataSink.  Data source cannot be null.  Data sink can 
     * be null and must be identified via setDataSink if DataPump
     * is to have any effect.
	 */
    public DataPump(IDataSource dataSource, IDataSink dataSink) {
        if(dataSource == null) throw new NullPointerException("Error: cannot construct data pump without a data source");
        this.dataSource = dataSource;
        dataSourceInfo = dataSource.getDataInfo();
        setDataSink(dataSink);
        putDataInfo(dataSource.getDataInfo());
	}
    
	/**
     * Transmits the data from the source to the sink. Before transmitting
     * the Data, this method will first check that the DataInfo from the source
     * is the same as it was last time this method was invoked.  If it has changed,
     * a call to putDataInfo in the sink will be invoked before passing along the Data.
	 */
	public void actionPerformed() {
        IData data = dataSource.getData();
        System.out.println(data.getValue(0));
        if (dataSourceInfo != dataSource.getDataInfo()) {
            dataSourceInfo = dataSource.getDataInfo();
            if (dataSink != null) {
                dataSink.putDataInfo(dataSourceInfo);
            }
        }
        putData(data);
    }
    
    /**
     * Returns the given Data.
     */
    public IData processData(IData inputData) {
        return inputData;
    }
    
    /**
     * Returns the given DataInfo.
     */
    public IDataInfo processDataInfo(IDataInfo inputDataInfo) {
        dataInfo = inputDataInfo.getFactory().makeDataInfo();
        dataInfo.addTag(tag);
        return dataInfo;
    }

    /**
     * @return Returns the dataSource.
     */
    public IDataSource getDataSource() {
        return dataSource;
    }

    protected IDataInfo dataSourceInfo;
    protected final IDataSource dataSource;
}
