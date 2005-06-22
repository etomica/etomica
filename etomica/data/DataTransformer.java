/*
 * History
 * Created on Jul 28, 2004 by kofke
 */
package etomica.data;

import etomica.Data;
import etomica.DataSink;

/**
 * Data pipe that processes given data and immediately transmits
 * processed data to other DataSink(s).
 */
public abstract class DataTransformer extends DataPipe {

    /**
     * Constructs DataTransformer with no DataSinks
     */
    public DataTransformer() {
        super();
    }

    /**
     * Constructs DataTransformer with the given DataSinks
     */
	public DataTransformer(DataSink[] dataSinks) {
        this();
		setDataSinks(dataSinks);
	}
    
    /**
     * Constructs DataTransformer with the a single DataSink.
     */
    public DataTransformer(DataSink dataSink) {
        this(new DataSink[] {dataSink});
    }
	
	/* (non-Javadoc)
	 * @see etomica.Integrator.IntervalListener#intervalAction(etomica.Integrator.IntervalEvent)
	 */
	public void putData(Data newData) {
        pushData(transformData(newData));
    }

    /**
     * Defines the transformation performed by this class.  Input and 
     * output data need not be the same length.
     */
    protected abstract Data transformData(Data data);

}
