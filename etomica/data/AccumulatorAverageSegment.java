package etomica.data;

import etomica.DataSink;
import etomica.DataSource;
import etomica.Integrator;
import etomica.data.AccumulatorAverage.Type;
import etomica.integrator.IntervalActionAdapter;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Feb 20, 2005 by kofke
 */
public class AccumulatorAverageSegment {

    /**
     * 
     */
    public AccumulatorAverageSegment(DataSource source, Integrator integrator, Type[] types, DataSink sink) {
        dataPump = new DataPump(source, (DataSink)null);
        intervalActionAdapter = new IntervalActionAdapter(dataPump, integrator);
        accumulator = new AccumulatorAverage();
        accumulator.setPushInterval(100);
        dataPump.addDataSink(accumulator);
        accumulator.makeDataPusher(types).addDataSink(sink);
    }
   
    /**
     * @return Returns the accumulator.
     */
    public AccumulatorAverage getAccumulator() {
        return accumulator;
    }
    /**
     * @return Returns the dataPump.
     */
    public DataPump getDataPump() {
        return dataPump;
    }
    /**
     * @return Returns the intervalActionAdapter.
     */
    public IntervalActionAdapter getIntervalActionAdapter() {
        return intervalActionAdapter;
    }
    
    private final DataPump dataPump;
    private final AccumulatorAverage accumulator;
    private final IntervalActionAdapter intervalActionAdapter;
}
