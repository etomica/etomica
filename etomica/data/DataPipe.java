package etomica.data;

import etomica.Data;
import etomica.DataSink;


/**
 * An object that serves as a DataSink, and which may hold one or more
 * other data sinks to which it can push its data.
 *
 */

/*
 * History
 * Created on Feb 19, 2005 by kofke
 */
public abstract class DataPipe extends DataPusher implements DataSink {

    public DataPipe() {
        super();
    }

    public abstract void putData(Data data);

}
