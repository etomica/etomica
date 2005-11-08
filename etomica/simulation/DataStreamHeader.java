package etomica.simulation;

import etomica.data.DataSource;

public class DataStreamHeader implements java.io.Serializable {

    DataStreamHeader(DataSource dataSource, Object object) {
        this.dataSource = dataSource;
        client = object;
    }
    
    public DataSource getDataSource() {
        return dataSource;
    }
    
    public void setClient(Object object) {
        client = object;
    }
    
    public Object getClient() {
        return client;
    }
    
    private final DataSource dataSource;
    private Object client;
}
