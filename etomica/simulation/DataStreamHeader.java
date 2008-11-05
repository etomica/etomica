package etomica.simulation;

import etomica.data.IEtomicaDataSource;

public class DataStreamHeader implements java.io.Serializable {

    DataStreamHeader(IEtomicaDataSource dataSource, Object object) {
        this.dataSource = dataSource;
        client = object;
    }
    
    public IEtomicaDataSource getDataSource() {
        return dataSource;
    }
    
    public void setClient(Object object) {
        client = object;
    }
    
    public Object getClient() {
        return client;
    }
    
    private final IEtomicaDataSource dataSource;
    private Object client;
}
