package etomica.data;

import etomica.util.Arrays;

public class DataStuff {

    public DataStuff(DataSource dataSource, Object object) {
        this.dataSource = dataSource;
        clients = new Object[]{object};
    }
    
    public DataSource getDataSource() {
        return dataSource;
    }
    
    public void setClients(Object[] objects) {
        clients = (Object[])objects.clone();
    }
    
    public Object[] getClients() {
        return (Object[])clients.clone();
    }
    
    public void removeClient(Object object) {
        clients = Arrays.removeObject(clients,object);
    }
    
    public void addClient(Object object) {
        for (int i=0; i<clients.length; i++) {
            if (clients[i] == object) {
                return;
            }
        }
        clients = Arrays.addObject(clients,object);
    }

    private final DataSource dataSource;
    private Object[] clients;
}
