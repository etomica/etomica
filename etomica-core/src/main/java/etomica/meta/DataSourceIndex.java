package etomica.meta;

import etomica.data.IDataSource;

public class DataSourceIndex extends ComponentIndex<IDataSource> {

    public DataSourceIndex(Class<IDataSource> componentClass) {
        super(componentClass);
    }

    public static void main(String[] args) {
        System.out.println(new DataSourceIndex(IDataSource.class).getComponentSet());
    }
}
