package etomica.meta;

import etomica.data.IEtomicaDataSource;

public class DataSourceIndex extends ComponentIndex<IEtomicaDataSource> {

    public DataSourceIndex(Class<IEtomicaDataSource> componentClass) {
        super(componentClass);
    }

    public static void main(String[] args) {
        System.out.println(new DataSourceIndex(IEtomicaDataSource.class).getComponentSet());
    }
}
