package etomica.meta.wrappers;


import java.util.List;

public abstract class CollectionWrapper<T> extends Wrapper<T> {

    public CollectionWrapper(T wrapped) {
        super(wrapped);
        this.properties.clear();
    }

    public abstract List<Wrapper<?>> getElements();
}
