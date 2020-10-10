package etomica.box.storage;

public final class ObjectStorage<V> extends AbstractObjectStorage<V> {
    private final Factory<V> factory;

    public ObjectStorage(Factory<V> factory) {
        super(factory.getVClass());
        this.factory = factory;
    }

    @Override
    protected V createObject(int idx) {
        return factory.create(idx);
    }

    @Override
    protected void destroy(V value, int idx) {
        this.factory.destroy(value, idx);
    }

    public interface Factory<V> {
        Class<? extends V> getVClass();

        V create(int idx);

        default void destroy(V value, int idx) {}
    }
}
