package etomica.meta.properties;

/**
 * Very simple property descriptor for an element in a Collection. Allows read of element but nothing else.
 */
public class CollectionElementProperty implements Property {

    private final Object instance;

    public CollectionElementProperty(Object instance) {
        this.instance = instance;
    }


    public Object invokeReader() {
        return instance;
    }

    public String getName() {
        return "dummy";
    }

    public Object invokeReader(int i) {
        return null;
    }

    public void invokeWriter(Object... params) {
        throw new RuntimeException();
    }

    public void invokeWriter(int i, Object... params) {
        throw new RuntimeException();
    }

    public void invokeAdder(Object o) {
        throw new RuntimeException();
    }

    public void invokeRemover(Object o) {
        throw new RuntimeException();
    }

    public int invokeCount() {
        return 0;
    }

    public boolean isIndexedProperty() {
        return false;
    }

    public boolean canRead() {
        return true;
    }

    public boolean canWrite() {
        return false;
    }

    public boolean canAdd() {
        return false;
    }

    public boolean canRemove() {
        return false;
    }

    public boolean canCount() {
        return false;
    }

    public Class<?> getPropertyType() {
        return instance.getClass();
    }
}
