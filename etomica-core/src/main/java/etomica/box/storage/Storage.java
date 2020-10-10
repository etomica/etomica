package etomica.box.storage;

public interface Storage {
    void add(int n);

    void addNull(int n);

    void ensureCapacity(int expectedAdditions);

    int size();

    void swapRemove(int i);
}
