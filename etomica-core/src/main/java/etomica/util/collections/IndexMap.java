package etomica.util.collections;

import java.util.Arrays;
import java.util.Collection;
import java.util.Map;
import java.util.Set;

/**
 * An specialized implementation of Map which is designed for mapping
 * an index to a value. It is backed by an array but allows inserting out
 * of index order without having to manage array resizing manually.
 * <p>
 * This should not be used as a general int-keyed map as the backing array will take space up to the
 * largest key. Use this when you have an expected maximum index (but not necessarily known
 * or constant) that you need to sporadically fill mappings up to.
 *
 * @param <V> the type of value contained in this map.
 */
public class IndexMap<V> implements Map<Integer, V> {
    private final int reservoirSize;
    private V[] values;
    private int size;
    private boolean isEmpty;

    /**
     * Creates a new map with the default initial capacity of 10 and
     * default reservoir size of 30.
     */
    public IndexMap() {
        this(10, 30);
    }

    /**
     * Creates a new map with the given initial capacity and default
     * reservoir size of 30.
     *
     * @param initialCapacity the base capacity to start with.
     */
    public IndexMap(int initialCapacity) {
        this(initialCapacity, 30);
    }

    /**
     * Creates a new map with the given initial capacity and reservoir size.
     * The backing array will be of size initialCapacity + reservoirSize, and
     * keys up to that value can be added without the array being resized.
     *
     * @param initialCapacity the base capacity to start with.
     * @param reservoirSize   the amount of padding to be included in the backing array.
     */
    @SuppressWarnings("unchecked")
    public IndexMap(int initialCapacity, int reservoirSize) {
        this.values = (V[]) new Object[initialCapacity + reservoirSize];
        this.size = 0;
        this.isEmpty = true;
        this.reservoirSize = reservoirSize;
    }

    @Override
    public int size() {
        return this.size;
    }

    @Override
    public boolean isEmpty() {
        return this.isEmpty;
    }

    @Override
    public boolean containsKey(Object key) {
        return key instanceof Number && this.containsKey(((Number) key).intValue());
    }

    @Override
    public boolean containsValue(Object value) {
        if(value == null) {
            return false;
        } else {
            for(V v : values) {
                if(v != null && v.equals(value)) {
                    return true;
                }
            }
        }

        return false;
    }

    @Override
    public V get(Object key) {
        if(key instanceof Number) {
            return this.get(((Number) key).intValue());
        } else {
            return null;
        }
    }

    @Override
    public V put(Integer key, V value) {
        return this.put(key.intValue(), value);
    }

    @Override
    public V remove(Object key) {
        if(key instanceof Number) {
            return this.remove(((Number) key).intValue());
        } else {
            return null;
        }
    }

    @Override
    public void putAll(Map<? extends Integer, ? extends V> m) {
        m.forEach(this::put);
    }

    @Override
    @SuppressWarnings("unchecked")
    public void clear() {
        this.values = (V[]) new Object[1 + reservoirSize];
        this.size = 0;
        this.isEmpty = true;
    }

    @Override
    public Set<Integer> keySet() {
        throw new UnsupportedOperationException();
    }

    @Override
    public Collection<V> values() {
        throw new UnsupportedOperationException();
    }

    @Override
    public Set<Entry<Integer, V>> entrySet() {
        throw new UnsupportedOperationException();
    }

    public boolean containsKey(int key) {
        return key < this.values.length && key >= 0 && this.values[key] != null;
    }

    public V get(int key) {
        if(key < this.values.length && key >= 0) {
            return values[key];
        } else {
            return null;
        }
    }

    public V put(int key, V value) {
        if(value == null) {
            throw new IllegalArgumentException("Value must not be null");
        }

        if(key >= this.values.length) {
            this.values = Arrays.copyOf(this.values, key + reservoirSize);
        }

        this.values[key] = value;
        this.isEmpty = false;
        this.size++;
        return value;
    }

    public V remove(int key) {
        if(key < this.values.length && key >= 0) {
            V removed = this.values[key];
            if(removed != null) {
                this.values[key] = null;
                this.size--;
                if(this.size == 0) {
                    this.isEmpty = true;
                }
                return removed;
            }
        }
        return null;
    }
}
