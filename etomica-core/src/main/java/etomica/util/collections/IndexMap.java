package etomica.util.collections;

import java.util.*;

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
public final class IndexMap<V> implements Map<Integer, V> {
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
        return new IndexMapValuesCollection();
    }

    /**
     * Do NOT use this. It is intended for the IDE debugger only. This will create a new HashMap for
     * each invocation and add all values.
     */
    @Override
    public Set<Entry<Integer, V>> entrySet() {
        Map<Integer, V> map = new HashMap<>();
        for (int i = 0; i < values.length; i++) {
            map.put(i, values[i]);
        }
        return map.entrySet();
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
//        if(value == null) {
//            throw new IllegalArgumentException("Value must not be null");
//        }

        if(key >= this.values.length) {
            this.values = Arrays.copyOf(this.values, key + reservoirSize);
        }

        if(values[key] == null) {
            this.values[key] = value;
            this.isEmpty = false;
            this.size++;
        } else {
            this.values[key] = value;
        }
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

    private class IndexMapValuesCollection extends AbstractCollection<V> {

        /**
         * Returns an iterator over the elements contained in this collection.
         *
         * @return an iterator over the elements contained in this collection
         */
        @Override
        public Iterator<V> iterator() {
            return new IndexMapIterator();
        }

        @Override
        public int size() {
            return IndexMap.this.size();
        }
    }

    private class IndexMapIterator implements Iterator<V> {
        private int cursor = 0;

        /**
         * Returns {@code true} if the iteration has more elements.
         * (In other words, returns {@code true} if {@link #next} would
         * return an element rather than throwing an exception.)
         *
         * @return {@code true} if the iteration has more elements
         */
        @Override
        public boolean hasNext() {
            while(cursor < values.length) {
                if(values[cursor] != null) {
                    return true;
                }
                cursor++;
            }
            return false;
        }

        /**
         * Returns the next element in the iteration.
         *
         * @return the next element in the iteration
         * @throws NoSuchElementException if the iteration has no more elements
         */
        @Override
        public V next() {
            cursor++;
            while(cursor - 1 < values.length) {
                if(values[cursor - 1] != null) {
                    return values[cursor - 1];
                }
                cursor++;
            }
            return null;
        }
    }
}
