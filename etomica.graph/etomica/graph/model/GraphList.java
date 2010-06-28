package etomica.graph.model;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Set;

import etomica.graph.model.comparators.ComparatorChain;

public class GraphList<E> implements Set<E> {

    protected final ArrayList<E> list;
    protected final Comparator<E> comparator;
    protected boolean sorted = true;

    public GraphList() {
        this((Comparator<E>)new ComparatorChain());
    }
    
    public GraphList(Comparator<E> comparator) {
        list = new ArrayList<E>();
        this.comparator = comparator;
    }

    public boolean add(E e) {
        list.add(e);
        sorted = false;
        return true;
    }

    public boolean addAll(Collection<? extends E> c) {
        list.addAll(c);
        sorted = false;
        return true;
    }
    
    protected void sort() {
        if (sorted) return;
        if (comparator != null) Collections.sort(list, comparator);
        sorted = true;
    }

    public void clear() {
        list.clear();
        sorted = true;
    }

    public boolean contains(Object o) {
        return list.contains(o);
    }

    public boolean containsAll(Collection<?> c) {
        return list.containsAll(c);
    }

    public boolean isEmpty() {
        return list.isEmpty();
    }

    public Iterator<E> iterator() {
        sort();
        return list.iterator();
    }

    public boolean remove(Object o) {
        return list.remove(o);
    }

    public boolean removeAll(Collection<?> c) {
        return list.removeAll(c);
    }

    public boolean retainAll(Collection<?> c) {
        return list.retainAll(c);
    }

    public int size() {
        return list.size();
    }

    public Object[] toArray() {
        sort();
        return list.toArray();
    }

    public <T> T[] toArray(T[] a) {
        sort();
        return list.toArray(a);
    }
}
