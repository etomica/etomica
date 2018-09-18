/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.model;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Set;

import etomica.graph.model.comparators.ComparatorChain;

public class GraphList implements Set<Graph> {

    protected final ArrayList<Graph> list;
    protected final Comparator<Graph> comparator;
    protected boolean sorted = true;

    public GraphList() {
        this(new ComparatorChain());
    }
    
    public GraphList(Comparator<Graph> comparator) {
        list = new ArrayList<>();
        this.comparator = comparator;
    }

    public boolean add(Graph e) {
        list.add(e);
        sorted = false;
        return true;
    }

    public boolean addAll(Collection<? extends Graph> c) {
        list.addAll(c);
        sorted = false;
        return true;
    }
    
    protected void sort() {
        if (sorted) return;
        if (comparator != null) list.sort(comparator);
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

    public Iterator<Graph> iterator() {
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

    public String toString() {
      if (list.size() == 0) return "";
      String out = list.get(0).toString();
      for (int i=1; i<list.size(); i++) {
        out += " "+list.get(i);
      }
      return out;
    }
}
