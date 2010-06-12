package etomica.graph.model;

import java.util.Collection;
import java.util.Iterator;
import java.util.Set;

import etomica.graph.iterators.ChainedIterator;

public class GraphSuperSet<E> implements Set<E> {

  protected final SubSetJudge judge;
  protected Set<Graph>[] subSets;
  
  public GraphSuperSet(SubSetJudge judge, Set<Graph>[] subSets) {
    this.judge = judge;
    this.subSets = subSets;
  }
  
  public boolean add(E e) {
    subSets[judge.getSetId((Graph)e)].add((Graph)e);
    return true;
  }

  public boolean addAll(Collection<? extends E> c) {
    for (E e : c) {
      add(e);
    }
    return true;
  }

  public void clear() {
    for (int i=0; i<subSets.length; i++) {
      subSets[i].clear();
    }
  }

  public boolean contains(Object o) {
    for (int i=0; i<subSets.length; i++) {
      if (subSets[i].contains(o)) return true;
    }
    return false;
  }

  public boolean containsAll(Collection<?> c) {
    return false;
  }

  public boolean isEmpty() {
    for (int i=0; i<subSets.length; i++) {
      if (!subSets[i].isEmpty()) return false;
    }
    return true;
  }
  
  public Iterator<E> iterator() {
    ChainedIterator iter = new ChainedIterator();
    for (int i=0; i<subSets.length; i++) {
      if (!subSets[i].isEmpty()) {
        iter.chainIterator(subSets[i].iterator());
      }
    }
    iter.start();
    return (Iterator<E>)iter;
  }

  public boolean remove(Object o) {
    throw new RuntimeException("nope");
  }

  public boolean removeAll(Collection<?> c) {
    throw new RuntimeException("nope");
  }

  public boolean retainAll(Collection<?> c) {
    throw new RuntimeException("nope");
  }

  @Override
  public int size() {
    int s = 0;
    for (int i=0; i<subSets.length; i++) {
      s += subSets[i].size();
    }
    return s;
  }
  
  public String toString() {
    if (isEmpty()) return "[]";
    String str = "";
    for (int i=0; i<subSets.length; i++) {
      if (!subSets[i].isEmpty()) {
        str += subSets[i].toString();
      }
    }
    return str;
  }

  public Object[] toArray() {
    return null;
  }

  public <T> T[] toArray(T[] a) {
    return null;
  }

  public interface SubSetJudge {
    public int getSetId(Graph g);
  }
}
