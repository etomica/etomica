package etomica.graph.test;

import java.util.Set;

import junit.framework.TestCase;

import etomica.graph.iterators.DefaultIterator;
import etomica.graph.iterators.IteratorToSet;
import etomica.graph.iterators.StoredIterator;
import etomica.graph.iterators.filters.IsomorphismFilter;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;
import etomica.graph.operations.Binary;
import etomica.graph.operations.Delete;


public class DeleteTest extends TestCase {

  public void testDelete1() {

    byte rangeBegin = 2;
    byte rangeEnd = 6;
    Binary delete = new Delete();
    IteratorToSet its = new IteratorToSet();
    for (byte nodeCount = rangeBegin; nodeCount <= rangeEnd; nodeCount++) {
      GraphIterator i1 = new IsomorphismFilter(new DefaultIterator(nodeCount));
      GraphIterator i2 = new StoredIterator(nodeCount);
      Set<Graph> s1 = its.getSet(i1);
      Set<Graph> s2 = its.getSet(i2);
      Set<Graph> r = delete.apply(s1, s2, null);
      assertTrue(r.size() == 0);
      show(s1, s2, r);
      r = delete.apply(s2, s1, null);
      assertTrue(r.size() == 0);
      show(s2, s1, r);
      r = delete.apply(s1, s1, null);
      assertTrue(r.size() == 0);
      show(s1, s1, r);
      r = delete.apply(s2, s2, null);
      assertTrue(r.size() == 0);
      show(s2, s2, r);
    }
  }

  public void testDelete2() {

    byte rangeBegin = 2;
    byte rangeEnd = 6;
    Binary delete = new Delete();
    IteratorToSet its = new IteratorToSet();
    for (byte nodeCount = rangeBegin; nodeCount <= rangeEnd; nodeCount++) {
      GraphIterator i1 = new DefaultIterator(nodeCount);
      GraphIterator i2 = new StoredIterator(nodeCount);
      Set<Graph> s1 = its.getSet(i1);
      Set<Graph> s2 = its.getSet(i2);
      Set<Graph> r = delete.apply(s2, s1, null);
      assertTrue(r.size() == 0);
      show(s2, s1, r);
      r = delete.apply(s1, s1, null);
      assertTrue(r.size() == 0);
      show(s1, s1, r);
      r = delete.apply(s2, s2, null);
      assertTrue(r.size() == 0);
      show(s2, s2, r);
    }
  }

  private void show(Set<Graph> s1, Set<Graph> s2, Set<Graph> result) {

    for (Graph g : result) {
      System.out.println(g);
    }
    System.out.println();
    System.out.println("s1 total.....: " + s1.size() + " graphs");
    System.out.println("s2 total.....: " + s2.size() + " graphs");
    System.out.println("del total....: " + result.size() + " graphs\n");
  }
}