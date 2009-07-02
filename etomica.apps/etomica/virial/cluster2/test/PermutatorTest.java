package etomica.virial.cluster2.test;

import java.util.Iterator;

import etomica.virial.cluster2.util.PermutationIterator;

public class PermutatorTest extends CustomTestCase {

  private void println(int[] vector, int base) {

    String result = "";
    for (int i = 0; i < vector.length; i++) {
      result += base - vector[i] - 1;
    }
    System.out.println(result);
  }

  public void testPermutator1() {

    Iterator<int[]> iter = new PermutationIterator(new int[] { 6, 3 });
    while (iter.hasNext()) {
      println(iter.next(), 2);
    }
  }

  public void testPermutator2() {

    Iterator<int[]> iter = new PermutationIterator(new int[] { 4, 3, 3 });
    while (iter.hasNext()) {
      println(iter.next(), 3);
    }
  }
}