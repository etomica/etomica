package etomica.graph.test;

import java.util.Arrays;

import etomica.graph.iterators.CartesianPermutator;
import etomica.graph.model.BytePermutator;

import junit.framework.TestCase;

public class CartesianPermutatorTest extends TestCase {

  public void testPermutator1x1x2x1() {

    BytePermutator p = new CartesianPermutator(new int[] { 1, 1 }, new int[] { 2, 1 });
    int count = 0;
    while (p.hasNext()) {
      byte[] b = p.next();
      count++;
      System.out.println(String.format("%d :: %s", count, Arrays.toString(b)));
    }
    System.out.println();
  }

  public void testPermutator1x1x2x2() {

    BytePermutator p = new CartesianPermutator(new int[] { 1, 1 }, new int[] { 2, 2 });
    int count = 0;
    while (p.hasNext()) {
      byte[] b = p.next();
      count++;
      System.out.println(String.format("%d :: %s", count, Arrays.toString(b)));
    }
    System.out.println();
  }

  public void testPermutator1x1x3x3() {

    BytePermutator p = new CartesianPermutator(new int[] { 1, 1 }, new int[] { 3, 3 });
    int count = 0;
    while (p.hasNext()) {
      byte[] b = p.next();
      count++;
      System.out.println(String.format("%d :: %s", count, Arrays.toString(b)));
    }
    System.out.println();
  }
}