/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph;

import java.util.Arrays;

import etomica.graph.iterators.CartesianPermutator;
import etomica.graph.model.Permutator;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Disabled;

@Disabled
public class CartesianPermutatorTest {

  @Test
  public void testPermutator1x1x2x1() {

    Permutator p = new CartesianPermutator(new int[] { 1, 1 }, new int[] { 2, 1 });
    int count = 0;
    while (p.hasNext()) {
      byte[] b = p.next();
      count++;
      System.out.println(String.format("%d :: %s", count, Arrays.toString(b)));
    }
    System.out.println();
  }

  @Test
  public void testPermutator1x1x2x2() {

    Permutator p = new CartesianPermutator(new int[] { 1, 1 }, new int[] { 2, 2 });
    int count = 0;
    while (p.hasNext()) {
      byte[] b = p.next();
      count++;
      System.out.println(String.format("%d :: %s", count, Arrays.toString(b)));
    }
    System.out.println();
  }

  @Test
  public void testPermutator1x1x3x3() {

    Permutator p = new CartesianPermutator(new int[] { 1, 1 }, new int[] { 3, 3 });
    int count = 0;
    while (p.hasNext()) {
      byte[] b = p.next();
      count++;
      System.out.println(String.format("%d :: %s", count, Arrays.toString(b)));
    }
    System.out.println();
  }
}
