/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.iterators;

import etomica.graph.model.Permutator;

public class CartesianPermutator implements Permutator {

  private int[] innerPartition;
  private int[] outerPartition;
  private byte[] outerPermutation = null;
  private DefaultPermutator outerPermutator = null;
  private DefaultPermutator innerPermutator = null;

  public CartesianPermutator(int[] outerPartition, int[] innerPartition) {

    this.innerPartition = innerPartition;
    this.outerPartition = outerPartition;
    bootstrap();
  }

  protected void bootstrap() {

    setOuterPermutator(createOuterPermutator());
    if (getOuterPermutator().hasNext()) {
      setOuterPermutation(getOuterPermutator().next());
    }
    setInnerPermutator(createInnerPermutator());
  }

  // we just computed a permutation in the cartesian; it is time to advance to the next
  // positions in the inner and outer permutators
  protected void advance() {

    // if the inner iteration is over, we must advance the outer iteration
    if (!getInnerPermutator().hasNext()) {
      setOuterPermutation(null);
      // if there is no next outer graph, then we are done;
      // otherwise, we get a new outer graph and a new inner iterator
      if (getOuterPermutator().hasNext()) {
        setOuterPermutation(getOuterPermutator().next());
        setInnerPermutator(createInnerPermutator());
      }
    }
  }

  protected byte[] combinePermutations(byte[] outer, byte[] inner) {

    byte[] result = new byte[outer.length + inner.length];
    System.arraycopy(outer, 0, result, 0, outer.length);
    System.arraycopy(inner, 0, result, outer.length, inner.length);
    return result;
  }

  public DefaultPermutator createInnerPermutator() {

    return new DefaultPermutator(innerPartition);
  }

  public DefaultPermutator createOuterPermutator() {

    return new DefaultPermutator(outerPartition);
  }

  public DefaultPermutator getInnerPermutator() {

    return innerPermutator;
  }

  public byte[] getOuterPermutation() {

    return outerPermutation;
  }

  public DefaultPermutator getOuterPermutator() {

    return outerPermutator;
  }

  // The cartesian product can continue for as long as there exists an outer graph
  // and an inner graph can be obtained from the inner iterator
  public boolean hasNext() {

    return getOuterPermutation() != null && getInnerPermutator().hasNext();
  }

  public byte[] next() {

    if (hasNext()) {
      byte[] result = combinePermutations(getOuterPermutation(), getInnerPermutator().next());
      advance();
      return result;
    }
    return null;
  }

  public void remove() {

    // no-op
  }

  protected void setInnerPermutator(DefaultPermutator permutator) {

    innerPermutator = permutator;
  }

  protected void setOuterPermutation(byte[] permutation) {

    outerPermutation = permutation;
  }

  protected void setOuterPermutator(DefaultPermutator permutator) {

    outerPermutator = permutator;
  }
}