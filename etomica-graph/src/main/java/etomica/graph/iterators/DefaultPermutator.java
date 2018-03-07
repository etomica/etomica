/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.iterators;

import java.util.NoSuchElementException;

import etomica.graph.model.Permutator;

/*
 * Consider the partition P of a given set S of objects, where each partition
 * Pi consists of |Pi| distinct objects of the same type. The total number of
 * elements in S is the sum of |Pi| for all i.
 *
 * This iterator produces all distinct permutations of these |S| objects, where
 * objects in each partition Pi are treated as equivalents. There are exactly
 *
 *   choose(|S|,|P1|) x choose(|S|-|P1|,|P2|) x ... x choose(|Pn|,|Pn|)
 *
 * permutations in the resulting set of objects.
 *
 * The constructor takes as input int[] psizes. We assume that i = Type(Pi) and
 * that psizes[i] = |Pi|. This means that we map each partition Pi to an integer
 * i and that we assume that the size of the partitions are given by psizes. The
 * iterator produces a vector consisting of |Pi| entries with value i. The client
 * code can map these integer values to any types of objects needed. Note that the
 * sum of sizes of the partitions must fit in a byte.
 *
 * @author Demian Lessa
 */
public class DefaultPermutator implements Permutator {

  private int size = 0;
  private int[][] groups;
  private int activeGroup = -1;
  private int curPermutation = 0;
  private long maxPermutations = 1;
  private int[] groupSizes;

  public DefaultPermutator(int[] psizes) {

    /**
     * Copy the partition sizes array into a local array to protect the iterator
     * from external changes to this array.
     */
    groupSizes = new int[psizes.length];
    System.arraycopy(psizes, 0, groupSizes, 0, psizes.length);
    /**
     * Compute the total size of the permutations we will generate.
     */
    for (int pid = 0; pid < groupSizes.length; pid++) {
      size += groupSizes[pid];
    }
    /**
     * A group consists of a number of cells in which we allow the permutations
     * of the partition nodes to orbit.
     */
    int groupSize = size;
    groups = new int[groupSizes.length][];
    for (int pid = 0; pid < groupSizes.length; pid++) {
      groups[pid] = new int[groupSize];
      for (int i = 0; i < groupSizes[pid]; i++) {
        groups[pid][i] = 1;
      }
      maxPermutations *= choose(groups[pid].length, groupSizes[pid]);
      groupSize -= groupSizes[pid];
    }
    activeGroup = groupSizes.length - 1;
  }

  // check if the current permutation is within the range of
  // permutations we iterate through
  public boolean hasNext() {

    return (curPermutation < maxPermutations);
  }

  // pop the current permutation then push the next one
  public byte[] next() {

    if (!hasNext()) {
      throw new NoSuchElementException();
    }
    byte[] result = currentPermutation();
    nextPermutation();
    return result;
  }

  // this iterator is read-only
  public void remove() {

    throw new UnsupportedOperationException();
  }

  private long rangeProduct(int minVal, int maxVal) {

    assert ((minVal >= 0) && (maxVal >= minVal));
    long result = 1;
    for (int val = minVal + (minVal == 0 ? 1 : 0); val <= maxVal; val++) {
      result *= val;
    }
    return result;
  }

//  private long factorial(int val) {
//
//    assert (val >= 0);
//    if (val <= 1) {
//      return 1;
//    }
//    else {
//      return rangeProduct(1, val);
//    }
//  }

  private long choose(int n, int k) {

    if ((k == 0) || (n - k == 0)) {
      return 1;
    }
    else if ((k == 1) || (n - k == 1)) {
      return n;
    }
    if (k < n / 2) {
      return rangeProduct(n - k + 1, n) / rangeProduct(1, k);
    }
    else {
      return rangeProduct(k + 1, n) / rangeProduct(1, n - k);
    }
  }

  private void nextPermutation() {

    curPermutation++;
    if (curPermutation >= maxPermutations) {
      return;
    }
    // if necessary, backtrack until we find the next group to update
    while ((activeGroup != 0) && !hasGroupPermutation()) {
      // reset the active group
      for (int i = 0; i < groups[activeGroup].length; i++) {
        groups[activeGroup][i] = (i < groupSizes[activeGroup] ? 1 : 0);
      }
      activeGroup--;
    }
    nextGroupPermutation();
  }

  private int tailCount() {

    int tailCount = 0;
    for (int i = groups[activeGroup].length - 1; i >= 0; i--) {
      if (groups[activeGroup][i] == 0) {
        break;
      }
      tailCount++;
    }
    return tailCount;
  }

  private int headCount() {

    int headCount = 0;
    for (int i = 0; i < groups[activeGroup].length; i++) {
      if (groups[activeGroup][i] == 0) {
        break;
      }
      headCount++;
    }
    return headCount;
  }

  /**
   * Checks whether the current group still has permutations to contribute.
   */
  private boolean hasGroupPermutation() {

    // Check if all 1 bits in this group are at the tail of the group.
    return tailCount() != groupSizes[activeGroup];
  }

  // the current permutation is the merge of ids of the
  // partitions that contribute each bit
  private byte[] currentPermutation() {

    byte[] result = new byte[size];
    for (int i = 0; i < size; i++) {
      result[i] = groupID(i);
    }
    return result;
  }

  /**
   * Computes the id of the group that represents the given bit.
   */
  private byte groupID(int bitIndex) {

    int mappedIndex = bitIndex;
    for (byte gid = 0; gid < groups.length; gid++) {
      if (groups[gid][mappedIndex] == 1) {
        return gid;
      }
      else {
        int zeroCount = 0;
        for (int index = 0; index < mappedIndex; index++) {
          if (groups[gid][index] == 0) {
            zeroCount++;
          }
        }
        mappedIndex = zeroCount;
      }
    }
    // the permutation index must refer to some permutation;
    // if we reach here, there is a problem...
    assert (false);
    return -1;
  }

  /**
   * The active group has not been exhausted by the iterator. Move to the next
   * permutation in the active group. Then, update the active group to the last
   * group.
   */
  private void nextGroupPermutation() {

    int tailCount = tailCount();
    int headCount = headCount();
    /*
     * In this particular case, all 1 bits are at the head and the tail
     * of the group, so it's time to restart with one less bit at the head.
     */
    if ((tailCount > 0) && (headCount > 0)
        && (tailCount + headCount == groupSizes[activeGroup])) {
      for (int i = 1; i < headCount - 1; i++) {
        groups[activeGroup][i - 1] = 1;
      }
      groups[activeGroup][headCount - 1] = 0;
      for (int i = headCount; i < headCount + tailCount + 1; i++) {
        groups[activeGroup][i] = 1;
      }
      for (int i = headCount + tailCount + 1; i < groups[activeGroup].length; i++) {
        groups[activeGroup][i] = 0;
      }
      // make sure we're now tracking the correct group
      activeGroup = groups.length - 1;
      return;
    }
    /*
     * This is the general case in which there exists at least
     * one 1 bit and one 0 bit between the head and tail.
     */
    int zeroIndex = -1;
    int oneIndex = -1;
    // find the rightmost 1 to the left of a 0
    for (int i = groups[activeGroup].length - 1; i >= 0; i--) {
      if ((zeroIndex == -1) && (groups[activeGroup][i] == 0)) {
        zeroIndex = i;
      }
      else if ((zeroIndex != -1) && (groups[activeGroup][i] == 1)) {
        oneIndex = i;
        zeroIndex = -1;
        break;
      }
    }
    // find the first zero to the right of the 1
    for (int i = oneIndex; i < groups[activeGroup].length; i++) {
      if (groups[activeGroup][i] == 0) {
        zeroIndex = i;
        break;
      }
    }
    // we must have found the 1 bit and the 0 bit
    assert ((oneIndex != -1) && (zeroIndex != -1));
    groups[activeGroup][oneIndex] = 0;
    // clear the tail bits, if necessary
    for (int i = 0; i < tailCount; i++) {
      groups[activeGroup][groups[activeGroup].length - i - 1] = 0;
    }
    // set the correct bits
    for (int i = 0; i < tailCount + 1; i++) {
      groups[activeGroup][zeroIndex + i] = 1;
    }
    // make sure we're now tracking the correct group
    activeGroup = groups.length - 1;
  }
}
