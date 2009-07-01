package etomica.virial.cluster2.util;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Consider the partition P = P1, ..., Pn of a set of objects S. Each partition
 * has an identifier (id) and a size. The total number of elements in the set of
 * objects is the sum of Pi.size for all Pi in P. This iterator produces all
 * distinct permutations of |S| objects such that a permutation consists of
 * Pi.size objects of type Pi.id for every i in [1..n].
 * 
 * @author Demian Lessa
 */
public class PermutationIterator implements Iterator<int[]> {

  private int size = 0;
  private int[][] groups;
  private int activeGroup = -1;
  private int curPermutation = 1;
  private int maxPermutations = 0;
  private int[] partitionSizes;

  public PermutationIterator(int[] psizes) {

    assert (partitionSizes != null);
    /**
     * Copy the partition sizes array into a local array to protect the iterator
     * from external changes to this array.
     */
    partitionSizes = new int[psizes.length];
    System.arraycopy(psizes, 0, partitionSizes, 0, psizes.length);
    /**
     * Compute the total size of the permutations we will generate.
     */
    for (int pid = 0; pid < partitionSizes.length; pid++) {
      size += partitionSizes[pid];
    }
    /**
     * A group consists of a number of cells in which we allow the permutations
     * of the partition nodes to orbit.
     */
    int groupSize = size;
    groups = new int[partitionSizes.length][];
    for (int pid = 0; pid < partitionSizes.length; pid++) {
      groups[pid] = new int[groupSize];
      for (int i = 0; i < partitionSizes[pid]; i++) {
        groups[pid][i] = 1;
      }
      maxPermutations += choose(groupSize, partitionSizes[pid]);
      groupSize -= partitionSizes[pid];
    }
    activeGroup = partitionSizes.length - 1;
  }

  private int factorial(int val) {

    assert (val >= 0);
    if (val > 1) {
      return val * factorial(val - 1);
    }
    return 1;
  }

  private int choose(int n, int k) {

    return factorial(n) / (factorial(n - k) * factorial(k));
  }

  private void nextPermutation() {

    curPermutation++;
    if (curPermutation > maxPermutations) {
      return;
    }
    while (!hasGroupPermutation()) {
      nextGroup();
    }
    nextGroupPermutation();
  }

  /**
   * The active group has not been exhausted by the iterator. Move to the next
   * permutation in the active group. Update the active group, if necessary, to
   * the last group.
   */
  private void nextGroupPermutation() {

    // TODO Auto-generated method stub
  }

  /**
   * Some group has been exhausted by the iterator. Backtrack and leave the
   * groups in a state where the next permutation is available.
   */
  private void nextGroup() {

    // TODO Auto-generated method stub
  }

  /**
   * Checks whether the current group still has permutations to contribute.
   */
  private boolean hasGroupPermutation() {

    /**
     * All we do is check if all 1 bits in this group's array are at the tail of
     * the group.
     */
    int tailCount = 0;
    for (int i = groups[activeGroup].length; i >= 0; i--) {
      if (groups[activeGroup][i] == 0) {
        break;
      }
      tailCount++;
    }
    return tailCount != partitionSizes[activeGroup];
  }

  private int[] currentPermutation() {

    int[] result = new int[size];
    for (int i = 0; i < size; i++) {
      result[i] = partitionID(i);
    }
    return result;
  }

  /**
   * Computes the partition id of the given permutation index.
   */
  private int partitionID(int permIndex) {

    int mappedIndex = permIndex;
    for (int pid = 0; pid < size; pid++) {
      if (groups[pid][mappedIndex] == 1) {
        return pid;
      }
      else {
        int zeroCount = 0;
        for (int index = 0; index < mappedIndex; index++) {
          if (groups[pid][index] == 0) {
            zeroCount++;
          }
        }
        mappedIndex = zeroCount;
      }
    }
    assert (false);
    return -1;
  }

  public boolean hasNext() {

    return (curPermutation <= maxPermutations);
  }

  public int[] next() {

    if (!hasNext()) {
      throw new NoSuchElementException();
    }
    int[] result = currentPermutation();
    nextPermutation();
    return result;
  }

  public void remove() {

    throw new UnsupportedOperationException();
  }
}