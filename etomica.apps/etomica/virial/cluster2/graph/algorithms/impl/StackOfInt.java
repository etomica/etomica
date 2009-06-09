package etomica.virial.cluster2.graph.algorithms.impl;

public class StackOfInt {

  private int[] stack;
  private int topIndex = -1;

  public StackOfInt(int capacity) {

    stack = new int[capacity];
  }

  public boolean isEmpty() {

    return topIndex == -1;
  }

  public int pop() {

    return stack[topIndex--];
  }

  public void push(int value) {

    stack[++topIndex] = value;
  }

  public int size() {

    return topIndex + 1;
  }

  public int top() {

    return stack[topIndex];
  }

  public int[] toIntArray() {
    
    int[] result = new int[size()];
    System.arraycopy(stack, 0, result, 0, size());
    return result;
  }
}