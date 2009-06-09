package etomica.virial.cluster2.graph.algorithms.impl;

public class QueueOfInt {

  private int[] queue;
  private int bottomIndex = -1;
  private int topIndex = -1;
  private int size = 0;

  public QueueOfInt(int capacity) {

    queue = new int[capacity];
  }

  public boolean isEmpty() {

    return size == 0;
  }

  public int bottom() {

    return queue[bottomIndex];
  }

  public int dequeue() {

    int result = queue[bottomIndex++];
    if (bottomIndex == queue.length) {
      bottomIndex = 0;
    }
    size--;
    return result;
  }

  public void enqueue(int value) {

    topIndex++;
    if (topIndex == queue.length) {
      topIndex = 0;
    }
    queue[topIndex] = value;
    size++;
  }

  public int size() {

    return size;
  }
}