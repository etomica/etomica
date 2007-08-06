package etomica.math.discrete;

/**
 * Generates int[] arrays such that the sum of (k+1) * array[k] over the elements
 * of array are equal to a given value. For example, for a target sum of 3 (as
 * given via setSum), this iterator generates {0,0,1}, {1,1,0}, and {3,0,0}
 * before expiring. For each case in this example, the indicated sum
 * [respectively (0*1 + 0*2 + 1*3), (1*1 + 1*2 + 0*3), (3*1 + 0*2 + 0*3)] the
 * indicated sum is equal to 3.  Number of elements in return array is set at construction.
 * 
 * @author David Kofke
 * 
 */
public class FixedSumIterator {

    public FixedSumIterator(int arrayLength) {
        array = new int[arrayLength];
        if (arrayLength != 1) {
            subIterator = new FixedSumIterator(arrayLength - 1);
        } else {
            subIterator = null;
        }
        index = -1;
    }

    /**
     * Sets the iterator to begin iteration.
     */
    public void reset() {
        index = targetSum / array.length;
        array[array.length - 1] = index;
        if (subIterator != null) {
            subIterator.setSum(targetSum - index * array.length);
            subIterator.reset();
        }
    }

    /**
     * Specifies the target sum that defines the iterates.  Iterator is
     * left in an unset state, so reset is required before invoking next().
     */
    public void setSum(int sum) {
        targetSum = sum;
        index = -1;
    }

    /**
     * Returns the next iterate, or null if iterator has expired.
     */
    public int[] next() {
        if (index < 0) {
            return null;
        }

        if (array.length == 1) {
            index = -1;
            return array;
        }

        int[] subNext = subIterator.next();
        if (subNext != null) {
            System.arraycopy(subNext, 0, array, 0, subNext.length);
            return array;
        }
        index--;
        array[array.length - 1] = index;
        if (index < 0) {
            return null;
        }
        subIterator.setSum(targetSum - index * array.length);
        subIterator.reset();
        return next();
    }

    public static void main(String[] args) {
        FixedSumIterator iterator = new FixedSumIterator(5);
        iterator.setSum(5);
        iterator.reset();
        for (int[] iterate = iterator.next(); iterate != null; iterate = iterator.next()) {
            System.out.println(etomica.util.Arrays.toString(iterate));
        }
    }

    private int targetSum;
    private int[] array;
    private final FixedSumIterator subIterator;
    private int index;
}
