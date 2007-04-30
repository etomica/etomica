package etomica.math.discrete;

import etomica.lattice.IndexIterator;
import etomica.math.SpecialFunctions;

/**
 * Iterator that returns different combinations of n integers taken k at a time.
 * For example, if constructed with n = 4 and k = 2, will return, in successive
 * calls to next(), the arrays {0,1},{0,2},{0,3},{1,2},{1,3},{2,3}. In this example
 * each array has k = 2 elements, selected from the n = 4 possible values {0,1,2,3}
 * with no duplication of elements in any individual iterate. 
 * Entries in each array will be in numerical order.
 * 
 * @author kofke
 * 
 */

public class CombinationIterator implements IndexIterator, java.io.Serializable {

    /**
     * Constructor for CombinationIterator. Iterator requires reset before use.
     * Iterator parameters are final and cannot be changed after construction.
     * 
     * @param n
     *            the number of integers used to generate combinations
     * @param k
     *            the number of integers in each combination, and thus the length of the returned array
     * 
     * @throws IllegalArgumentException
     *             if values do not satisfy 0 <= k <= n.
     */
    public CombinationIterator(int n, int k) {
        super();
        if (n < 0 || k > n) {
            throw new IllegalArgumentException("CombinationIterator must have 0 <= k <= n");
        }
        this.n = n;
        this.k = k;
        index = new int[k];
        indexCopy = new int[k];
    }

    public void reset() {
        for (int i = 0; i < k; i++) {
            index[i] = i;
        }
        if (k != 0) {
            index[k - 1]--;
        } else {
            k0unset = false;
        }
    }

    public boolean hasNext() {
        if (k == 0) {
            return !k0unset;
        }
        for (int j = 0; j < k; j++) {
            if (index[j] < n - k + j) {
                return true;
            }
        }
        return false;
    }

    public int[] next() {
        if (k == 0) {
            k0unset = true;
            return index;
        }
        int j = k - 1;
        while (index[j] == n - k + j) {
            j--;
        }
        index[j]++;
        for (int i = j + 1; i < k; i++) {
            index[i] = index[i - 1] + 1;
        }
        System.arraycopy(index, 0, indexCopy, 0, k);
        return index;
    }
    
    public int getD() {
        return n;
    }

    public static void main(String[] args) {
        // int n = Integer.parseInt(args[0]);
        // int k = Integer.parseInt(args[1]);
        int n = 5;
        int k = 3;
        CombinationIterator p = new CombinationIterator(n, k);
        p.reset();
        int count = 1;
        while (p.hasNext()) {
            int[] a = p.next();
            System.out.print(count++ + ". {");
            for (int i = 0; i < a.length - 1; i++)
                System.out.print(a[i] + ", ");
            if (a.length > 0)
                System.out.print(a[a.length - 1]);
            System.out.println("}");
        }
        int s = SpecialFunctions.factorial(n) / SpecialFunctions.factorial(k) / SpecialFunctions.factorial(n - k);
        System.out.println("Expected total:" + s);
    }

    private final int[] index;
    private final int[] indexCopy;
    private final int n, k;
    private boolean k0unset = true;
    private static final long serialVersionUID = 1L;

}
