package etomica.veos;

public interface BetaSource {

    /**
     * Returns k such that beta[j] == 0 for all j > k. Largest index for a non-zero beta. 
     */
    public int maxIndex();

    public double beta(int k);

}
