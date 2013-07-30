package etomica.virial;


public interface CoordinatePairSet {

    /**
     * Returns atom pair for ith and jth atoms in set.
     */
    public double getr2(int i, int j);

    /**
     * Informs the CoordinatePairSet that the configuration has changed and that it
     * has a new ID
     */
    public void reset(long cPairID);

    public long getID();

}