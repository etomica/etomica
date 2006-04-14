package etomica.virial;


public interface CoordinatePairSet {

    /**
     * Returns atom pair for ith and jth atoms in set.
     */
    public double getr2(int i, int j);

    public void reset();

    public int getID();

}