package etomica.units;

/**
 * Class to form a dimension from ratio of two dimensions.
 * Used primarily to construct intensive quantities, such as energy/volume.
 */
public final class DimensionRatio extends Dimension {
    private final Dimension nDimension;
    private final Dimension dDimension;
    private final double[] signature = new double[4];
    public DimensionRatio(Dimension numerator, Dimension denominator) {
        nDimension = numerator;
        dDimension = denominator;
        for(int i=0; i<numerator.signature().length; i++) {
            signature[i] = numerator.signature()[i] - denominator.signature()[i];
        }
    }
    public double[] signature() {return signature;}
    public Dimension nDimension() {return nDimension;}
    public Dimension dDimension() {return dDimension;}
    public Unit defaultIOUnit() {return new UnitRatio(nDimension.defaultIOUnit(),dDimension.defaultIOUnit());}
    public Class baseUnit() {return UnitRatio.class;}
    public String toString() {return nDimension.toString() + "/" + dDimension.toString();}
}
    
