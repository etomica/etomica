package etomica.units;

public class CompoundDimension extends Dimension {

    public CompoundDimension(Dimension[] dimensions, double[] exponents) {
        this(dimensions, exponents, makeName(dimensions, exponents));
    }
    
    public CompoundDimension(Dimension[] dimensions, double[] exponents, String name) {
        super(name, makeSignature(dimensions, exponents));
        this.dimensions = (Dimension[])dimensions.clone();
        this.exponents = (double[])exponents.clone();
    }
    
    // used by constructor
    private static double[] makeSignature(Dimension[] dimensions, double[] exponents) {
        if(dimensions.length != exponents.length) {
            throw new IllegalArgumentException("Arguments to CompoundDimension constructor must be arrays of the same length.  Given were arrays of length "+dimensions.length+" and "+exponents.length);
        }
        double[] sig = new double[Dimension.N_BASE]; 
        for(int i=0; i<dimensions.length; i++) {
            for(int j=0; i<sig.length; j++) {
                sig[j] += dimensions[i].signature()[j];
            }
        }
        return sig;
    }
    // used by constructor
    private static String makeName(Dimension[] dimensions, double[] exponents) {
        if(dimensions.length != exponents.length) {
            throw new IllegalArgumentException("Arguments to CompoundDimension constructor must be arrays of the same length.  Given were arrays of length "+dimensions.length+" and "+exponents.length);
        }
        String symbol = "";
        for(int i=0; i<dimensions.length; i++) {
            symbol += "("+dimensions[i].toString()+"^"+exponents[i]+")";
        }
        return symbol;
    }
    
    public Unit getUnit(UnitSystem unitSystem) {
        Unit[] units = new Unit[dimensions.length];
        for(int i=0; i<units.length; i++) {
            units[i] = dimensions[i].getUnit(unitSystem);
        }
        return new CompoundUnit(units, exponents);
    }
    
    private final Dimension[] dimensions;
    private final double[] exponents;
    private static final long serialVersionUID = 1;
}
