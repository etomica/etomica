package etomica.units;

public class CompoundUnit implements Unit {

    public CompoundUnit(Unit[] units, double[] exponents) {
        this(units, exponents, makeName(units, exponents), makeSymbol(units, exponents));
    }
    
    public CompoundUnit(Unit[] units, double[] exponents, String name, String symbol) {
        if(units.length != exponents.length) {
            throw new IllegalArgumentException("Arguments to CompoundUnit constructor must be arrays of the same length.  Given were arrays of length "+units.length+" and "+exponents.length);
        }
        this.symbol = symbol;
        double f = 1.0;
        Dimension[] dimensions = new Dimension[units.length];
        for(int i=0; i<units.length; i++) {
            f *= Math.pow(units[i].fromSim(1.0),exponents[i]);
            dimensions[i] = units[i].dimension();
        }
        from = f;
        dimension = new CompoundDimension(dimensions, exponents);
    }
    
    // used by constructor
    private static String makeName(Unit[] units, double[] exponents) {
        if(units.length != exponents.length) {
            throw new IllegalArgumentException("Arguments to CompoundUnit constructor must be arrays of the same length.  Given were arrays of length "+units.length+" and "+exponents.length);
        }
        String name = "";
        for(int i=0; i<units.length; i++) {
            name += "("+units[i].toString();
            if (exponents[i] != 1.0) {
                name += "^"+exponents[i];
            }
            name += ")";
        }
        return name;
    }
    
    // used by constructor
    private static String makeSymbol(Unit[] units, double[] exponents) {
        if(units.length != exponents.length) {
            throw new IllegalArgumentException("Arguments to CompoundUnit constructor must be arrays of the same length.  Given were arrays of length "+units.length+" and "+exponents.length);
        }
        String symbol = "";
        for(int i=0; i<units.length; i++) {
            symbol += units[i].symbol();
            if (exponents[i] != 1.0) {
                symbol += "^"+exponents[i];
            }
            if(i != units.length-1) symbol += "-";
        }
        return symbol;
    }
    
    public Dimension dimension() {
        return dimension;
    }

    public double toSim(double x) {
        return x / from;
    }

    public double fromSim(double x) {
        return x * from;
    }

    public String symbol() {
        return symbol;
    }

    public boolean prefixAllowed() {
        return false;
    }

    private final double from;
    private final String symbol;
    private final Dimension dimension;
    private static final long serialVersionUID = 1;
}
