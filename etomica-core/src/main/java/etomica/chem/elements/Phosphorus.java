package etomica.chem.elements;

public class Phosphorus extends ElementChemical{
    protected Phosphorus(String symbol){this(symbol,30.974);}

    public Phosphorus(String symbol, double mass) {
        super(symbol, mass, 12);
    }

    private static final long serialVersionUID = 1L;
    public static final Phosphorus INSTANCE = new Phosphorus("P");
}
