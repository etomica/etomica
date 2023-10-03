package etomica.chem.elements;

public class Carbon_1p extends ElementChemical{
    protected Carbon_1p(String symbol) {
        this(symbol, 12.0107);
    }

    protected Carbon_1p(String symbol, double mass) {
        super(symbol, mass, 6);
    }

    public static final Carbon_1p INSTANCE = new Carbon_1p("C_1p");

    private static final long serialVersionUID = 1L;
}
