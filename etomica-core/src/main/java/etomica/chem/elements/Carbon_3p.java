package etomica.chem.elements;

public class Carbon_3p extends ElementChemical{
    protected Carbon_3p(String symbol) {
        this(symbol, 12.0107);
    }

    protected Carbon_3p(String symbol, double mass) {
        super(symbol, mass, 6);
    }

    public static final Carbon_3p INSTANCE = new Carbon_3p("C_3p");

    private static final long serialVersionUID = 1L;
}
