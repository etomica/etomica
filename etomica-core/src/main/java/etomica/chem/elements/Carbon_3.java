package etomica.chem.elements;

public class Carbon_3 extends ElementChemical{
    protected Carbon_3(String symbol) {
        this(symbol, 12.0107);
    }

    protected Carbon_3(String symbol, double mass) {
        super(symbol, mass, 6);
    }

    public static final Carbon_3 INSTANCE = new Carbon_3("C_3");

    private static final long serialVersionUID = 1L;
}
