package etomica.chem.elements;

public class Carbon_1 extends ElementChemical{
    protected Carbon_1(String symbol) {
        this(symbol, 12.0107);
    }

    protected Carbon_1(String symbol, double mass) {
        super(symbol, mass, 6);
    }

    public static final Carbon_1 INSTANCE = new Carbon_1("C_1");

    private static final long serialVersionUID = 1L;
}
