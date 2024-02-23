package etomica.chem.elements;

public class Carbon_2 extends ElementChemical{
    protected Carbon_2(String symbol) {
        this(symbol, 12.0107);
    }

    protected Carbon_2(String symbol, double mass) {
        super(symbol, mass, 6);
    }

    public static final Carbon_2 INSTANCE = new Carbon_2("C_2");

    private static final long serialVersionUID = 1L;
}
