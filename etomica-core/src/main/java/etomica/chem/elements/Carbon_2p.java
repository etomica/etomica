package etomica.chem.elements;

public class Carbon_2p extends ElementChemical{
    protected Carbon_2p(String symbol) {
        this(symbol, 12.0107);
    }

    protected Carbon_2p(String symbol, double mass) {
        super(symbol, mass, 6);
    }

    public static final Carbon_2p INSTANCE = new Carbon_2p("C_2p");

    private static final long serialVersionUID = 1L;
}
