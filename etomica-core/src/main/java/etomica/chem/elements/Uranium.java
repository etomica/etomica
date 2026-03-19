package etomica.chem.elements;

public class Uranium extends ElementChemical {

    protected Uranium(String symbol) {
        this(symbol, 238.02);
    }

    protected Uranium(String symbol, double mass) {
        super(symbol, mass, 92);
    }

    public static final Uranium INSTANCE = new Uranium("U");

    private static final long serialVersionUID = 1L;
}
