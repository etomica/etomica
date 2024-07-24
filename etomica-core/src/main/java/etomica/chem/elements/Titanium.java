package etomica.chem.elements;

public class Titanium extends ElementChemical {

    protected Titanium(String symbol) {
        this(symbol, 47.867);
    }

    protected Titanium(String symbol, double mass) {
        super(symbol, mass, 22);
    }

    public static final Titanium INSTANCE = new Titanium("Ti");

    private static final long serialVersionUID = 1L;
}
