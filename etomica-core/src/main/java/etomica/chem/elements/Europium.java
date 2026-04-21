package etomica.chem.elements;

public class Europium extends ElementChemical {
    protected Europium(String symbol) {
        this(symbol, 151.96);
    }
    protected Europium(String symbol, double mass) {
        super(symbol, mass, 63);
    }

    public static final Europium INSTANCE = new Europium("Eu");

    private static final long serialVersionUID = 1L;
}
