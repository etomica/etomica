package etomica.chem.elements;

public class Dysprosium extends ElementChemical {

    protected Dysprosium(String symbol) {
        this(symbol, 162.5);
    }

    protected Dysprosium(String symbol, double mass) {
        super(symbol, mass, 66);
    }

    public static final Dysprosium INSTANCE = new Dysprosium("Dy");

    private static final long serialVersionUID = 1L;
}
