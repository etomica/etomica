package etomica.chem.elements;

public class Thulium extends ElementChemical {

    protected Thulium(String symbol) {
        this(symbol, 168.93);
    }

    protected Thulium(String symbol, double mass) {
        super(symbol, mass, 69);
    }

    public static final Thulium INSTANCE = new Thulium("Tm");

    private static final long serialVersionUID = 1L;
}
