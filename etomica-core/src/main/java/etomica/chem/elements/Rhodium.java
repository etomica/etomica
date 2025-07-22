package etomica.chem.elements;

public class Rhodium extends ElementChemical {

    protected Rhodium(String symbol) {
        this(symbol, 102.9);
    }

    protected Rhodium(String symbol, double mass) {
        super(symbol, mass, 45);
    }

    public static final Rhodium INSTANCE = new Rhodium("Rh");

    private static final long serialVersionUID = 1L;
}
