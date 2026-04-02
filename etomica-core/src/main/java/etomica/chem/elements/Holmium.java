package etomica.chem.elements;

public class Holmium extends ElementChemical {

    protected Holmium(String symbol) {
        this(symbol, 164.93);
    }

    protected Holmium(String symbol, double mass) {
        super(symbol, mass, 67);
    }

    public static final Holmium INSTANCE = new Holmium("Ho");

    private static final long serialVersionUID = 1L;
}
