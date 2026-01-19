package etomica.chem.elements;

public class Ytterbium extends ElementChemical {

    protected Ytterbium(String symbol) {
        this(symbol, 173.04);
    }

    protected Ytterbium(String symbol, double mass) {
        super(symbol, mass, 70);
    }

    public static final Ytterbium INSTANCE = new Ytterbium("Yb");

    private static final long serialVersionUID = 1L;
}
