package etomica.chem.elements;

public class Erbium extends ElementChemical {

    protected Erbium(String symbol) {
        this(symbol, 167.26);
    }

    protected Erbium(String symbol, double mass) {
        super(symbol, mass, 68);
    }

    public static final Erbium INSTANCE = new Erbium("Er");

    private static final long serialVersionUID = 1L;
}

