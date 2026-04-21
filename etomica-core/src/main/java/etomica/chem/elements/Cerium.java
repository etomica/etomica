package etomica.chem.elements;

public class Cerium extends ElementChemical {

    protected Cerium(String symbol) {
        this(symbol, 140.12);
    }

    protected Cerium(String symbol, double mass) {
        super(symbol, mass, 58);
    }

    public static final Cerium INSTANCE = new Cerium("Ce");

    private static final long serialVersionUID = 1L;
}

