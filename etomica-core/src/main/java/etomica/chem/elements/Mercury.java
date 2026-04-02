package etomica.chem.elements;

public class Mercury extends ElementChemical {

    protected Mercury(String symbol) {
        this(symbol, 200.59);
    }

    protected Mercury(String symbol, double mass) {
        super(symbol, mass, 80);
    }

    public static final Mercury INSTANCE = new Mercury("Hg");

    private static final long serialVersionUID = 1L;
}

