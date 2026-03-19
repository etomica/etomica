package etomica.chem.elements;

public class Yttrium extends ElementChemical {
    protected Yttrium(String symbol) {
        this(symbol, 88.906);
    }
    protected Yttrium(String symbol, double mass) {
        super(symbol, mass, 39);
    }

    public static final Yttrium INSTANCE = new Yttrium("Y");

    private static final long serialVersionUID = 1L;
}
