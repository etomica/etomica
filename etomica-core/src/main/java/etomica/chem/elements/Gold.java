package etomica.chem.elements;

public class Gold extends ElementChemical {

    protected Gold(String symbol) {
        this(symbol, 196.966);
    }

    public Gold(String symbol, double mass) {
        super(symbol, mass, 79);
    }

    private static final long serialVersionUID = 1L;
    public static final Gold INSTANCE = new Gold("Au");
}
