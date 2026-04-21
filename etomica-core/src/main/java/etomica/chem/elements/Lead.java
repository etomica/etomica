package etomica.chem.elements;

public class Lead extends ElementChemical {

    protected Lead(String symbol) {
        this(symbol, 207);
    }

    protected Lead(String symbol, double mass) {
        super(symbol, mass, 82);
    }

    public static final Lead INSTANCE = new Lead("Pb");

    private static final long serialVersionUID = 1L;
}
