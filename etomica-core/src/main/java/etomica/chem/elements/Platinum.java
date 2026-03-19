package etomica.chem.elements;


public class Platinum extends ElementChemical {

    protected Platinum(String symbol) {
        this(symbol, 195.08);
    }

    protected Platinum(String symbol, double mass) {
        super(symbol, mass, 78);
    }

    public static final Platinum INSTANCE = new Platinum("Pt");

    private static final long serialVersionUID = 1L;
}
