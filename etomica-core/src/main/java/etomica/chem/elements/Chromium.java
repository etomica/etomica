package etomica.chem.elements;

public class Chromium extends ElementChemical {

	protected Chromium(String symbol) {this(symbol, 52.00);}

    protected Chromium(String symbol, double mass) {
        super(symbol, mass, 24);
    }

    public static final Chromium INSTANCE = new Chromium("Cr");

    private static final long serialVersionUID = 1L;
}
