package etomica.chem.elements;

public class Krypton extends ElementChemical{
    protected Krypton(String symbol) {
        this(symbol, 83.798);
    }

    protected Krypton(String symbol, double mass) {
        super(symbol, mass, 36);
    }

    public static final Krypton INSTANCE = new Krypton("Kr");

    private static final long serialVersionUID = 1L;
}
