package etomica.chem.elements;

public class Samarium extends ElementChemical {

    protected Samarium(String symbol) {
        this(symbol, 150.36);
    }

    protected Samarium(String symbol, double mass) {
        super(symbol, mass, 62);
    }

    public static final Samarium INSTANCE = new Samarium("Sm");

    private static final long serialVersionUID = 1L;
}