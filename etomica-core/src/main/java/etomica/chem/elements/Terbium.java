package etomica.chem.elements;

public class Terbium extends ElementChemical {

    protected Terbium(String symbol) {
        this(symbol, 158.93);
    }

    protected Terbium(String symbol, double mass) {
        super(symbol, mass, 65);
    }

    public static final Platinum INSTANCE = new Platinum("Tb");

    private static final long serialVersionUID = 1L;
}