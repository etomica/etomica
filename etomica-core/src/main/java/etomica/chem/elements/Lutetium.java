package etomica.chem.elements;

public class Lutetium extends ElementChemical {

    protected Lutetium(String symbol) {
        this(symbol, 174.97);
    }

    protected Lutetium(String symbol, double mass) {
        super(symbol, mass, 71);
    }

    public static final Lutetium INSTANCE = new Lutetium("Lu");

    private static final long serialVersionUID = 1L;
}
