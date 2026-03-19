package etomica.chem.elements;

public class Indium extends ElementChemical {

    protected Indium(String symbol) {
        this(symbol, 114.818);
    }

    protected Indium(String symbol, double mass) {
        super(symbol, mass, 49);
    }

    public static final Indium INSTANCE = new Indium("In");

    private static final long serialVersionUID = 1L;
}
