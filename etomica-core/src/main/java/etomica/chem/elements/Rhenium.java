package etomica.chem.elements;

public class Rhenium  extends ElementChemical {

    protected Rhenium(String symbol) {
        this(symbol, 186.21);
    }

    protected Rhenium(String symbol, double mass) {
        super(symbol, mass, 75);
    }

    public static final Rhenium INSTANCE = new Rhenium("Re");

    private static final long serialVersionUID = 1L;
}
