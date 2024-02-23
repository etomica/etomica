package etomica.chem.elements;

public class Ruthenium extends ElementChemical {

    protected Ruthenium(String symbol) {
        this(symbol, 101.10);
    }

    protected Ruthenium(String symbol, double mass) {
        super(symbol, mass, 44);
    }

    public static final Ruthenium INSTANCE = new Ruthenium("Ru");

    private static final long serialVersionUID = 1L;
}
