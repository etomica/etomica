package etomica.chem.elements;

public class Calcium extends ElementChemical {

    protected Calcium(String symbol) {
        this(symbol, 40.078);
    }

    protected Calcium(String symbol, double mass) {
        super(symbol, mass, 20);
    }

    public static final Calcium INSTANCE = new Calcium("Ca");

}
