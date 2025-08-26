package etomica.chem.elements;

public class Magnesium extends ElementChemical{
    protected Magnesium(String symbol){this(symbol,24.305);}

    public Magnesium(String symbol, double mass) {
        super(symbol, mass, 12);
    }

    private static final long serialVersionUID = 1L;
    public static final Magnesium INSTANCE = new Magnesium("Mg");
}
