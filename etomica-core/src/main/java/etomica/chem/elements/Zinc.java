package etomica.chem.elements;

public class Zinc extends ElementChemical {

	protected Zinc(String symbol) {this(symbol, 65.39);}

    protected Zinc(String symbol, double mass) {
        super(symbol, mass, 30);
    }

    public static final Zinc INSTANCE = new Zinc("Zn");

    private static final long serialVersionUID = 1L;
}
