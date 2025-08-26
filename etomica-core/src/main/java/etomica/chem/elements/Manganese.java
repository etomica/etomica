package etomica.chem.elements;

public class Manganese extends ElementChemical {

	protected Manganese(String symbol) {this(symbol, 54.94);}

    protected Manganese(String symbol, double mass) {
        super(symbol, mass, 25);
    }

    public static final Manganese INSTANCE = new Manganese("Mn");

    private static final long serialVersionUID = 1L;
}
