package etomica.chem.elements;

public class Iridium extends ElementChemical {

	protected Iridium(String symbol) {this(symbol, 192.2);}

    protected Iridium(String symbol, double mass) {
        super(symbol, mass, 77);
    }

    public static final Iridium INSTANCE = new Iridium("Ir");

    private static final long serialVersionUID = 1L;
}
