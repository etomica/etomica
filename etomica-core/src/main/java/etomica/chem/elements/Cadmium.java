package etomica.chem.elements;

public class Cadmium extends ElementChemical {
    protected Cadmium(String symbol) {this(symbol, 112.414);}
    protected Cadmium(String symbol, double mass) {
        super(symbol, mass, 48);
    }

    public static final Cadmium INSTANCE = new Cadmium("Cd");

}
