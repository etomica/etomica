package etomica.chem.elements;

public class Selenium extends ElementChemical {

    protected Selenium(String symbol) {this(symbol, 78.96);}
    protected Selenium(String symbol, double mass) {
        super(symbol, mass, 34);
    }

    public static final Selenium INSTANCE = new Selenium("Se");

    private static final long serialVersionUID = 1L;
}
