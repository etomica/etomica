package etomica.chem.elements;

public class Cesium extends ElementChemical {

    protected Cesium(String symbol) {
        this(symbol, 132.905);
    }

    protected Cesium(String symbol, double mass) {
        super(symbol, mass, 55);
    }

    public static final Cesium INSTANCE = new Cesium("Cs");

    private static final long serialVersionUID = 1L;
}
