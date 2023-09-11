package etomica.chem.elements;

public class Silicon extends ElementChemical{
    protected Silicon(String symbol){this(symbol,28.086);}

    public Silicon(String symbol, double mass) {
        super(symbol, mass, 12);
    }

    private static final long serialVersionUID = 1L;
    public static final Silicon INSTANCE = new Silicon("Si");
}
