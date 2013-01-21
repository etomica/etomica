package etomica.chem.elements;

/**
 * Class for the Fluoride element.
 * reference: Atomic weights of the elements 2009 (IUPAC Technical Report)
 * http://iupac.org/publications/pac/83/2/0359/
 * 
 * @author Shu
 */
public class Fluoride extends ElementChemical {

    protected Fluoride(String symbol) {
        this(symbol, 18.9984);
    }
    
    public Fluoride(String symbol, double mass) {
        super(symbol, mass, 9);
    }

    private static final long serialVersionUID = 1L;
    public static final Fluoride INSTANCE = new Fluoride("F");
}
