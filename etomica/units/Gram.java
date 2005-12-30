package etomica.units;
import java.io.ObjectStreamException;

import etomica.util.Constants;

public final class Gram extends SimpleUnit {

  /**
   * Singleton instance of this unit.
   */
    public static final Gram UNIT = new Gram();
    
    private Gram() {
        super(Mass.DIMENSION,
        	Constants.AVOGADRO, //conversion from grams to Daltons
        	"grams", "g", Prefix.ALLOWED
        	);   
    }
    
    /**
     * Required to guarantee singleton when deserializing.
     * 
     * @return the singleton UNIT
     */
    private Object readResolve() throws ObjectStreamException {
        return UNIT;
    }
    
    private static final long serialVersionUID = 1;

}