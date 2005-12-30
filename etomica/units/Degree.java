package etomica.units;

import java.io.ObjectStreamException;

public final class Degree extends SimpleUnit {

  /**
   * Single instance of this unit.
   */
    public static final Degree UNIT = new Degree();
    
    private Degree() {
        super(Angle.DIMENSION,
        	Math.PI/180., //conversion from degrees to radians
        	"degrees","\u00B0", Prefix.ALLOWED //unicode for the degree symbol
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