package etomica.units;

public final class Second extends BaseUnit.Time {

  /**
   * Singleton instance of this unit.
   */
    public static final Second UNIT = new Second();
    
    private Second() {
        super(
        	1e+12, //conversion from seconds to picoseconds
        	"seconds",
        	"s"
        	);   
    }
}