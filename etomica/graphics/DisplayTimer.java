package etomica.graphics;
import etomica.DataSource;
import etomica.data.DataSourceCountTime;
import etomica.integrator.IntegratorMD;
import etomica.units.Picosecond;
import etomica.units.PrefixedUnit;

/** 
 * DisplayBox to present the elapsed time in a
 * molecular dynamics simulation.
 */
public class DisplayTimer extends DisplayBox {

    public DisplayTimer(IntegratorMD integrator) {
        this(new DataSourceCountTime(integrator));
    }
    public DisplayTimer(DataSource dataSource) {
        super(dataSource);
        this.setUnit(new PrefixedUnit(Picosecond.UNIT));
        this.setPrecision(7);
        graphic().setSize(100,60);
    }
}
