package etomica.graphics;
import etomica.integrator.IntegratorMD;
import etomica.units.Picosecond;
import etomica.units.PrefixedUnit;

/** 
 * DisplayBox to present the elapsed time in a
 * molecular dynamics simulation.
 */
public class DisplayTimer extends DisplayBox {

    public DisplayTimer(IntegratorMD integrator) {
        this(integrator.chronoMeter());
    }
    public DisplayTimer(IntegratorMD.ChronoMeter meter) {
        super(meter.integrator().simulation());
        this.setMeter(meter);
        this.setUnit(new PrefixedUnit(Picosecond.UNIT));
        this.setPrecision(7);
        graphic().setSize(100,60);
    }
}
