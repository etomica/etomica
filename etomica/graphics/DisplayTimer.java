package etomica.graphics;
import etomica.data.DataPump;
import etomica.data.DataSourceCountTime;
import etomica.integrator.IntegratorMD;
import etomica.units.Prefix;
import etomica.units.PrefixedUnit;
import etomica.units.Second;

/** 
 * DisplayBox to present the elapsed time in a
 * molecular dynamics simulation.  Designed for use with
 * a single integrator.
 */
public class DisplayTimer extends DisplayBox {

    public DisplayTimer(IntegratorMD integrator) {
        this(integrator, new DataSourceCountTime(integrator));
    }
        
    private DisplayTimer(IntegratorMD integrator, DataSourceCountTime timer) {
        super(timer.getDataInfo());
        this.timer = timer;
        this.integrator = integrator;
        dataPump = new DataPump(timer, this);
        integrator.addIntervalAction(dataPump);
        setUpdateInterval(100);
        this.setUnit(new PrefixedUnit(Second.UNIT));
        ((PrefixedUnit)unit).setPrefix(Prefix.PICO);
        this.setPrecision(7);
        graphic().setSize(100,60);
    }
    
    /**
     * Sets the period for updating the display.  Number of integrator 
     * interval events between updates.  Does not have any effect on the
     * value displayed; affects only how often it is updated.
     */
    public void setUpdateInterval(int interval) {
        integrator.setActionInterval(dataPump, interval);
    }

    /**
     * Unhooks the DisplayTimer from the integrator
     */
    public void dispose() {
        integrator.removeIntervalListener(dataPump);
    }
    
    /**
     * Returns the data source used to count the time, to permit
     * access to its methods for reset, etc.  
     * @return
     */
    public DataSourceCountTime getTimer() {
        return timer;
    }
    
    protected final DataPump dataPump;
    protected final IntegratorMD integrator;
    protected final DataSourceCountTime timer;
}
