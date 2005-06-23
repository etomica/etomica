package etomica.graphics;
import etomica.data.DataPump;
import etomica.data.DataSourceCountTime;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntervalActionAdapter;
import etomica.units.Picosecond;
import etomica.units.PrefixedUnit;

/** 
 * DisplayBox to present the elapsed time in a
 * molecular dynamics simulation.  Designed for use with
 * a single integrator.
 */
public class DisplayTimer extends DisplayBox {

    public DisplayTimer(IntegratorMD integrator) {
        timer = new DataSourceCountTime();
        DataPump dataPump = new DataPump(timer, this);
        intervalActionAdapter = new IntervalActionAdapter(dataPump, integrator);
        setUpdateInterval(100);
        this.setUnit(new PrefixedUnit(Picosecond.UNIT));
        this.setPrecision(7);
        graphic().setSize(100,60);
    }
    
    /**
     * Sets the period for updating the display.  Number of integrator 
     * interval events between updates.  Does not have any effect on the
     * value displayed; affects only how often it is updated.
     */
    public void setUpdateInterval(int interval) {
        intervalActionAdapter.setActionInterval(interval);
    }

    /**
     * Returns the data source used to count the time, to permit
     * access to its methods for reset, etc.  
     * @return
     */
    public DataSourceCountTime getTimer() {
        return timer;
    }
    
    private final IntervalActionAdapter intervalActionAdapter;
    private final DataSourceCountTime timer;
}
