package etomica.data.meter;

import etomica.EtomicaInfo;
import etomica.data.DataSourceScalar;
import etomica.box.Box;
import etomica.units.Dimension;
import etomica.units.Temperature;

/**
 * Meter for measurement of the temperature based on kinetic-energy
 * equipartition
 */

public class MeterTemperature extends DataSourceScalar {

    public MeterTemperature() {
		super("Temperature", Temperature.DIMENSION);
		meterKE = new MeterKineticEnergy();
	}

	public static EtomicaInfo getEtomicaInfo() {
		EtomicaInfo info = new EtomicaInfo(
				"Records temperature as given via kinetic energy");
		return info;
	}

	public double getDataAsScalar() {
        if (box == null) throw new IllegalStateException("must call setBox before using meter");
		return (2. / (box.atomCount() * box.getSpace().D()))
				* meterKE.getDataAsScalar();
	}

	public Dimension getDimension() {
		return Temperature.DIMENSION;
	}

    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }
    /**
     * @param box The box to set.
     */
    public void setBox(Box box) {
        this.box = box;
        meterKE.setBox(box);
    }

    private static final long serialVersionUID = 1L;
    protected Box box;
	protected final MeterKineticEnergy meterKE;
}