package etomica.units.systems;
import etomica.Default;
import etomica.units.BaseUnit;
import etomica.units.Count;
import etomica.units.Decimal;
import etomica.units.Prefix;
import etomica.units.Radian;
import etomica.units.Unit;
import etomica.units.UnitSystem;

/**
 * Lennard-Jones system of units, such that all quantities are made
 * dimensionless with respect to a characteristic size (sigma), energy (epsilon)
 * and mass. Values of sigma, epsilon and mass may be set to any value using the
 * static accessor methods (set/get).  Changes are immediately propagated to the
 * sim- >LJ unit conversion factors.
 */
 
 /* History
  * 03/12/04 (DAK) removed most static features, allowing for different
  * instances
  */
  
public class LJ extends UnitSystem implements java.io.Serializable {
    
    private double sigma = Default.ATOM_SIZE;
    private double epsilon = Default.POTENTIAL_WELL;
    private double mass = Default.ATOM_MASS;
    
    private Mass massUnit = new Mass();
    private Length lengthUnit = new Length();
    private Time timeUnit = new Time();
    private Charge chargeUnit = new Charge();
    private Dipole dipoleUnit = new Dipole();
    private Energy energyUnit = new Energy();
    private Temperature temperatureUnit = new Temperature();
    private Pressure2D pressure2DUnit = new Pressure2D();
    private Pressure pressureUnit = new Pressure();
    private Volume2D volume2DUnit = new Volume2D();
    private Volume volumeUnit = new Volume();
        
    public double getSigma() {return sigma;}
    public void setSigma(double s) {sigma = s; update();}
    public double getEpsilon() {return epsilon;}
    public void setEpsilon(double e) {epsilon = e; update();}
    public double getMass() {return mass;}
    public void setMass(double m) {mass = m; update();}
    
	public Unit quantity() {return Count.UNIT;}
    public Unit fraction() {return Decimal.UNIT;}
	public Unit mass() {return massUnit;}
	public Unit length() {return lengthUnit;}
	public Unit time() {return timeUnit;}
	public Unit angle() {return Radian.UNIT;}
	public Unit charge() {return chargeUnit;}  
	public Unit dipole() {return dipoleUnit;}  
	public Unit energy() {return energyUnit;}
	public Unit temperature() {return energyUnit;}
	public Unit pressure(int D) {return (D==2) ? (Unit)pressure2DUnit : (Unit)pressureUnit;}  
	public Unit volume(int D) {return (D==2) ? (Unit)volume2DUnit : (Unit)volumeUnit;} 

        
    private void update() {
        massUnit.updateConversions();
        lengthUnit.updateConversions();
        timeUnit.updateConversions();
        chargeUnit.updateConversions();
        dipoleUnit.updateConversions();
        energyUnit.updateConversions();
        temperatureUnit.updateConversions();
        pressureUnit.updateConversions();
        pressure2DUnit.updateConversions();
        volumeUnit.updateConversions();
        volume2DUnit.updateConversions();
    }
    
    //  \u03B5 is unicode for epsilon
    //  \u03C3 is unicode for sigma
    //  \u00BD is unicode for "1/2"
    
    public final class Mass extends BaseUnit.Mass {
        private Mass() {
        	super(1.0, "LJ units", "m", Prefix.NOT_ALLOWED);
        	updateConversions();
        }
        private void updateConversions() {
            setToSimConversionFactor(mass); //conversion from LJ unit to simulation unit
        }
    }

    public final class Length extends BaseUnit.Length {
        private Length() {
			super(1.0, "LJ units", "\u03C3", Prefix.NOT_ALLOWED);
			updateConversions();
        }
        private void updateConversions() {
			setToSimConversionFactor(sigma); //conversion from LJ unit to simulation unit
        }
    }

    public final class Time extends BaseUnit.Time {
        private Time() {
			super(1.0, "LJ units", "(\u03B5m\u03C3^2)^\u00BD", Prefix.NOT_ALLOWED);
            updateConversions();
        }
        private void updateConversions() {
			setToSimConversionFactor(1.0/Math.sqrt(epsilon * mass * sigma*sigma)); //conversion from LJ unit to simulation unit
        }
    }

    public final class Charge extends BaseUnit.Charge {
        private Charge() {
			super(1.0, "LJ units", "(\u03B5\u03C3)^(-\u00BD)", Prefix.NOT_ALLOWED);
            updateConversions();
        }
        private void updateConversions() {
			setToSimConversionFactor(Math.sqrt(epsilon * sigma)); //conversion from LJ unit to simulation unit
        }
    }

    public final class Dipole extends BaseUnit.Dipole {
        private Dipole() {
			super(1.0, "LJ units", "(\u03B5\u03C3^3)^(-\u00BD)", Prefix.NOT_ALLOWED);
            updateConversions();
        }
        private void updateConversions() {
			setToSimConversionFactor(Math.sqrt(epsilon * sigma*sigma*sigma)); //conversion from LJ unit to simulation unit
        }
    }

    public final class Energy extends BaseUnit.Energy {
        private Energy() {
			super(1.0, "LJ units", "\u03B5", Prefix.NOT_ALLOWED);
            updateConversions();
        }
        private void updateConversions() {
			setToSimConversionFactor(epsilon); //conversion from LJ unit to simulation unit
        }
    }


    public final class Temperature extends BaseUnit.Temperature {
        private Temperature() {
			super(1.0, "LJ units", "\u03B5", Prefix.NOT_ALLOWED);
            updateConversions();
        }
        private void updateConversions() {
			setToSimConversionFactor(epsilon); //conversion from LJ unit to simulation unit
        }
    }
    public final class Pressure extends BaseUnit.Pressure {
        private Pressure() {
			super(1.0, "LJ units", "\u03B5/\u03C3^3", Prefix.NOT_ALLOWED);
            updateConversions();
        }
        private void updateConversions() {
			setToSimConversionFactor(epsilon/(sigma*sigma*sigma)); //conversion from LJ unit to simulation unit
        }
    }

    public final class Pressure2D extends BaseUnit.Pressure2D {
        private Pressure2D() {
			super(1.0, "LJ units", "\u03B5/\u03C3^2", Prefix.NOT_ALLOWED);
            updateConversions();
        }
        private void updateConversions() {
			setToSimConversionFactor(epsilon/(sigma*sigma)); //conversion from LJ unit to simulation unit
        }
    }

    public final class Volume extends BaseUnit.Volume {
        private Volume() {
			super(1.0, "LJ units", "\u03C3^3", Prefix.NOT_ALLOWED);
            updateConversions();
        }
        private void updateConversions() {
			setToSimConversionFactor(sigma*sigma*sigma); //conversion from LJ unit to simulation unit
        }
    }

    public final class Volume2D extends BaseUnit.Volume2D {
        private Volume2D() {
			super(1.0, "LJ units", "\u03C3^2", Prefix.NOT_ALLOWED);
            updateConversions();
        }
        private void updateConversions() {
			setToSimConversionFactor(sigma*sigma); //conversion from LJ unit to simulation unit
        }
    }
}