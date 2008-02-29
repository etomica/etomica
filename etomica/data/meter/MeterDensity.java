package etomica.data.meter;

import etomica.EtomicaInfo;
import etomica.data.DataSourceScalar;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.units.DimensionRatio;
import etomica.units.Quantity;
import etomica.units.Volume;

/**
 * Meter for measurement of the total molecule number density in a box
 * Molecule number density is defined (number of molecules)/(volume of box)
 */
public class MeterDensity extends DataSourceScalar {
    
    public MeterDensity(Space space) {
        super("Number Density",new DimensionRatio(Quantity.DIMENSION, Volume.dimension(space.D())));
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Number density (molecules/volume) in a box");
        return info;
    }
    
    public void setSpecies(ISpecies s) {
        species = s;
    }
    public ISpecies getSpecies() {
    	return species;
    }

    public double getDataAsScalar() {
        if (box == null) throw new IllegalStateException("must call setBox before using meter");
        return (species == null ? 
        			box.moleculeCount() : 
        			box.getNMolecules(species))
				/box.volume();
    }
    
    /**
     * @return Returns the box.
     */
    public IBox getBox() {
        return box;
    }
    /**
     * @param box The box to set.
     */
    public void setBox(Box box) {
        this.box = box;
    }

    private static final long serialVersionUID = 1L;
    private Box box;
    private ISpecies species;
}
