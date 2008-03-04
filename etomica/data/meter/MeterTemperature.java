package etomica.data.meter;

import etomica.EtomicaInfo;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeMoleculeOriented;
import etomica.atom.AtomTypeOrientedSphere;
import etomica.atom.IMolecule;
import etomica.atom.MoleculeOrientedDynamic;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.simulation.ISimulation;
import etomica.species.ISpecies;
import etomica.units.Dimension;
import etomica.units.Temperature;

/**
 * Meter for measurement of the temperature based on kinetic-energy
 * equipartition.  The class uses a MeterKineticEnergy by default to calculate
 * the kinetic energy, but any DataSourceScalar can be used for this purpose by
 * calling setKineticEnergyMeter.
 * 
 * If the ISimulation is not given, the class will assume that all atoms have
 * only translational degrees of freedom.  If the ISimulation is given, this
 * class will examine the ISpecies and calculate the actual number of degrees
 * of freedom (more for oriented atoms or molecules).
 * 
 * @author Andrew Schultz
 */
public class MeterTemperature extends DataSourceScalar {

    public MeterTemperature(Box box, int D) {
        this(null, box, D);
    }

    public MeterTemperature(ISimulation sim, Box box, int D) {
		super("Temperature", Temperature.DIMENSION);
		dim = D;
		meterKE = new MeterKineticEnergy();
		((MeterKineticEnergy)meterKE).setBox(box);
		this.sim = sim;
		this.box = box;
	}
    
    public void setKineticEnergyMeter(DataSourceScalar meterKineticEnergy) {
        meterKE = meterKineticEnergy;
    }

	public static EtomicaInfo getEtomicaInfo() {
		EtomicaInfo info = new EtomicaInfo(
				"Records temperature as given via kinetic energy");
		return info;
	}

	public double getDataAsScalar() {
	    int totalD = box.atomCount() * dim;
	    if (sim != null) {
	        totalD = 0;
	        ISpecies[] species = sim.getSpeciesManager().getSpecies();
	        for (int i=0; i<species.length; i++) {
	            int nMolecules = box.getNMolecules(species[i]);
	            if (nMolecules > 0) {
	                IMolecule molecule = (IMolecule)box.getMoleculeList(species[i]).getAtom(0);
	                if (molecule instanceof MoleculeOrientedDynamic) {
	                    if (Double.isInfinite(((AtomTypeMoleculeOriented)species[i].getMoleculeType()).getMass())) {
	                        continue;
	                    }
                        totalD += 6*nMolecules;
	                }
	                else {
	                    AtomSet children = molecule.getChildList();
	                    if (children.getAtomCount() == 0 || 
	                        Double.isInfinite(((AtomTypeLeaf)children.getAtom(0).getType()).getMass())) {
	                        continue;
	                    }
	                    if (children.getAtom(0).getType() instanceof AtomTypeOrientedSphere) {
	                        // oriented sphere at this point corresponds to cylindrical symmetry
	                        totalD += 5*nMolecules*children.getAtomCount();
	                    }
	                    else {
	                        totalD += 3*nMolecules*children.getAtomCount();
	                    }
	                }
	            }
	        }
	    }
		return (2. / totalD) * meterKE.getDataAsScalar();
	}

	public Dimension getDimension() {
		return Temperature.DIMENSION;
	}

    private static final long serialVersionUID = 1L;
    protected Box box;
	protected DataSourceScalar meterKE;
	protected final ISimulation sim;
	protected final int dim;
}