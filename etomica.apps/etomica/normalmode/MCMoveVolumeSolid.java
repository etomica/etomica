package etomica.normalmode;

import etomica.action.BoxInflate;
import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.space.ISpace;
import etomica.units.Dimension;
import etomica.units.Pressure;
import etomica.util.Function;

/**
 * Monte Carlo volume-change move for simulations of crystalline solids in the
 * NPT ensemble.  When changing the volume, atom coordinates are scaled away
 * from or toward their lattice in order to improve the likelihood of
 * acceptance.  
 *
 * @author David Kofke, Andrew Schultz
 */
public class MCMoveVolumeSolid extends MCMoveBoxStep {
    
    private static final long serialVersionUID = 2L;
    protected double pressure;
    private MeterPotentialEnergy energyMeter;
    protected final BoxInflate inflate;
    private final int D;
    private IRandom random;
    protected final AtomIteratorLeafAtoms affectedAtomIterator;
    protected double temperature;
    protected final CoordinateDefinition coordinateDefinition;
    protected final IVectorMutable dr;
    protected final IVectorMutable nominalBoxSize;

    private transient double uOld, vOld, vNew, vScale;
    private transient double uNew = Double.NaN, latticeScale;
    
    protected Function uLatFunction = uLat0;

    /**
     * @param potentialMaster an appropriate PotentialMaster instance for calculating energies
     * @param space the governing space for the simulation
     */
    public MCMoveVolumeSolid(IPotentialMaster potentialMaster, CoordinateDefinition coordinateDefinition, IRandom random,
    		            ISpace _space, double pressure) {
        super(potentialMaster);
        this.coordinateDefinition = coordinateDefinition;
        this.random = random;
        this.D = _space.D();
        nominalBoxSize = _space.makeVector();
        nominalBoxSize.E(coordinateDefinition.getBox().getBoundary().getBoxSize());
        dr = _space.makeVector();
        inflate = new BoxInflate(_space);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(0.1);
        setStepSizeMin(0.0);
        setStepSize(0.01);
        setPressure(pressure);
        energyMeter.setIncludeLrc(true);
        affectedAtomIterator = new AtomIteratorLeafAtoms();
    }
    
    public void setBox(IBox p) {
        super.setBox(p);
        energyMeter.setBox(p);
        inflate.setBox(p);
        affectedAtomIterator.setBox(p);
    }

    /**
     * Sets the temperature being sampled.  The temperature is used to
     * determine how to scale the atom coordinates.
     * 
     * In actuality, only P/kT is important, but we'll keep the methods
     * separate. 
     */
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }

    /**
     * Sets a function that returns the lattice energy for a given density.
     * If not set, the default lattice energy is taken to be 0 for all densities.
     */
    public void setULatFunction(Function newULatFunction) {
        uLatFunction = newULatFunction;
    }
    
    public Function getULatFunction() {
        return uLatFunction;
    }
    
    public boolean doTrial() {
        vOld = box.getBoundary().volume();
        uOld = energyMeter.getDataAsScalar();
        vScale = (2.*random.nextDouble()-1.)*stepSize;
        vNew = vOld * Math.exp(vScale); //Step in ln(V)
        double rScale = Math.exp(vScale/D);
        inflate.setScale(rScale);
        inflate.actionPerformed();

        IAtomList leafList = box.getLeafList();
        int nAtoms = leafList.getAtomCount();
        double uLatOld = nAtoms*uLatFunction.f(nAtoms/vOld);
        double uLatNew = nAtoms*uLatFunction.f(nAtoms/vNew);

        latticeScale = Math.exp((pressure*(vNew-vOld)+(uLatNew-uLatOld))/(nAtoms*temperature*D))/rScale;

        IVector boxSize = box.getBoundary().getBoxSize();
        for (int i=0; i<leafList.getAtomCount(); i++) {
            IAtom atomi = leafList.getAtom(i);
            IVector site = coordinateDefinition.getLatticePosition(atomi);
            dr.E(site);
            dr.TE(boxSize);
            dr.DE(nominalBoxSize);
            // dr is now the new lattice site
            dr.TE(1-latticeScale);
            atomi.getPosition().TE(latticeScale);
            atomi.getPosition().PE(dr);
        }
        
        uNew = energyMeter.getDataAsScalar();
        return true;
    }
    
    public double getA() {
        return vNew/vOld;
    }
    
    public double getB() {
        int nAtoms = box.getLeafList().getAtomCount();
        double uLatOld = nAtoms*uLatFunction.f(nAtoms/vOld);
        double uLatNew = nAtoms*uLatFunction.f(nAtoms/vNew);
        return -((uNew-uLatNew) - (uOld-uLatOld));
    }
    
    public void acceptNotify() {
        /* do nothing */
    }
    
    public void rejectNotify() {
        IAtomList leafList = box.getLeafList();
        IVector boxSize = box.getBoundary().getBoxSize();
        latticeScale = 1.0 / latticeScale;
        for (int i=0; i<leafList.getAtomCount(); i++) {
            IAtom atomi = leafList.getAtom(i);
            IVector site = coordinateDefinition.getLatticePosition(atomi);
            dr.E(site);
            dr.TE(boxSize);
            dr.DE(nominalBoxSize);
            // dr is now the new lattice site
            dr.TE(1.0-latticeScale);
            atomi.getPosition().TE(latticeScale);
            atomi.getPosition().PE(dr);
        }
        
        inflate.undo();
    }

    public double energyChange() {return uNew - uOld;}
    
    public AtomIterator affectedAtoms() {
        return affectedAtomIterator;
    }

    public void setPressure(double p) {pressure = p;}
    public final double getPressure() {return pressure;}
    public Dimension getPressureDimension() {return Pressure.DIMENSION;}

    /**
     * Nominal function for lattice energy
     */
    public static Function uLat0 = new Function() {
        public double f(double x) {return 0;}
    };
}