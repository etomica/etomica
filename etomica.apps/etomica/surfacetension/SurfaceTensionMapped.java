package etomica.surfacetension;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.DataPipe;
import etomica.data.DataProcessor;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterProfileByVolume;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray;
import etomica.integrator.Integrator.Forcible;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.PotentialCalculationForceSum;
import etomica.space.ISpace;
import etomica.units.Null;


/**
 * DataProcessor that adds in mapped average contributions to the surface
 * tension.
 * 
 * @author Andrew Schultz
 */
public class SurfaceTensionMapped extends DataProcessor implements AgentSource<Forcible> {

    protected final DataInfoDouble dataInfo = new DataInfoDouble("surface tension", Null.DIMENSION);
    protected final DataDouble data = new DataDouble();
    protected final FitTanh fit;
    protected final MeterProfileByVolume densityProfileMeter;
    protected final IPotentialMaster potentialMaster;
    protected final PotentialCalculationForceSum pcForce;
    protected final AtomLeafAgentManager<Forcible> forceAgentManager;
    protected final ISpace space;
    protected final IBox box;
    protected final IteratorDirective allAtoms;

    public SurfaceTensionMapped(ISpace space, IBox box, ISpecies species, IPotentialMaster potentialMaster) {
        this.space = space;
        this.box = box;
        densityProfileMeter = new MeterProfileByVolume(space);
        densityProfileMeter.setBox(box);
        densityProfileMeter.getXDataSource().setNValues(400);
        MeterNMolecules meterNMolecules = new MeterNMolecules();
        meterNMolecules.setSpecies(species);
        densityProfileMeter.setDataSource(meterNMolecules);

        fit = new FitTanh();
        
        this.potentialMaster = potentialMaster;
        pcForce = new PotentialCalculationForceSum();
        forceAgentManager = new AtomLeafAgentManager<Forcible>(this, box, Forcible.class);
        pcForce.setAgentManager(forceAgentManager);
        allAtoms = new IteratorDirective();
    }
    
    public DataPipe getDataCaster(IEtomicaDataInfo inputDataInfo) {
        return null;
    }

    protected IEtomicaDataInfo processDataInfo(IEtomicaDataInfo inputDataInfo) {
        return dataInfo;
    }

    protected IData processData(IData inputData) {
        double st = inputData.getValue(0);
        
        double[] x = ((DataDoubleArray)densityProfileMeter.getXDataSource().getData()).getData();
        double[] y = ((DataDoubleArray)densityProfileMeter.getData()).getData();
        double[] param = fit.doFit(x, y);
//        System.out.println("fit: "+Arrays.toString(param));
        IAtomList atoms = box.getLeafList();
        double rho = atoms.getAtomCount()/box.getBoundary().volume();
        double L = box.getBoundary().getBoxSize().getX(0);
        double rhoL2tanh = (Math.tanh(2 * (L/2 - param[3]) / param[2]) + Math.tanh(2 * (L/2 + param[3]) / param[2]));
        double rhoL2 = param[0] + 0.5 * (param[1] - param[0]) * rhoL2tanh;
        double dzsdL = (rho-rhoL2)/(param[1]-param[0])/rhoL2tanh;

        potentialMaster.calculate(box, allAtoms, pcForce);
        double mapSum = 0;
        for (int i=0; i<atoms.getAtomCount(); i++) {
            IAtom atom = atoms.getAtom(i);
            IVector f = forceAgentManager.getAgent(atom).force();
            double px = atom.getPosition().getX(0);
            double t1 = Math.tanh(2 * (px + param[3]) / param[2]);
            double t2 = Math.tanh(2 * (px - param[3]) / param[2]);
            double irho = param[0] + 0.5 * (param[1] - param[0]) * (t1 - t2);
            double foo = param[0]*px/L - (param[1]-param[0])*param[2]/(4*L) * Math.log(Math.cosh(2*(px-param[3])/param[2])/Math.cosh(2*(px-param[3])/param[2]));
            double bar = dzsdL * 0.5 * (param[1]-param[0]) * (t1 + t2);
            double xL = (foo - bar)/irho; 
            mapSum += (xL - px/L)*f.getX(0);
            
            // still needs Jacobian
        }
        
        data.x = st + mapSum/box.getBoundary().volume();
        return data;
    }

    public Forcible makeAgent(IAtom a) {
        return new IntegratorVelocityVerlet.MyAgent(space);
    }

    public void releaseAgent(Forcible agent, IAtom atom) {}

}
