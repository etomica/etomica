package etomica.modules.interfacial;

import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IPotential;
import etomica.api.IVector;
import etomica.integrator.IntegratorBox;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialSoft;
import etomica.space.ISpace;
import etomica.space.Tensor;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.  Additionally, this class has the potential
 * calculate the pressureTensor (which can be done efficiently during the
 * gradient calculation).
 */
public class PotentialCalculationForcePressureBinSum extends PotentialCalculationForceSum {

    public PotentialCalculationForcePressureBinSum(ISpace space) {
        this.space = space;
        virialTensor = space.makeTensor();
        virialTensorTot = space.makeTensor();
        dataCountDown = 10;
    }
    
    public void setCutoff(double newCutoff) {
        cutoff = newCutoff;
    }
    
    public double getCutoff() {
        return cutoff;
    }
    
    public void setBinSize(double newBinSize) {
        binSize = newBinSize;
        if (box != null) {
            double l = box.getBoundary().getDimensions().x(0);
            nBins = (int)(l/binSize);
            binSize = l / nBins;
            densityProfile = new double[nBins];
            virialTensorProfile = new Tensor[nBins];
            uLrc = new double[(nBins+1)/2];
            fLrc = new double[(nBins+1)/2];
            xVirialLrc = new double[(nBins+1)/2];
            yzVirialLrc = new double[(nBins+1)/2];
            for (int i=0; i<nBins; i++) {
                virialTensorProfile[i] = space.makeTensor();
            }
        }
    }
    
    public void setBox(IBox newBox) {
        box = newBox;
        if (binSize != 0) {
            double l = box.getBoundary().getDimensions().x(0);
            nBins = (int)(l/binSize);
            binSize = l / nBins;
            densityProfile = new double[nBins];
            uLrc = new double[(nBins+1)/2];
            fLrc = new double[(nBins+1)/2];
            xVirialLrc = new double[(nBins+1)/2];
            yzVirialLrc = new double[(nBins+1)/2];
            virialTensorProfile = new Tensor[nBins];
            for (int i=0; i<nBins; i++) {
                virialTensorProfile[i] = space.makeTensor();
            }
        }
    }

    /**
     * Zeros out the pressureTensor.  This method should be called before
     * invoking potentialMaster.calculate so that the pressureTensor is
     * correct at the end of the calculation.
     */
    public void reset() {
        super.reset();
        virialTensorTot.E(0);
        if (dataCountDown == 0) {
            dataCountDown = 10;
        }
        dataCountDown--;

        IVector dim = box.getBoundary().getDimensions();
        L = dim.x(0);
        halfL = 0.5 * dim.x(0);
        Li = 0.5 / halfL;
        nBins = (int)(L/binSize);

        if (dataCountDown == 0 || nBins != densityProfile.length) {
            if (nBins != densityProfile.length || uLrc[0] == 0) {
                binSize = L / nBins;
                densityProfile = new double[nBins];
                virialTensorProfile = new Tensor[nBins];
                for (int i=0; i<nBins; i++) {
                    virialTensorProfile[i] = space.makeTensor();
                }
                uLrc = new double[(nBins+1)/2];
                fLrc = new double[(nBins+1)/2];
                xVirialLrc = new double[(nBins+1)/2];
                yzVirialLrc = new double[(nBins+1)/2];
                double iCutoff2 = 1.0/(cutoff*cutoff);
                double iCutoff4 = iCutoff2*iCutoff2;
                double iCutoff6 = iCutoff4*iCutoff2;
                double iCutoff10 = iCutoff6*iCutoff4;
                double iCutoff12 = iCutoff6*iCutoff6;
                for (int i=0; i<uLrc.length; i++) {
                    double dx = binSize * i;
                    if (dx < cutoff) {
                        uLrc[i] = 4.0*Math.PI*(0.2*iCutoff10 - 0.5*iCutoff4);
                        fLrc[i] = 8.0*Math.PI*(iCutoff12 - iCutoff6);
                        xVirialLrc[i] = 4.0*Math.PI*dx*dx*(iCutoff12 - iCutoff6);
                        yzVirialLrc[i] = 2.0*Math.PI*((1.2*cutoff*cutoff-dx*dx)*iCutoff12 - (1.5*cutoff*cutoff-dx*dx)*iCutoff6);
                    }
                    else {
                        uLrc[i] = 4.0*Math.PI*(0.2/Math.pow(dx, 10.0) - 0.5/Math.pow(dx, 4.0));
                        fLrc[i] = 8.0*Math.PI*(1.0/Math.pow(dx, 12.0) - 1.0/Math.pow(dx, 6.0));
                        xVirialLrc[i] = 4.0*Math.PI*dx*dx*(1.0/Math.pow(dx, 12) - 1.0/Math.pow(dx, 6.0));
                        yzVirialLrc[i] = 2.0*Math.PI*(0.2/Math.pow(dx, 10) - 0.5/Math.pow(dx, 4));
                    }
                }
            }
            else {
                for (int i=0; i<nBins; i++) {
                    densityProfile[i] = 0;
                    virialTensorProfile[i].E(0);
                }
            }
            IAtomSet leafAtoms = box.getLeafList();

            double dV = binSize;
            for (int i=1; i<dim.getD(); i++) {
                dV *= dim.x(i);
            }
            double dVi = 1.0 / dV;
            for (int i=0; i<leafAtoms.getAtomCount(); i++) {
                double x = ((IAtomPositioned)leafAtoms.getAtom(i)).getPosition().x(0);
                // wrap around PBC
                x -= Math.round(x*Li) * L;
                int iBin = (int) ((x + halfL) / binSize);
                densityProfile[iBin] += dVi;
            }
        }
    }
    
    public double[] getDensityProfile() {
        return densityProfile;
    }

    /**
	 * Adds forces due to given potential acting on the atoms produced by the iterator.
	 * Implemented for only 1- and 2-body potentials.
	 */
	public void doCalculation(IAtomSet atoms, IPotential potential) {
		PotentialSoft potentialSoft = (PotentialSoft)potential;
		int nBody = potential.nBody();
		virialTensor.E(0);
		IVector[] f = potentialSoft.gradient(atoms, virialTensor);
		virialTensorTot.PE(virialTensor);
		virialTensor.TE(1.0/binSize);
		boolean takeData = dataCountDown == 0;
		switch(nBody) {
		    case 0:
		        break;
			case 1:
			    // this is our LRC potential
			    ((IntegratorBox.Forcible)integratorAgentManager.getAgent(atoms.getAtom(0))).force().ME(f[0]);
                if (takeData) {
                    double x0 = ((IAtomPositioned)atoms.getAtom(0)).getPosition().x(0);
                    // wrap around PBC
                    x0 -= Math.round(x0*Li) * L;
                    int iBin0 = (int) ((x0 + halfL) / binSize);
                    virialTensorProfile[iBin0].PE(virialTensor);
                }
			    break;
			case 2:
                ((IntegratorBox.Forcible)integratorAgentManager.getAgent(atoms.getAtom(0))).force().ME(f[0]);
                ((IntegratorBox.Forcible)integratorAgentManager.getAgent(atoms.getAtom(1))).force().ME(f[1]);
                if (takeData) {
                    double x0 = ((IAtomPositioned)atoms.getAtom(0)).getPosition().x(0);
                    // wrap around PBC
                    x0 -= Math.round(x0*Li) * L;
                    int iBin0 = (int) ((x0 + halfL) / binSize);
                    double x1 = ((IAtomPositioned)atoms.getAtom(1)).getPosition().x(0);
                    // wrap around PBC
                    x1 -= Math.round(x1*Li) * L;
                    int iBin1 = (int) ((x1 + halfL) / binSize);
                    if (iBin0 == iBin1) {
                        virialTensorProfile[iBin0].PE(virialTensor);
                    }
                    else if (Math.abs(iBin1-iBin0) > nBins/2) {
                        if (iBin0 > nBins/2) {
                            //swap so iBin1 is on the right side, iBin0 on the left
                            int foo = iBin0;
                            iBin0 = iBin1;
                            iBin1 = foo;
                        }
                        int nBetween = (nBins-1)-iBin1+1 + iBin0+1;
                        virialTensor.TE(1.0/nBetween);
                        for (int i=iBin1; i<nBins; i++) {
                            virialTensorProfile[i].PE(virialTensor);
                        }
                        for (int i=0; i<iBin0+1; i++) {
                            virialTensorProfile[i].PE(virialTensor);
                        }
                    }
                    else {
                        if (iBin0 > iBin1) {
                            //swap so iBin1 is on the right side, iBin0 on the left
                            int foo = iBin0;
                            iBin0 = iBin1;
                            iBin1 = foo;
                        }
                        int nBetween = iBin1 - iBin0 + 1;
                        virialTensor.TE(1.0/nBetween);
                        for (int i=iBin0; i<iBin1+1; i++) {
                            virialTensorProfile[i].PE(virialTensor);
                        }
                    }
                }
		 		break;
            default:
                throw new RuntimeException("I didn't expect a sort of Spanish inquisition");
		}
	}

    /**
     * Returns the virial portion of pressure tensor calculated during the last
     * potential calculation.  In order to be valid, reset() must be called
     * before invoking potentialMaster.calculate.  The given tensor has not
     * been normalized by the system volume.
     */
    public Tensor getPressureTensor() {
        return virialTensorTot;
    }
    
    public Tensor[] getPressureTensorProfile() {
        return virialTensorProfile;
    }

    private static final long serialVersionUID = 1L;
    protected final ISpace space;
    protected final Tensor virialTensor, virialTensorTot;
    protected Tensor[] virialTensorProfile;
    protected double binSize;
    protected IBox box;
    protected double[] densityProfile;
    protected double[] uLrc, fLrc, xVirialLrc, yzVirialLrc;
    protected double L, halfL, Li;
    protected int nBins;
    protected int dataCountDown;
    protected double cutoff;
}
