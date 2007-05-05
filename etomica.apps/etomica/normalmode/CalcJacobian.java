package etomica.normalmode;

import etomica.atom.AtomAgentManager;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.normalmode.CoordinateDefinition.SiteSource;
import etomica.phase.Phase;
import etomica.space.IVector;

/**
 * Class that calculates the dq/dx Jacobian.
 * @author Andrew Schultz
 */
public class CalcJacobian {

    public CalcJacobian(int coordinateDim) {
        this.coordinateDim = coordinateDim;
        iterator = new AtomIteratorAllMolecules();
    }
    
    public double[][] getJacobian() {
        int l = coordinateDim*iterator.size();
        double[][] tempJacobian;
        boolean doFull = false;
        if (doFull) {
            tempJacobian = new double[l][l];
        }
        else {
            tempJacobian = new double[l-coordinateDim][l];
        }
        
        int vectorPos = 0;
        if (doFull) {
            for (int j=0; j<coordinateDim; j++) {
                for (int i=0; i<l/coordinateDim; i++) {
                    tempJacobian[j][i*coordinateDim+j] = 1;
                }
            }
            vectorPos = 1;
        }
        for (int iVector = 0; iVector < waveVectors.length; iVector++) {
            iterator.reset();
            // sum T over atoms
            int atomCount = 0;
            double phaseAngle = Double.NaN;
            for (IAtom atom = iterator.nextAtom(); atom != null;
                 atom = iterator.nextAtom()) {
                IVector latticePosition = (IVector)siteManager.getAgent(atom);
                double kR = waveVectors[iVector].dot(latticePosition);
                double coskR = Math.cos(kR);
                double sinkR = Math.sin(kR);
                for (int iDim = 0; iDim < coordinateDim; iDim++) {
                    if (waveVectorCoefficients[iVector] == 1) {
                        tempJacobian[vectorPos*coordinateDim+iDim][atomCount*coordinateDim+iDim] = coskR;
                        tempJacobian[(vectorPos+1)*coordinateDim+iDim][atomCount*coordinateDim+iDim] = -sinkR;
                    }
                    else {
                        // single degree of freedom.
                        // All kR must be 100% in-phase or 100% out of phase. kR-phaseAngle will be 0 or pi
                        if (!Double.isNaN(phaseAngle) && Math.abs(Math.sin(kR-phaseAngle)) > 1.e-10) {
                            throw new RuntimeException("oops "+iVector+" "+waveVectors[iVector]+" "+phaseAngle+" "+kR);
                        }
                        phaseAngle = kR;
                        // either one works so long as it's not 0
                        if (Math.abs(coskR) > 0.1) {
                            tempJacobian[vectorPos*coordinateDim+iDim][atomCount*coordinateDim+iDim] = coskR > 0 ? 1 : -1;
                        }
                        else {
                            tempJacobian[vectorPos*coordinateDim+iDim][atomCount*coordinateDim+iDim] = sinkR > 0 ? 1 : -1;
                        }
                    }
                }
                atomCount++;
            }
            for (int iDim=0; iDim<coordinateDim; iDim++) {
                // zero out the last matrix element if it's statistically 0
                if (Math.abs(tempJacobian[vectorPos*coordinateDim+iDim][l-coordinateDim+iDim]) < 1.e-14) {
                    tempJacobian[vectorPos*coordinateDim+iDim][l-coordinateDim+iDim] = 0;
                }
                if (waveVectorCoefficients[iVector] == 1) {
                    if (Math.abs(tempJacobian[(vectorPos+1)*coordinateDim+iDim][l-coordinateDim+iDim]) < 1.e-14) {
                        tempJacobian[(vectorPos+1)*coordinateDim+iDim][l-coordinateDim+iDim] = 0;
                    }
                }
                for (int j = 0; j < l - coordinateDim; j+=coordinateDim) {
                    // subtract the element corresonding to the last atom from the other atoms
                    if (!doFull) {
                        tempJacobian[vectorPos*coordinateDim+iDim][j+iDim] -= tempJacobian[vectorPos*coordinateDim+iDim][l-coordinateDim+iDim];
                    }
                    if (Math.abs(tempJacobian[vectorPos*coordinateDim+iDim][j+iDim]) < 1.e-14) {
                        // zero out the element if it's statistically 0
                        tempJacobian[vectorPos*coordinateDim+iDim][j+iDim] = 0;
                    }
                    if (waveVectorCoefficients[iVector] == 1) {
                        // if both real and imaginary are important, handle the imaginary part
                        if (!doFull) {
                            tempJacobian[(vectorPos+1)*coordinateDim+iDim][j+iDim] -= tempJacobian[vectorPos*coordinateDim+iDim+1][l-coordinateDim+iDim];
                        }
                        if (Math.abs(tempJacobian[(vectorPos+1)*coordinateDim+iDim][j+iDim]) < 1.e-14) {
                            tempJacobian[(vectorPos+1)*coordinateDim+iDim][j+iDim] = 0;
                        }
                    }
                }
            }
            vectorPos += Math.round(waveVectorCoefficients[iVector]*2);
        }
        if (doFull) {
            return tempJacobian;
        }
        double[][] jacobian = new double[l-coordinateDim][l-coordinateDim];
        for (int i=0; i<l-coordinateDim; i++) {
            for (int j=0; j<l-coordinateDim; j++) {
                jacobian[i][j] = tempJacobian[i][j];
            }
        }
        return jacobian;
    }

    /**
     * Sets the object that defines the normal-coordinate wave vectors.
     */
    public void setWaveVectorFactory(WaveVectorFactory newWaveVectorFactory) {
        waveVectorFactory = newWaveVectorFactory;
    }

    /**
     * @return the WaveVectorFactory last given via the set methods.
     */
    public WaveVectorFactory getWaveVectorFactory() {
        return waveVectorFactory;
    }
    
    public void setPhase(Phase newPhase) {

        realT = new double[coordinateDim];
        imaginaryT = new double[coordinateDim];

        waveVectorFactory.makeWaveVectors(newPhase);
        waveVectors = waveVectorFactory.getWaveVectors();
        waveVectorCoefficients = waveVectorFactory.getCoefficients();

        siteManager = new AtomAgentManager(new SiteSource(newPhase.getSpace()), newPhase);
        iterator.setPhase(newPhase);
        iterator.reset();
    }
    
    protected void setWaveVectors(IVector[] newWaveVectors, double[] coefficients) {
        waveVectors = newWaveVectors;
        waveVectorCoefficients = coefficients;
    }
    
    private static final long serialVersionUID = 1L;
    protected int coordinateDim;
    protected double[] realT, imaginaryT;
    protected IVector[] waveVectors;
    protected double[] waveVectorCoefficients;
    protected NormalModes normalModes;
    private final AtomIteratorAllMolecules iterator;
    private AtomAgentManager siteManager;
    private WaveVectorFactory waveVectorFactory;
}
