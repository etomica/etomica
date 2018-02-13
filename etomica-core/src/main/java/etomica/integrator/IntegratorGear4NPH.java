/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// includes a main method

package etomica.integrator;

import etomica.action.BoxInflate;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.meter.MeterTemperature;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.modifier.ModifierBoolean;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Kelvin;
import etomica.units.dimensions.*;
import etomica.util.random.IRandom;

/**
 * Gear 4th-order predictor-corrector integrator for constant enthalphy, pressure.
 *
 * @author Ed Maginn
 * @author David Kofke
 */
public class IntegratorGear4NPH extends IntegratorGear4 {

    private static final long serialVersionUID = 1L;
    double vol1, vol2, vol3, vol4;
    protected /*final*/ ForceSumNPH forceSumNPH;//MeterTPH won't permit this to be final (?)
    private final IteratorDirective allAtoms = new IteratorDirective();
    protected final BoxInflate inflate;
    double targetH;
    double targetP;
    double targetT = Kelvin.UNIT.toSim(300.);
    private double rrp = 300.;
    private double rrh = 300.;
    private double kp, kh;
    protected int D;
    protected MeterTemperature meterTemperature;
    
    public IntegratorGear4NPH(Simulation sim, PotentialMaster potentialMaster, Space _space, Box box) {
        this(potentialMaster, sim.getRandom(),0.05, 1.0, _space, box);
    }
    
    public IntegratorGear4NPH(PotentialMaster potentialMaster, IRandom random,
                              double timeStep, double temperature, Space _space, Box box) {
        super(potentialMaster, random, timeStep, temperature, _space, box);
        kp = 1.0/rrp/getTimeStep();
        kh = 1.0/rrh/getTimeStep();
        D = space.D();
        setIsothermal(true);
        forceSumNPH = new ForceSumNPH(space);
        inflate = new BoxInflate(space);
    }

    public void setTimeStep(double t) {
        super.setTimeStep(t);
        kp = 1.0/(rrp*t);
        kh = 1.0/(rrh*t);
    }
    
    public void setRelaxationRateP(double value) {
        rrp = value;
        setTimeStep(timeStep);
    }
    public double getRelaxationRateP() {return rrp;}
    public Dimension getRelaxationRatePDimension() {return Null.DIMENSION;}
    
    public void setRelaxationRateH(double value) {
        rrh = value;
        setTimeStep(timeStep);
    }
    public double getRelaxationRateH() {return rrh;}
    public Dimension getRelaxationRateHDimension() {return Null.DIMENSION;}
    
    public synchronized void setTargetH(double value) {targetH = value;}
    public double getTargetH() {return targetH;}
    public Dimension getTargetHDimension() {return Energy.DIMENSION;}
    
    public synchronized void setTargetP(double value) {targetP = value;}
    public double getTargetP() {return targetP;}
    public Dimension getTargetPDimension() {return Pressure.DIMENSION;}
    
    public synchronized void setTargetT(double value) {targetT = value;}
    public double getTargetT() {return targetT;}
    public Dimension getTargetTDimension() {return Temperature.DIMENSION;}

    public void setBox(Box box) {
        super.setBox(box);
        inflate.setBox(this.box);
        meterTemperature = new MeterTemperature(this.box, D);
        forceSumNPH.setBox(this.box);
        forceSumNPH.setAgentManager(forces);
    }
    
    
//--------------------------------------------------------------
// steps all particles across time interval tStep

    protected void doStepInternal() {
        super.doStepInternal();
        predictor();
        calculateForces();
        if(isothermal) drivePT();
        else drivePH();
        corrector();
    }//end of doStep
    
    public void drivePT() {
        double kineticT = meterTemperature.getDataAsScalar();
        double mvsq = kineticT * D * box.getLeafList().getAtomCount();
        double volume = box.getBoundary().volume();
        double density = box.getMoleculeList().getMoleculeCount() / box.getBoundary().volume();
        double pCurrent = density*kineticT - forceSumNPH.w/(D*volume);
        double pDot = kp*(targetP - pCurrent);
        double kDot = kp*(targetT - kineticT)*box.getLeafList().getAtomCount();
        chi = ( - forceSumNPH.rvx - D*pDot*volume)/
                    ( forceSumNPH.x + D*D*pCurrent*volume);
        zeta = (forceSumNPH.vf - kDot)/mvsq - chi;
    }
    
    public void drivePH() {
        double kineticT = meterTemperature.getDataAsScalar();
        double mvsq = kineticT * D * box.getLeafList().getAtomCount();
        double volume = box.getBoundary().volume();
        double density = box.getMoleculeList().getMoleculeCount() / box.getBoundary().volume();
        double pCurrent = density*kineticT - forceSumNPH.w/(D*volume);
        double hCurrent = 0.5*mvsq + forceSumNPH.u + pCurrent*volume;
        double pDot = kp*(targetP - pCurrent);
        double hDot = kh*(targetH - hCurrent);
        zeta = (pDot*volume - hDot)/mvsq;
        chi = (2.0*forceSumNPH.vf -forceSumNPH.rvx - D*pDot*volume - 2*zeta*mvsq)/
                    (2.0*mvsq + forceSumNPH.x + D*D*pCurrent*volume);
    }
    
    protected void calculateForces() {
        //Compute all forces

        //zero forces on all atoms
        forceSumNPH.reset();

        forceSumNPH.u = 0.0;
        forceSumNPH.w = 0.0;
        forceSumNPH.x = 0.0;
        forceSumNPH.rvx = 0.0;
        forceSumNPH.vf = 0.0;

        potentialMaster.calculate(box, allAtoms, forceSumNPH);
    }//end of calculateForces
    
    protected void corrector() {
        super.corrector();
        double volOld = box.getBoundary().volume();
        // voi (dV/dt... voi???) = Dr^(D-1) dr/dt, chi = dr/dt
        double voi = D*volOld*chi/box.getBoundary().getBoxSize().getX(0);
        double corvol = voi - vol1;
        double volNew = volOld + c0*corvol;
        vol1 = voi;
        vol2 += c2*corvol;
        vol3 += c3*corvol;
        vol4 += c4*corvol;
        double rScale = Math.pow(volNew/volOld,1.0/D);
        inflate.setScale(rScale);
        inflate.actionPerformed();
    }//end of corrector
        
    protected void predictor() {
        super.predictor();
        double volOld = box.getBoundary().volume();
        double volNew = volOld + p1*vol1 + p2*vol2 + p3*vol3 + p4*vol4;
        vol1 += p1*vol2 + p2*vol3 + p3*vol4;
        vol2 += p1*vol3 + p2*vol4;
        vol3 += p1*vol4;
        double rScale = Math.pow(volNew/volOld,1.0/D);
        inflate.setScale(rScale);
        inflate.actionPerformed();
    }
//--------------------------------------------------------------


    public void reset() {
        super.reset();
        vol1 = 0.0;
        vol2 = 0.0;
        vol3 = 0.0;
        vol4 = 0.0;
    }
              
    //inner class used to toggle between NPT and NPH ensembles
    public static class EnsembleToggler implements ModifierBoolean, java.io.Serializable {
        public EnsembleToggler(IntegratorGear4NPH integrator, int D) {
            this.integrator = integrator;
            this.dim = D;
        }
        public void setBoolean(boolean isothermal) {
            integrator.setIsothermal(isothermal);
            if(!isothermal) {
                integrator.calculateForces();
                Box box = integrator.getBox();
                double kineticT = integrator.getMeterTemperature().getDataAsScalar();
                double mvsq = kineticT * dim * box.getLeafList().getAtomCount();
                double volume = box.getBoundary().volume();
                double density = box.getMoleculeList().getMoleculeCount() / box.getBoundary().volume();
                double pCurrent = density*kineticT - integrator.forceSumNPH.w/(dim*volume);
                double hCurrent = 0.5*mvsq + integrator.forceSumNPH.u + pCurrent*volume;
                integrator.setTargetH(hCurrent);
            }
        }
        private static final long serialVersionUID = 1L;
        public boolean getBoolean() {return integrator.isIsothermal();}
        private IntegratorGear4NPH integrator;
        private final int dim;
    }

    //meter for enthalpy, obtaining values from
    //most recent call to the ForceSumNPH instance
    public static final class MeterEnthalpy implements IDataSource, java.io.Serializable {
        
        public MeterEnthalpy(IntegratorGear4NPH integrator, int D) {
            data = new DataDouble();
            dataInfo = new DataInfoDouble("Enthalpy", Energy.DIMENSION);
            this.integrator = integrator;
            tag = new DataTag();
            dataInfo.addTag(tag);
            dim = D;
        }
        
        public IDataInfo getDataInfo() {
            return dataInfo;
        }
        
        public DataTag getTag() {
            return tag;
        }
        
        public IData getData() {
            Box box = integrator.getBox();
            double kineticT = integrator.getMeterTemperature().getDataAsScalar();
            double mvsq = kineticT* dim * box.getLeafList().getAtomCount();
            double volume = box.getBoundary().volume();
            double density = box.getMoleculeList().getMoleculeCount() / box.getBoundary().volume();
            double pCurrent = density*kineticT - integrator.forceSumNPH.w/(dim*volume);
            double hCurrent = 0.5*mvsq + integrator.forceSumNPH.u + pCurrent*volume;
            data.x = hCurrent/box.getMoleculeList().getMoleculeCount();
            return data;
        }
        
        private static final long serialVersionUID = 1L;
        private DataDouble data;
        private IntegratorGear4NPH integrator;
        private IDataInfo dataInfo;
        private DataTag tag;
        private final int dim;
    }
        
    public static final class ForceSumNPH extends PotentialCalculationForceSum {
        
        private static final long serialVersionUID = 1L;
        double u; //energy sum
        double w; //virial sum
        double x; //hypervirial sum
        double rvx; 
        double vf;
        private final Vector dr;
        private final Vector dv;
        private Boundary nearestImageTransformer;

        public ForceSumNPH(Space space) {
            dr = space.makeVector();
            dv = space.makeVector();
        }
        
        public void setBox(Box box) {
            nearestImageTransformer = box.getBoundary();
        }
        
        //pair
        public void doCalculation(IAtomList pair, IPotential potential2) {
            Potential2Soft potentialSoft = (Potential2Soft)potential2;
            IAtomKinetic atom0 = (IAtomKinetic)pair.getAtom(0);
            IAtomKinetic atom1 = (IAtomKinetic)pair.getAtom(1);
            dv.Ev1Mv2(atom1.getVelocity(), atom0.getVelocity());
            
            dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
            nearestImageTransformer.nearestImage(dr);

            double r2 = dr.squared();
            u += potentialSoft.energy(pair);
            w += potentialSoft.virial(pair);
            double hv = potentialSoft.hyperVirial(pair);
            x += hv;
            rvx += hv * dr.dot(dv)/r2;
            Vector[] f = potentialSoft.gradient(pair);
            vf += dv.dot(f[0]); //maybe should be (-)?
            integratorAgentManager.getAgent((IAtom)atom0).ME(f[0]);
            integratorAgentManager.getAgent((IAtom)atom1).ME(f[1]);
        }
    }
}
