/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.action.AtomActionTranslateBy;
import etomica.action.BoxImposePbc;
import etomica.action.MoleculeChildAtomAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.*;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationLattice;
import etomica.data.meter.MeterKineticEnergyRigid;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet.MyAgent;
import etomica.lattice.LatticeCubicFcc;
import etomica.listener.IntegratorListenerAction;
import etomica.math.function.Function;
import etomica.models.water.DipoleSourceWater;
import etomica.models.water.OrientationCalcWater3P;
import etomica.models.water.P2WaterSPCSoft;
import etomica.models.water.SpeciesWater3POriented;
import etomica.molecule.*;
import etomica.molecule.MoleculeAgentManager.MoleculeAgentSource;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.*;
import etomica.space3d.IOrientationFull3D;
import etomica.space3d.OrientationFull3D;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.ISpeciesOriented;
import etomica.units.Electron;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.Debug;

import java.awt.*;
import java.io.Serializable;

public class IntegratorRigidIterative extends IntegratorMD implements AgentSource<IntegratorVelocityVerlet.MyAgent>, SpeciesAgentManager.AgentSource, MoleculeAgentSource {

    private static final long serialVersionUID = 2L;
    protected PotentialCalculationTorqueSum torqueSum;
    private final IteratorDirective allAtoms;
    protected final Tensor pressureTensor;
    protected final Tensor workTensor;
    protected final RotationTensor3D rotationTensor;
    protected final RotationTensor3D tempRotationTensor;
    protected final Vector xWork;
    protected final SpeciesAgentManager typeAgentManager;
    protected final Vector tempAngularVelocity;
    protected final MoleculePositionCOM atomPositionCOM;
    protected final AtomActionTranslateBy translateBy;
    protected final MoleculeChildAtomAction translator;
    protected final OrientationFull3D tempOrientation;
    public int printInterval;
    protected int maxIterations;
    protected double omegaTolerance;
    protected final boolean storeAngularMomentum = false;

    protected final Simulation sim;
    protected AtomLeafAgentManager<IntegratorVelocityVerlet.MyAgent> leafAgentManager;
    protected MoleculeAgentManager moleculeAgentManager;

    public IntegratorRigidIterative(Simulation sim, PotentialMaster potentialMaster, Space _space) {
        this(sim, potentialMaster, 0.05, 1.0, _space);
    }
    
    public IntegratorRigidIterative(Simulation sim, PotentialMaster potentialMaster,
                                    double timeStep, double temperature, Space _space) {
        super(potentialMaster,sim.getRandom(),timeStep,temperature, _space);
        this.sim = sim;
        // if you're motivated to throw away information earlier, you can use 
        // PotentialCalculationForceSum instead.
        torqueSum = new PotentialCalculationTorqueSum();
        allAtoms = new IteratorDirective();
        // allAtoms is used only for the force calculation, which has no LRC
        // but we're also calculating the pressure tensor, which does have LRC.
        // things deal with this OK.
        allAtoms.setIncludeLrc(true);
        pressureTensor = space.makeTensor();
        workTensor = space.makeTensor();
        rotationTensor = (RotationTensor3D)space.makeRotationTensor();
        tempRotationTensor = (RotationTensor3D)space.makeRotationTensor();
        xWork = space.makeVector();
        typeAgentManager = new SpeciesAgentManager(this, sim);
        tempAngularVelocity = space.makeVector();
        tempOrientation = new OrientationFull3D(space);
        atomPositionCOM = new MoleculePositionCOM(space);
        translateBy = new AtomActionTranslateBy(space);
        translator = new MoleculeChildAtomAction(translateBy);
        maxIterations = 20;
        omegaTolerance = 1.e-30;
        meterKE = new MeterKineticEnergyRigid(space, sim);
    }

    public void setBox(Box box) {
        if (this.box != null) {
            // allow agentManager to de-register itself as a BoxListener
            leafAgentManager.dispose();
            moleculeAgentManager.dispose();
            leafAgentManager = null;
            moleculeAgentManager = null;
        }
        super.setBox(box);
        leafAgentManager = new AtomLeafAgentManager<IntegratorVelocityVerlet.MyAgent>(this, box,IntegratorVelocityVerlet.MyAgent.class);
        moleculeAgentManager = new MoleculeAgentManager(sim, box, this);
        torqueSum.setAgentManager(leafAgentManager);
        torqueSum.setMoleculeAgentManager(moleculeAgentManager);
        ((MeterKineticEnergyRigid)meterKE).setBox(box);
    }

    public void setTorqueSum(PotentialCalculationTorqueSum pc) {
        torqueSum = pc;
        if (leafAgentManager != null) {
            torqueSum.setAgentManager(leafAgentManager);
            torqueSum.setMoleculeAgentManager(moleculeAgentManager);
        }
    }

    public void setOrientationCalc(ISpecies moleculeType, OrientationCalc calcer) {
        typeAgentManager.setAgent(moleculeType, calcer);
    }
    
    public void setMaxIterations(int newMaxIterations) {
        maxIterations = newMaxIterations;
    }
    
//--------------------------------------------------------------
// steps all particles across time interval tStep

    // assumes one box
    public void doStepInternal() {
        super.doStepInternal();
        currentKineticEnergy = 0;
        int iterationsTotal = 0;
        int numRigid = 0;
        if (Debug.ON && Debug.DEBUG_NOW) {
            IAtomList pair = Debug.getAtoms(box);
            if (pair != null) {
                Vector dr = space.makeVector();
                dr.Ev1Mv2(pair.getAtom(1).getPosition(), pair.getAtom(0).getPosition());
                System.out.println(pair+" dr "+dr);
            }
        }
        IMoleculeList moleculeList = box.getMoleculeList();
        int nMolecules = moleculeList.getMoleculeCount();
        for (int iMolecule = 0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = moleculeList.getMolecule(iMolecule);
            IAtomList children = molecule.getChildList();
            OrientationCalc calcer = (OrientationCalc)typeAgentManager.getAgent(molecule.getType());
            if (calcer == null) {
                for (int iLeaf=0; iLeaf<children.getAtomCount(); iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic)children.getAtom(iLeaf);
                    MyAgent agent = leafAgentManager.getAgent(a);
                    Vector r = a.getPosition();
                    Vector v = a.getVelocity();
                    if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                        System.out.println("first "+a+" r="+r+", v="+v+", f="+agent.force);
                    }
                    v.PEa1Tv1(0.5*timeStep* a.getType().rm(),agent.force);  // p += f(old)*dt/2
                    r.PEa1Tv1(timeStep,v);         // r += p*dt/m
                }
                continue;
            }
            
            MoleculeAgent agent = (MoleculeAgent)moleculeAgentManager.getAgent(molecule);
            IMoleculeOrientedKinetic orientedMolecule = (IMoleculeOrientedKinetic)molecule;
            IOrientationFull3D orientation = (IOrientationFull3D)orientedMolecule.getOrientation();
            Vector moment = ((ISpeciesOriented)molecule.getType()).getMomentOfInertia();

            // use the angular velocity field to store angular momentum during the time step  :(
            Vector angularMomentum = orientedMolecule.getAngularVelocity();
//            System.out.println("angular velocity "+angularMomentum);
            // transform to body-fixed, multiply by moment of inertia, transform back
            tempAngularVelocity.E(orientedMolecule.getAngularVelocity());
            if (!storeAngularMomentum) {
                rotationTensor.setOrientation(orientation);
                rotationTensor.transform(angularMomentum);
                angularMomentum.TE(moment);
                rotationTensor.invert();
                rotationTensor.transform(angularMomentum);
            }
            else {
                rotationTensor.setOrientation(orientation);
                rotationTensor.transform(tempAngularVelocity);
                tempAngularVelocity.DE(moment);
                rotationTensor.invert();
                rotationTensor.transform(tempAngularVelocity);
            }
            
            // xWork = angular velocity from the previous iteration
            // save angular velocity to see if we've converged
            xWork.E(tempAngularVelocity);
//            System.out.println("initial guess half-timestep angular velocity "+angularVelocity);
            angularMomentum.PEa1Tv1(0.5*timeStep, agent.torque);

            double dtheta = Math.sqrt(tempAngularVelocity.squared());
            if (dtheta > 0) {
                tempAngularVelocity.TE(1/dtheta);
            }

            numRigid++;
            for (int i = 0; i<maxIterations; i++) {
                iterationsTotal++;
                tempOrientation.E(orientation);
                if (dtheta > 0) {
                    // estimate orientation at half timestep, we need this to get the moment of inertia at half timestep
                    tempOrientation.rotateBy(0.5*dtheta*timeStep, tempAngularVelocity);
                }

                tempAngularVelocity.E(angularMomentum);

                rotationTensor.setOrientation(tempOrientation);

                //find body-fixed angular momentum
                rotationTensor.transform(tempAngularVelocity);
                //now divide out moment of inertia to get body-fixed angular velocity
                tempAngularVelocity.DE(moment);
                rotationTensor.invert();
                //now rotate back to get space-fixed angular velocity
                rotationTensor.transform(tempAngularVelocity);
//                System.out.println("iteration "+i+" half-timestep angular velocity "+angularVelocity);
                
                xWork.ME(tempAngularVelocity);
                double omegaError = xWork.squared() / tempAngularVelocity.squared();
                xWork.E(tempAngularVelocity);
                dtheta = Math.sqrt(tempAngularVelocity.squared());
                if (dtheta > 0) {
                    tempAngularVelocity.TE(1.0/dtheta);
                }

                if (omegaError < omegaTolerance) {
                    break;
                }
                if (i == maxIterations-1) {
                    System.err.println("omegaError still "+omegaError+" after "+maxIterations+" iterations");
//                    throw new RuntimeException("omegaError still "+omegaError+" after "+maxIterations+" iterations");
                }
            }

            if (dtheta > 0) {
                orientation.rotateBy(timeStep*dtheta, tempAngularVelocity);
            }
                
//            System.out.println("o "+orientation.getDirection()+" "+orientation.getSecondaryDirection());
            calcer.setOrientation(molecule, orientation);
            //advance linear velocity to half timestep
            orientedMolecule.getVelocity().PEa1Tv1(0.5*timeStep/((ISpeciesOriented)molecule.getType()).getMass(), agent.force);
            
            //advance position to full timestep
            Vector transVec = ((AtomActionTranslateBy)translator.getAtomAction()).getTranslationVector();
            transVec.Ea1Tv1(timeStep, orientedMolecule.getVelocity());
            orientedMolecule.getPosition().PE(transVec);
            translator.actionPerformed(molecule);
        }

//        if(forceSum instanceof PotentialCalculationForcePressureSum){
//            pressureTensor.E(((PotentialCalculationForcePressureSum)forceSum).getPressureTensor());
//        }
        
        torqueSum.reset();
        //Compute forces and torques on each atom at t+dt
        potentialMaster.calculate(box, allAtoms, torqueSum);
        
        for (int iMolecule = 0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = moleculeList.getMolecule(iMolecule);
            IAtomList children = molecule.getChildList();
            if (typeAgentManager.getAgent(molecule.getType()) == null) {
                // unimolecular or at least not rigid
                //Finish integration step
                for (int iLeaf=0; iLeaf<children.getAtomCount(); iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic)children.getAtom(iLeaf);
//                    System.out.println("force: "+((MyAgent)a.ia).force.toString());
                    Vector velocity = a.getVelocity();
                    workTensor.Ev1v2(velocity,velocity);
                    workTensor.TE(a.getType().getMass());
                    pressureTensor.PE(workTensor);
                    if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                        System.out.println("second "+a+" v="+velocity+", f="+leafAgentManager.getAgent(a).force);
                    }
                    velocity.PEa1Tv1(0.5*timeStep* a.getType().rm(),leafAgentManager.getAgent(a).force);  //p += f(new)*dt/2
                    currentKineticEnergy += velocity.squared()* a.getType().getMass();
                }
                // skip the rotational stuff
                continue;
            }
            
            IMoleculeOrientedKinetic orientedMolecule = (IMoleculeOrientedKinetic)molecule;
            MoleculeAgent agent = (MoleculeAgent)moleculeAgentManager.getAgent(molecule);
            //calc torque and linear force
            for (int i=0; i<children.getAtomCount(); i++) {
                IAtom atom = children.getAtom(i);
                Vector atomForce = leafAgentManager.getAgent(atom).force;
                if (atomForce.isZero()) {
                    continue;
                }

                agent.force.PE(atomForce);

                xWork.Ev1Mv2(atom.getPosition(), orientedMolecule.getPosition());
                xWork.XE(atomForce);
                agent.torque.PE(xWork);
            }
            
            //advance linear velocity to full timestep
            double mass = ((ISpeciesOriented)molecule.getType()).getMass();
            orientedMolecule.getVelocity().PEa1Tv1(0.5*timeStep/mass, agent.force);

            //advance momentum to full timestep
            Vector angularVelocity = orientedMolecule.getAngularVelocity();
            Vector moment = ((ISpeciesOriented)molecule.getType()).getMomentOfInertia();

            // we actually stored the half-timestep angular momentum in this field...
            // advance to full timestep
            angularVelocity.PEa1Tv1(0.5*timeStep, agent.torque);
            tempAngularVelocity.E(angularVelocity);
            // we actually had this before, but we didn't save it...
            // transform to body-fixed, multiply by moment of inertia, transform back
            // ==> angular momentum
            rotationTensor.setOrientation((IOrientationFull3D)orientedMolecule.getOrientation());
            rotationTensor.transform(tempAngularVelocity);
            tempAngularVelocity.DE(moment);
            xWork.E(tempAngularVelocity);
            if (!storeAngularMomentum) {
                rotationTensor.invert();
                rotationTensor.transform(tempAngularVelocity);
                angularVelocity.E(tempAngularVelocity);
            }

            currentKineticEnergy += mass * orientedMolecule.getVelocity().squared();
            xWork.TE(xWork);
            currentKineticEnergy += xWork.dot(moment);
        }
        pressureTensor.TE(1/box.getBoundary().volume());

        currentKineticEnergy *= 0.5;
        if (printInterval > 0 && stepCount%printInterval == 0) {
            double PE = meterPE.getDataAsScalar();
            int moleculeCount = box.getMoleculeList().getMoleculeCount();
            double fac = Joule.UNIT.fromSim(1.0/moleculeCount)*Constants.AVOGADRO;
            System.out.println(currentTime+" "+(iterationsTotal/(double)numRigid)+" "+Kelvin.UNIT.fromSim(currentKineticEnergy/moleculeCount/3)+" "
                              +fac*currentKineticEnergy+" "+fac*PE+" "+fac*(PE+currentKineticEnergy));
        }

        if (isothermal) {
            doThermostatInternal();
        }
    }

    /**
     * Returns the pressure tensor based on the forces calculated during the
     * last time step.
     */
    public Tensor getPressureTensor() {
        return pressureTensor;
    }

    public void scaleMomenta() {
        double KE = 0;
        IMoleculeList moleculeList = box.getMoleculeList();
        int nMolecules = moleculeList.getMoleculeCount();
        int D = 0;
        momentum.E(0);
        double totalMass = 0;
        for (int iMolecule = 0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = moleculeList.getMolecule(iMolecule);
            IAtomList children = molecule.getChildList();
            if (typeAgentManager.getAgent(molecule.getType()) == null) {
                for (int iLeaf=0; iLeaf<children.getAtomCount(); iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic)children.getAtom(iLeaf);
                    double mass = a.getType().getMass();
                    momentum.PEa1Tv1(mass, a.getVelocity());
                    totalMass += mass;
                }
                continue;
            }
            IAtomOrientedKinetic orientedMolecule = (IAtomOrientedKinetic)molecule;
            double mass = ((ISpeciesOriented)((IMolecule)orientedMolecule).getType()).getMass();
            momentum.PEa1Tv1(mass, orientedMolecule.getVelocity());
            totalMass += mass;
        }
        momentum.TE(1.0/totalMass);
        for (int iMolecule = 0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = moleculeList.getMolecule(iMolecule);
            IAtomList children = molecule.getChildList();
            if (typeAgentManager.getAgent(molecule.getType()) == null) {
                // unimolecular or at least not rigid
                //Finish integration step
                for (int iLeaf=0; iLeaf<children.getAtomCount(); iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic)children.getAtom(iLeaf);
//                    System.out.println("force: "+((MyAgent)a.ia).force.toString());
                    Vector velocity = a.getVelocity();
                    velocity.ME(momentum);
                    KE += velocity.squared() * a.getType().getMass();
                    D += 3;
                }
                continue;
            }
            IAtomOrientedKinetic orientedMolecule = (IAtomOrientedKinetic)molecule;
            orientedMolecule.getVelocity().ME(momentum);
            KE += orientedMolecule.getVelocity().squared() * ((ISpeciesOriented)molecule.getType()).getMass();

            Vector moment = ((ISpeciesOriented)molecule.getType()).getMomentOfInertia();

            tempAngularVelocity.E(orientedMolecule.getAngularVelocity());
            rotationTensor.setOrientation((IOrientationFull3D)orientedMolecule.getOrientation());
            rotationTensor.transform(tempAngularVelocity);
            tempAngularVelocity.TE(moment);

            tempAngularVelocity.TE(tempAngularVelocity);
            tempAngularVelocity.DE(moment);
            KE += tempAngularVelocity.getX(0) + tempAngularVelocity.getX(1)+ tempAngularVelocity.getX(2);
            D += 6; //3;
        }
        if (temperature == 0 && KE == 0) {
            return;
        }
        double scale = Math.sqrt(temperature*D/KE);
        currentKineticEnergy = 0.5*KE*scale*scale;
//        System.out.println("initial KE "+currentKineticEnergy);
        for (int iMolecule = 0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = moleculeList.getMolecule(iMolecule);
            IAtomList children = molecule.getChildList();
            if (typeAgentManager.getAgent(molecule.getType()) == null) {
                // unimolecular or at least not rigid
                //Finish integration step
                for (int iLeaf=0; iLeaf<children.getAtomCount(); iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic)children.getAtom(iLeaf);
//                    System.out.println("force: "+((MyAgent)a.ia).force.toString());
                    a.getVelocity().TE(scale);
                }
                continue;
            }
            IAtomOrientedKinetic orientedMolecule = (IAtomOrientedKinetic)molecule;
            orientedMolecule.getVelocity().TE(scale);

            Vector angularVelocity = orientedMolecule.getAngularVelocity();
            rotationTensor.setOrientation((IOrientationFull3D)orientedMolecule.getOrientation());
            rotationTensor.transform(angularVelocity);
            angularVelocity.TE(scale);
            rotationTensor.invert();
            rotationTensor.transform(angularVelocity);
        }
    }

    public void randomizeMomenta() {
        IMoleculeList moleculeList = box.getMoleculeList();
        int nMolecules = moleculeList.getMoleculeCount();
//        System.out.println("rerandomize");
        for (int iMolecule = 0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = moleculeList.getMolecule(iMolecule);
            OrientationCalc calcer = (OrientationCalc)typeAgentManager.getAgent(molecule.getType());
            if (calcer == null) {
                IAtomList children = molecule.getChildList();
                for (int i=0; i<children.getAtomCount(); i++) {
                    super.randomizeMomentum((IAtomKinetic)children.getAtom(i));
                }
                continue;
            }

            MoleculeOrientedDynamic orientedMolecule = (MoleculeOrientedDynamic)molecule;
            Vector velocity = orientedMolecule.getVelocity();
            Vector moment = ((ISpeciesOriented)molecule.getType()).getMomentOfInertia();

            double mass = ((ISpeciesOriented)molecule.getType()).getMass();
            int D = velocity.getD();
            for(int i=0; i<D; i++) {
                velocity.setX(i,random.nextGaussian());
            }
            velocity.TE(Math.sqrt(temperature/mass));
    
            Vector angularVelocity = orientedMolecule.getAngularVelocity();
            for(int i=0; i<D; i++) {
                angularVelocity.setX(i,random.nextGaussian());
            }

            tempAngularVelocity.Ea1Tv1(temperature, moment);
            tempAngularVelocity.map(new Function.Sqrt());
            tempAngularVelocity.TE(angularVelocity);
            angularVelocity.E(tempAngularVelocity);
            //angularVelocity is now the correct body-fixed angular momentum
            if (!storeAngularMomentum) {
                angularVelocity.DE(moment);
            }
            
            calcer.calcOrientation(orientedMolecule, (IOrientationFull3D)orientedMolecule.getOrientation());
            rotationTensor.setOrientation((IOrientationFull3D)orientedMolecule.getOrientation());
            // body-fixed to space-fixed, so invert
            rotationTensor.invert();
            // transform to space-fixed angular momentum
            rotationTensor.transform(angularVelocity);
        }
        
    }
    
    public void randomizeMomentum(IAtomKinetic atom) {
        if (atom instanceof Atom) {
            super.randomizeMomentum(atom);
        }
        else {
//            System.out.println("rerandomize 1");
            OrientationCalc calcer = (OrientationCalc)typeAgentManager.getAgent(((IMolecule)atom).getType());
            if (calcer == null) {
                super.randomizeMomentum(atom);
                return;
            }

            MoleculeOrientedDynamic orientedMolecule = (MoleculeOrientedDynamic)atom;
            Vector velocity = orientedMolecule.getVelocity();
            double mass = ((ISpeciesOriented)orientedMolecule.getType()).getMass();
            int D = velocity.getD();
            for(int i=0; i<D; i++) {
                velocity.setX(i,random.nextGaussian());
            }
            velocity.TE(Math.sqrt(temperature/mass));

            Vector angularVelocity = orientedMolecule.getAngularVelocity();
            Vector moment = ((ISpeciesOriented)orientedMolecule.getType()).getMomentOfInertia();
            for(int i=0; i<D; i++) {
                angularVelocity.setX(i,random.nextGaussian());
            }
            tempAngularVelocity.Ea1Tv1(temperature, moment);
            tempAngularVelocity.map(new Function.Sqrt());
            tempAngularVelocity.TE(angularVelocity);
            angularVelocity.E(tempAngularVelocity);
            //angularVelocity is now the correct body-fixed angular momentum
            angularVelocity.DE(moment);
            
            calcer.calcOrientation(orientedMolecule, (IOrientationFull3D)orientedMolecule.getOrientation());
            rotationTensor.setOrientation((IOrientationFull3D)orientedMolecule.getOrientation());
            // body-fixed to space-fixed, so invert
            rotationTensor.invert();
            // transform to space-fixed angular momentum
            rotationTensor.transform(angularVelocity);
        }
    }
    
    public void reset() {
        super.reset();
        if (Debug.ON && Debug.DEBUG_NOW) {
            IAtomList pair = Debug.getAtoms(box);
            if (pair != null) {
                Vector dr = space.makeVector();
                dr.Ev1Mv2(pair.getAtom(1).getPosition(), pair.getAtom(0).getPosition());
                System.out.println(pair+" dr "+dr);
            }
        }
        
        IMoleculeList moleculeList = box.getMoleculeList();
        int nMolecules = moleculeList.getMoleculeCount();

        for (int iMolecule=0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = moleculeList.getMolecule(iMolecule);
            OrientationCalc calcer = (OrientationCalc)typeAgentManager.getAgent(molecule.getType());
            if (calcer == null) {
                continue;
            }

            MoleculeOrientedDynamic orientedMolecule = (MoleculeOrientedDynamic)molecule;
            calcer.calcOrientation(molecule, (IOrientationFull3D)orientedMolecule.getOrientation());
            orientedMolecule.getPosition().E(atomPositionCOM.position(molecule));
        }

        torqueSum.reset();
        //Compute forces on each atom at t+dt
        potentialMaster.calculate(box, allAtoms, torqueSum);

        for (int iMolecule=0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = moleculeList.getMolecule(iMolecule);
            if (typeAgentManager.getAgent(molecule.getType()) == null) {
                continue;
            }

            MoleculeAgent agent = (MoleculeAgent)moleculeAgentManager.getAgent(molecule);
            IMoleculeOrientedKinetic orientedMolecule = (IMoleculeOrientedKinetic)molecule;
            //calc angular velocities
            IAtomList children = molecule.getChildList();
            for (int i=0; i<children.getAtomCount(); i++) {
                IAtom atom = children.getAtom(i);
                Vector force = leafAgentManager.getAgent(atom).force;
                if (force.isZero()) {
                    continue;
                }

                xWork.Ev1Mv2(atom.getPosition(), orientedMolecule.getPosition());

                xWork.XE(force);
                agent.torque.PE(xWork);
                agent.force.PE(force);
            }
        }
    }
    
//--------------------------------------------------------------

    public Class getMoleculeAgentClass() {
        return MoleculeAgent.class;
    }

    public final Object makeAgent(IMolecule a) {
        return new MoleculeAgent(space);
    }

    public void releaseAgent(Object agent, IMolecule atom) {}

    public Class getAgentClass() {
        return MyAgent.class;
    }

    public final IntegratorVelocityVerlet.MyAgent makeAgent(IAtom a, Box agentBox) {
        return new MyAgent(space);
    }

    public void releaseAgent(IntegratorVelocityVerlet.MyAgent agent, IAtom atom, Box agentBox) {}

    public static class MoleculeAgent implements Integrator.Torquable, Integrator.Forcible, Serializable {  //need public so to use with instanceof
        private static final long serialVersionUID = 1L;
        public final Vector torque;
        public final Vector force;

        public MoleculeAgent(Space space) {
            torque = space.makeVector();
            force = space.makeVector();
        }
        
        public Vector torque() {return torque;}
        public Vector force() {return force;}
    }

    public Class getSpeciesAgentClass() {
        return OrientationCalc.class;
    }
    
    public Object makeAgent(ISpecies type) {
        return null;
    }

    public void releaseAgent(Object agent, ISpecies type) {}
    
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        boolean periodic = true;
        Box box;
        if (periodic) {
            box = new Box(new BoundaryRectangularPeriodic(sim.getSpace(), 10), space);
        }
        else {
            box = new Box(new BoundaryRectangularNonperiodic(sim.getSpace()), space);
        }
        sim.addBox(box);
        SpeciesWater3POriented species = new SpeciesWater3POriented(sim.getSpace(), true);
        sim.addSpecies(species);
        int numMolecules = 256;
        box.setNMolecules(species, numMolecules);
        box.setDensity(1/18.0*Constants.AVOGADRO/1E24);
        if (true) {
            new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        }
        else {
            new ConfigurationFile("water"+numMolecules+(periodic ? "pbc":"")).initializeCoordinates(box);
        }

        PotentialMaster potentialMaster = new PotentialMaster();
        double timeInterval = 0.001;
        int maxIterations = 20;
        IntegratorRigidIterative integrator = new IntegratorRigidIterative(sim, potentialMaster, timeInterval, 1, space);
        integrator.printInterval = 10;
        integrator.setMaxIterations(maxIterations);
        integrator.setBox(box);
        OrientationCalcWater3P calcer = new OrientationCalcWater3P(sim.getSpace());
        integrator.setOrientationCalc(species, calcer);
        integrator.setTemperature(Kelvin.UNIT.toSim(298));
        ActivityIntegrate ai = new ActivityIntegrate(integrator);
        sim.getController().addAction(ai);

        P2WaterSPCSoft p2Water = new P2WaterSPCSoft(sim.getSpace());
        
        if (periodic) {
            BoxImposePbc pbc = new BoxImposePbc(box, space);
            pbc.setApplyToMolecules(true);
            integrator.getEventManager().addListener(new IntegratorListenerAction(pbc));

            double boxlength = box.getBoundary().getBoxSize().getX(0);

            DipoleSourceWater dipoleSource = new DipoleSourceWater(sim.getSpace());
            dipoleSource.setDipoleStrength(2*Electron.UNIT.toSim(0.41)*Math.cos(109.5/2.0*Math.PI/180));
            IMoleculePositionDefinition positionDefinition = new MoleculePositionCOM(space) ;
            P2ReactionFieldDipole pNRF = new P2ReactionFieldDipole(sim.getSpace(),positionDefinition);
            pNRF.setDipoleSource(dipoleSource);
            pNRF.setRange(boxlength*0.5);
            pNRF.setDielectric(78.4);

            potentialMaster.addPotential(new P2MoleculeSoftTruncatedSwitched(pNRF, boxlength*0.49, space), new ISpecies[]{species, species});
            potentialMaster.lrcMaster().addPotential(pNRF.makeP0());

            potentialMaster.addPotential(new P2MoleculeSoftTruncatedSwitched(p2Water, boxlength*0.49, space), new ISpecies[]{species,species});
        }
        else {
            potentialMaster.addPotential(p2Water, new ISpecies[]{species,species});
        }
        
        if (true) {
            ai.setMaxSteps(1000000);
            ai.actionPerformed();
        }
        else {
            ai.setSleepPeriod(2);
            SimulationGraphic graphic = new SimulationGraphic(sim, "Rigid", 1, space, sim.getController());
            ((ColorSchemeByType)graphic.getDisplayBox(box).getColorScheme()).setColor(species.getHydrogenType(), Color.WHITE);
            ((ColorSchemeByType)graphic.getDisplayBox(box).getColorScheme()).setColor(species.getOxygenType(), Color.RED);
            graphic.makeAndDisplayFrame();
        }
    }
}
