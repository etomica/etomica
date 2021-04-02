/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.action.AtomActionTranslateBy;
import etomica.action.BoxImposePbc;
import etomica.action.MoleculeChildAtomAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.*;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationLattice;
import etomica.data.meter.MeterKineticEnergyRigid;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.lattice.LatticeCubicFcc;
import etomica.math.function.Function;
import etomica.models.water.DipoleSourceMolecularWater;
import etomica.models.water.OrientationCalcWater3P;
import etomica.models.water.P2WaterSPCSoft;
import etomica.models.water.SpeciesWater3P;
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
import etomica.species.SpeciesAgentManager;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesManager;
import etomica.units.Electron;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.Debug;
import etomica.util.random.IRandom;

import java.awt.*;
import java.io.Serializable;

public class IntegratorRigidIterative extends IntegratorMD implements SpeciesAgentManager.AgentSource<OrientationCalc>, MoleculeAgentSource<IntegratorRigidIterative.MoleculeAgent>, AtomLeafAgentManager.AgentSource<IntegratorRigidIterative.AtomForceAgent> {

    protected PotentialCalculationTorqueSum torqueSum;
    private final IteratorDirective allAtoms;
    protected final Tensor pressureTensor;
    protected final Tensor workTensor;
    protected final RotationTensor3D rotationTensor;
    protected final RotationTensor3D tempRotationTensor;
    protected final Vector xWork;
    protected final SpeciesAgentManager<OrientationCalc> typeAgentManager;
    protected final Vector tempAngularVelocity;
    protected final MoleculePositionCOM atomPositionCOM;
    protected final AtomActionTranslateBy translateBy;
    protected final MoleculeChildAtomAction translator;
    protected final OrientationFull3D tempOrientation;
    public int printInterval;
    protected int maxIterations;
    protected double omegaTolerance;
    protected final boolean storeAngularMomentum = false;

    protected AtomLeafAgentManager<AtomForceAgent> leafAgentManager;
    protected MoleculeAgentManager<MoleculeAgent> moleculeAgentManager;

    public IntegratorRigidIterative(SpeciesManager sm, IRandom random, PotentialMaster potentialMaster, Box box) {
        this(sm, random, potentialMaster, 0.05, 1.0, box);
    }

    public IntegratorRigidIterative(SpeciesManager sm, IRandom random, PotentialMaster potentialMaster,
                                    double timeStep, double temperature, Box box) {
        super(potentialMaster, random, timeStep, temperature, box);
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
        rotationTensor = (RotationTensor3D) space.makeRotationTensor();
        tempRotationTensor = (RotationTensor3D)space.makeRotationTensor();
        xWork = space.makeVector();
        typeAgentManager = new SpeciesAgentManager<>(this, sm);
        tempAngularVelocity = space.makeVector();
        tempOrientation = new OrientationFull3D(space);
        atomPositionCOM = new MoleculePositionCOM(space);
        translateBy = new AtomActionTranslateBy(space);
        translator = new MoleculeChildAtomAction(translateBy);
        maxIterations = 40;
        omegaTolerance = 1.e-30;
        meterKE = new MeterKineticEnergyRigid(space, sm);

        leafAgentManager = new AtomLeafAgentManager<>(this, box);
        moleculeAgentManager = new MoleculeAgentManager<>(sm, box, this);
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
    protected void doStepInternal() {
        super.doStepInternal();
        currentKineticEnergy = 0;
        int iterationsTotal = 0;
        int numRigid = 0;
        if (Debug.ON && Debug.DEBUG_NOW) {
            IAtomList pair = Debug.getAtoms(box);
            if (pair != null) {
                Vector dr = space.makeVector();
                dr.Ev1Mv2(pair.get(1).getPosition(), pair.get(0).getPosition());
                System.out.println(pair+" dr "+dr);
            }
        }
        IMoleculeList moleculeList = box.getMoleculeList();
        int nMolecules = moleculeList.size();
        for (int iMolecule = 0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = moleculeList.get(iMolecule);
            IAtomList children = molecule.getChildList();
            OrientationCalc calcer = (OrientationCalc)typeAgentManager.getAgent(molecule.getType());
            if (calcer == null) {
                for (int iLeaf = 0; iLeaf<children.size(); iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic)children.get(iLeaf);
                    Vector force = leafAgentManager.getAgent(a).force();
                    Vector r = a.getPosition();
                    Vector v = a.getVelocity();
                    if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                        System.out.println("first "+a+" r="+r+", v="+v+", f="+force);
                    }
                    v.PEa1Tv1(0.5*timeStep* a.getType().rm(),force);  // p += f(old)*dt/2
                    r.PEa1Tv1(timeStep,v);         // r += p*dt/m
                }
                continue;
            }
            
            MoleculeAgent agent = (MoleculeAgent)moleculeAgentManager.getAgent(molecule);
            IMoleculeOrientedKinetic orientedMolecule = (IMoleculeOrientedKinetic)molecule;
            IOrientationFull3D orientation = (IOrientationFull3D)orientedMolecule.getOrientation();
            Vector moment = molecule.getType().getMomentOfInertia();

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
            orientedMolecule.getVelocity().PEa1Tv1(0.5*timeStep/molecule.getType().getMass(), agent.force);
            
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
            IMolecule molecule = moleculeList.get(iMolecule);
            IAtomList children = molecule.getChildList();
            if (typeAgentManager.getAgent(molecule.getType()) == null) {
                // unimolecular or at least not rigid
                //Finish integration step
                for (int iLeaf = 0; iLeaf<children.size(); iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic)children.get(iLeaf);
//                    System.out.println("force: "+((MyAgent)a.ia).force.toString());
                    Vector velocity = a.getVelocity();
                    workTensor.Ev1v2(velocity,velocity);
                    workTensor.TE(a.getType().getMass());
                    pressureTensor.PE(workTensor);
                    if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                        System.out.println("second "+a+" v="+velocity+", f="+leafAgentManager.getAgent(a));
                    }
                    velocity.PEa1Tv1(0.5 * timeStep * a.getType().rm(), leafAgentManager.getAgent(a).force());  //p += f(new)*dt/2
                    currentKineticEnergy += velocity.squared()* a.getType().getMass();
                }
                // skip the rotational stuff
                continue;
            }
            
            IMoleculeOrientedKinetic orientedMolecule = (IMoleculeOrientedKinetic)molecule;
            MoleculeAgent agent = (MoleculeAgent)moleculeAgentManager.getAgent(molecule);
            //calc torque and linear force
            for (int i = 0; i<children.size(); i++) {
                IAtom atom = children.get(i);
                Vector atomForce = leafAgentManager.getAgent(atom).force();
                if (atomForce.isZero()) {
                    continue;
                }

                agent.force.PE(atomForce);

                xWork.Ev1Mv2(atom.getPosition(), orientedMolecule.getPosition());
                xWork.XE(atomForce);
                agent.torque.PE(xWork);
            }
            
            //advance linear velocity to full timestep
            double mass = molecule.getType().getMass();
            orientedMolecule.getVelocity().PEa1Tv1(0.5*timeStep/mass, agent.force);

            //advance momentum to full timestep
            Vector angularVelocity = orientedMolecule.getAngularVelocity();
            Vector moment = molecule.getType().getMomentOfInertia();

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
            int moleculeCount = box.getMoleculeList().size();
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
        int nMolecules = moleculeList.size();
        int D = 0;
        momentum.E(0);
        double totalMass = 0;
        for (int iMolecule = 0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = moleculeList.get(iMolecule);
            IAtomList children = molecule.getChildList();
            if (typeAgentManager.getAgent(molecule.getType()) == null) {
                for (int iLeaf = 0; iLeaf<children.size(); iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic)children.get(iLeaf);
                    double mass = a.getType().getMass();
                    momentum.PEa1Tv1(mass, a.getVelocity());
                    totalMass += mass;
                }
                continue;
            }
            IMoleculeOrientedKinetic orientedMolecule = (IMoleculeOrientedKinetic)molecule;
            double mass = ((IMolecule)orientedMolecule).getType().getMass();
            momentum.PEa1Tv1(mass, orientedMolecule.getVelocity());
            totalMass += mass;
        }
        momentum.TE(1.0/totalMass);
        for (int iMolecule = 0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = moleculeList.get(iMolecule);
            IAtomList children = molecule.getChildList();
            if (typeAgentManager.getAgent(molecule.getType()) == null) {
                // unimolecular or at least not rigid
                //Finish integration step
                for (int iLeaf = 0; iLeaf<children.size(); iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic)children.get(iLeaf);
//                    System.out.println("force: "+((MyAgent)a.ia).force.toString());
                    Vector velocity = a.getVelocity();
                    velocity.ME(momentum);
                    KE += velocity.squared() * a.getType().getMass();
                    D += 3;
                }
                continue;
            }
            IMoleculeOrientedKinetic orientedMolecule = (IMoleculeOrientedKinetic)molecule;
            orientedMolecule.getVelocity().ME(momentum);
            KE += orientedMolecule.getVelocity().squared() * molecule.getType().getMass();

            Vector moment = molecule.getType().getMomentOfInertia();

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
            IMolecule molecule = moleculeList.get(iMolecule);
            IAtomList children = molecule.getChildList();
            if (typeAgentManager.getAgent(molecule.getType()) == null) {
                // unimolecular or at least not rigid
                //Finish integration step
                for (int iLeaf = 0; iLeaf<children.size(); iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic)children.get(iLeaf);
//                    System.out.println("force: "+((MyAgent)a.ia).force.toString());
                    a.getVelocity().TE(scale);
                }
                continue;
            }
            IMoleculeOrientedKinetic orientedMolecule = (IMoleculeOrientedKinetic)molecule;
            orientedMolecule.getVelocity().TE(scale);

            Vector angularVelocity = orientedMolecule.getAngularVelocity();
            rotationTensor.setOrientation((IOrientationFull3D)orientedMolecule.getOrientation());
            rotationTensor.transform(angularVelocity);
            angularVelocity.TE(scale);
            rotationTensor.invert();
            rotationTensor.transform(angularVelocity);
        }
    }

    public void shiftMomenta() {
        momentum.E(0);
        IMoleculeList moleculeList = box.getMoleculeList();
        int nMolecules = moleculeList.size();
        if (nMolecules == 0) return;
        if (nMolecules > 1) {
            double totalMass = 0;
            for (int iLeaf = 0; iLeaf < nMolecules; iLeaf++) {
                IMolecule m = moleculeList.get(iLeaf);
                if (!(m instanceof IMoleculeKinetic)) {
                    IAtomList children = m.getChildList();
                    for (int i = 0; i < children.size(); i++) {
                        IAtomKinetic a = (IAtomKinetic)children.get(i);
                        momentum.PEa1Tv1(a.getType().getMass(), a.getVelocity());
                        totalMass += a.getType().getMass();
                    }
                    continue;
                }
                double mass = m.getType().getMass();
                if (mass != Double.POSITIVE_INFINITY) {
                    momentum.PEa1Tv1(mass, ((IMoleculeKinetic) m).getVelocity());
                    totalMass += mass;
                }
            }
            if (totalMass == 0) return;
            momentum.TE(1.0 / totalMass);
            //momentum is now net velocity
            //set net momentum to 0
            for (int iLeaf = 0; iLeaf < nMolecules; iLeaf++) {
                IMolecule m = moleculeList.get(iLeaf);
                if (!(m instanceof IMoleculeKinetic)) {
                    IAtomList children = m.getChildList();
                    for (int i = 0; i < children.size(); i++) {
                        IAtomKinetic a = (IAtomKinetic)children.get(i);
                        double rm = a.getType().rm();
                        if (rm != 0 && rm != Double.POSITIVE_INFINITY) {
                            a.getVelocity().ME(momentum);
                        }
                    }
                    continue;
                }
                double mass = m.getType().getMass();
                if (mass != Double.POSITIVE_INFINITY) {
                    ((IMoleculeKinetic)m).getVelocity().ME(momentum);
                }
            }
            if (Debug.ON) {
                momentum.E(0);
                for (int iLeaf = 0; iLeaf < nMolecules; iLeaf++) {
                    IMolecule m = moleculeList.get(iLeaf);
                    if (!(m instanceof IMoleculeKinetic)) {
                        IAtomList children = m.getChildList();
                        for (int i = 0; i < children.size(); i++) {
                            IAtomKinetic a = (IAtomKinetic)children.get(i);
                            momentum.PEa1Tv1(a.getType().getMass(), a.getVelocity());
                        }
                        continue;
                    }
                    double mass = m.getType().getMass();
                    if (mass != Double.POSITIVE_INFINITY) {
                        momentum.PEa1Tv1(mass, ((IMoleculeKinetic) m).getVelocity());
                    }
                }
                momentum.TE(1.0 / totalMass);
                if (Math.sqrt(momentum.squared()) > 1.e-10) {
                    System.out.println("Net momentum per leaf atom is " + momentum + " but I expected it to be 0");
                }
            }
            momentum.E(0);
        }

    }

    public void randomizeMomenta() {
        IMoleculeList moleculeList = box.getMoleculeList();
        int nMolecules = moleculeList.size();
//        System.out.println("rerandomize");
        for (int iMolecule = 0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = moleculeList.get(iMolecule);
            OrientationCalc calcer = (OrientationCalc)typeAgentManager.getAgent(molecule.getType());
            if (calcer == null) {
                IAtomList children = molecule.getChildList();
                for (int i = 0; i<children.size(); i++) {
                    super.randomizeMomentum((IAtomKinetic)children.get(i));
                }
                continue;
            }

            MoleculeOrientedDynamic orientedMolecule = (MoleculeOrientedDynamic)molecule;
            Vector velocity = orientedMolecule.getVelocity();
            Vector moment = molecule.getType().getMomentOfInertia();

            double mass = molecule.getType().getMass();
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
            double mass = orientedMolecule.getType().getMass();
            int D = velocity.getD();
            for(int i=0; i<D; i++) {
                velocity.setX(i,random.nextGaussian());
            }
            velocity.TE(Math.sqrt(temperature/mass));

            Vector angularVelocity = orientedMolecule.getAngularVelocity();
            Vector moment = orientedMolecule.getType().getMomentOfInertia();
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
                dr.Ev1Mv2(pair.get(1).getPosition(), pair.get(0).getPosition());
                System.out.println(pair+" dr "+dr);
            }
        }
        
        IMoleculeList moleculeList = box.getMoleculeList();
        int nMolecules = moleculeList.size();

        for (int iMolecule=0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = moleculeList.get(iMolecule);
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
            IMolecule molecule = moleculeList.get(iMolecule);
            if (typeAgentManager.getAgent(molecule.getType()) == null) {
                continue;
            }

            MoleculeAgent agent = (MoleculeAgent)moleculeAgentManager.getAgent(molecule);
            IMoleculeOrientedKinetic orientedMolecule = (IMoleculeOrientedKinetic)molecule;
            //calc angular velocities
            IAtomList children = molecule.getChildList();
            for (int i = 0; i<children.size(); i++) {
                IAtom atom = children.get(i);
                Vector force = leafAgentManager.getAgent(atom).force();
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

    public final MoleculeAgent makeAgent(IMolecule a) {
        return new MoleculeAgent(space);
    }

    public void releaseAgent(MoleculeAgent agent, IMolecule atom) {
    }

    @Override
    public AtomForceAgent makeAgent(IAtom a, Box agentBox) {
        return new AtomForceAgent(space);
    }

    @Override
    public void releaseAgent(AtomForceAgent agent, IAtom atom, Box agentBox) {
    }

    public static class AtomForceAgent implements Integrator.Forcible, Serializable {  //need public so to use with instanceof
        public final Vector force;

        public AtomForceAgent(Space space) {
            force = space.makeVector();
        }

        public Vector force() {
            return force;
        }
    }

    public static class MoleculeAgent implements Integrator.Torquable, Integrator.Forcible, Serializable {  //need public so to use with instanceof
        public final Vector torque;
        public final Vector force;

        public MoleculeAgent(Space space) {
            torque = space.makeVector();
            force = space.makeVector();
        }
        
        public Vector torque() {return torque;}
        public Vector force() {return force;}
    }

    public OrientationCalc makeAgent(ISpecies type) {
        return null;
    }

    public void releaseAgent(OrientationCalc agent, ISpecies type) {
    }
    
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        SpeciesGeneral species = SpeciesWater3P.create(true, true);
        sim.addSpecies(species);
        boolean periodic = true;
        Box box;
        if (periodic) {
            box = new Box(new BoundaryRectangularPeriodic(sim.getSpace(), 10), space);
        }
        else {
            box = new Box(new BoundaryRectangularNonperiodic(sim.getSpace()), space);
        }
        sim.addBox(box);
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
        IntegratorRigidIterative integrator = new IntegratorRigidIterative(sim.getSpeciesManager(), sim.getRandom(), potentialMaster, timeInterval, 1, box);
        integrator.printInterval = 10;
        integrator.setMaxIterations(maxIterations);
        OrientationCalcWater3P calcer = new OrientationCalcWater3P(sim.getSpace());
        integrator.setOrientationCalc(species, calcer);
        integrator.setTemperature(Kelvin.UNIT.toSim(298));

        P2WaterSPCSoft p2Water = new P2WaterSPCSoft(sim.getSpace());
        
        if (periodic) {
            BoxImposePbc pbc = new BoxImposePbc(box, space);
            pbc.setApplyToMolecules(true);
            integrator.getEventManager().addListener(new IntegratorListenerAction(pbc));

            double boxlength = box.getBoundary().getBoxSize().getX(0);

            DipoleSourceMolecularWater dipoleSource = new DipoleSourceMolecularWater(sim.getSpace());
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
            sim.getController().runActivityBlocking(new ActivityIntegrate(integrator, 1000000));
        }
        else {
            sim.getController().setSleepPeriod(2);
            sim.getController().addActivity(new ActivityIntegrate(integrator));
            SimulationGraphic graphic = new SimulationGraphic(sim, "Rigid", 1);
            ((ColorSchemeByType) graphic.getDisplayBox(box).getColorScheme()).setColor(species.getTypeByName("H"), Color.WHITE);
            ((ColorSchemeByType) graphic.getDisplayBox(box).getColorScheme()).setColor(species.getTypeByName("O"), Color.RED);
            graphic.makeAndDisplayFrame();
        }
    }
}
