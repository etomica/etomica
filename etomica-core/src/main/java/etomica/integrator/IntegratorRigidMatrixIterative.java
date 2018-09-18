/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.action.AtomActionTranslateBy;
import etomica.action.IAction;
import etomica.action.MoleculeChildAtomAction;
import etomica.atom.*;
import etomica.box.Box;
import etomica.data.meter.MeterKineticEnergyRigid;
import etomica.math.function.Function;
import etomica.molecule.*;
import etomica.molecule.MoleculeAgentManager.MoleculeAgentSource;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationTorqueSum;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.IOrientationFull3D;
import etomica.space3d.OrientationFull3D;
import etomica.space3d.RotationTensor3D;
import etomica.species.ISpecies;
import etomica.species.ISpeciesOriented;
import etomica.species.SpeciesAgentManager;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.Debug;

/**
 * Integrator implementation of Omelyan's leapfrog integrator for rotational
 * motion.  Molecular Simulation, 22 (1999) 213-236.
 * 
 * @author Andrew Schultz
 */
public class IntegratorRigidMatrixIterative extends IntegratorMD implements SpeciesAgentManager.AgentSource, MoleculeAgentSource {

    private static final long serialVersionUID = 2L;
    protected PotentialCalculationTorqueSum torqueSum;
    private final IteratorDirective allAtoms;
    protected final Tensor pressureTensor;
    protected final Tensor workTensor;
    protected final RotationTensor3D rotationTensor;
    protected final Vector xWork, yWork;
    protected final SpeciesAgentManager typeAgentManager;
    protected final Vector tempAngularVelocity;
    protected final MoleculePositionCOM atomPositionCOM;
    protected final AtomActionTranslateBy translateBy;
    protected final MoleculeChildAtomAction translator;
    protected final OrientationFull3D tempOrientation;
    public int printInterval;
    protected int maxIterations;
    protected double omegaTolerance, rotationTolerance;
    protected final RotationTensor3D axesTensor;
    protected final RotationTensor3D tempAxesTensor;
    protected final Tensor omegaTensor;
    protected final Simulation sim;
    
    protected AtomLeafAgentManager<Vector> forces;
    protected MoleculeAgentManager moleculeAgentManager;

    public IntegratorRigidMatrixIterative(Simulation sim, PotentialMaster potentialMaster, Box box) {
        this(sim, potentialMaster, 0.05, 1.0, box);
    }
    
    public IntegratorRigidMatrixIterative(Simulation sim, PotentialMaster potentialMaster,
                                          double timeStep, double temperature, Box box) {
        super(potentialMaster,sim.getRandom(),timeStep,temperature, box);
        this.sim = sim;
        torqueSum = new PotentialCalculationTorqueSum();
        allAtoms = new IteratorDirective();
        // allAtoms is used only for the force calculation, which has no LRC
        // but we're also calculating the pressure tensor, which does have LRC.
        // things deal with this OK.
        allAtoms.setIncludeLrc(true);
        pressureTensor = this.space.makeTensor();
        workTensor = this.space.makeTensor();
        rotationTensor = (RotationTensor3D)this.space.makeRotationTensor();
        xWork = this.space.makeVector();
        yWork = this.space.makeVector();
        typeAgentManager = new SpeciesAgentManager(this, sim);
        tempAngularVelocity = this.space.makeVector();
        tempOrientation = new OrientationFull3D(this.space);
        atomPositionCOM = new MoleculePositionCOM(this.space);
        translateBy = new AtomActionTranslateBy(this.space);
        translator = new MoleculeChildAtomAction(translateBy);
        maxIterations = 20;
        omegaTolerance = 1.e-30;
        rotationTolerance = 1.e-30;
        axesTensor = (RotationTensor3D)this.space.makeRotationTensor();
        tempAxesTensor = (RotationTensor3D)this.space.makeRotationTensor();
        omegaTensor = this.space.makeTensor();
        meterKE = new MeterKineticEnergyRigid(space, sim);

        forces = new AtomLeafAgentManager<>(a -> space.makeVector(), box);
        moleculeAgentManager = new MoleculeAgentManager(sim, this.box, this);
        torqueSum.setAgentManager(forces);
        torqueSum.setMoleculeAgentManager(moleculeAgentManager);
        ((MeterKineticEnergyRigid)meterKE).setBox(box);
    }

    public void setOrientationCalc(ISpecies species, OrientationCalc calcer) {
        typeAgentManager.setAgent(species, calcer);
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
            OrientationCalc calcer = (OrientationCalc)typeAgentManager.getAgent(molecule.getType());
            IAtomList children = molecule.getChildList();
            if (calcer == null) {
                for (int iLeaf = 0; iLeaf<children.size(); iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic)children.get(iLeaf);
                    Vector force = forces.getAgent(a);
                    Vector r = a.getPosition();
                    Vector v = a.getVelocity();
                    if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                        System.out.println("first "+a+" r="+r+", v="+v+", f="+ force);
                    }
                    v.PEa1Tv1(0.5*timeStep* a.getType().rm(), force);  // p += f(old)*dt/2
                    r.PEa1Tv1(timeStep,v);         // r += p*dt/m
                }
                continue;
            }

            IntegratorRigidIterative.MoleculeAgent agent = (IntegratorRigidIterative.MoleculeAgent) moleculeAgentManager.getAgent(molecule);
            IMoleculeOrientedKinetic orientedMolecule = (IMoleculeOrientedKinetic)molecule;
            Vector moment = ((ISpeciesOriented)molecule.getType()).getMomentOfInertia();
            double mass = ((ISpeciesOriented)molecule.getType()).getMass();
            IOrientationFull3D orientation = (IOrientationFull3D)orientedMolecule.getOrientation();

            // use the angular velocity field to store angular momentum during the time step  :(
            Vector angularMomentum = orientedMolecule.getAngularVelocity();

            // transform to body-fixed, multiply by moment of inertia, transform back
            rotationTensor.setOrientation(orientation);
            axesTensor.E(rotationTensor);
            rotationTensor.transform(angularMomentum);
            
            omegaTensor.E(0);
            omegaTensor.setComponent(0, 1, angularMomentum.getX(2));
            omegaTensor.setComponent(0, 2, -angularMomentum.getX(1));
            omegaTensor.setComponent(1, 0, -angularMomentum.getX(2));
            omegaTensor.setComponent(1, 2, angularMomentum.getX(0));
            omegaTensor.setComponent(2, 0, angularMomentum.getX(1));
            omegaTensor.setComponent(2, 1, -angularMomentum.getX(0));
            
            // xWork = angular velocity from the previous iteration
            // save angular velocity to see if we've converged
            xWork.E(angularMomentum);
            
            angularMomentum.TE(moment);
            rotationTensor.invert();
            rotationTensor.transform(angularMomentum);
            rotationTensor.invert();

            angularMomentum.PEa1Tv1(0.5*timeStep, agent.torque);

            numRigid++;
            for (int i = 0; i<maxIterations; i++) {
                iterationsTotal++;
                tempAxesTensor.E(omegaTensor);
                tempAxesTensor.TE(rotationTensor);
                tempAxesTensor.TE(timeStep);
                tempAxesTensor.TE(0.5);
                tempAxesTensor.PE(axesTensor);
                // tempAxesTensor is now our estimate of A at the half timestep
                
                //now find half-timestep angular velocity
                tempAngularVelocity.E(angularMomentum);
                // invert so we can go space -> body
                tempAxesTensor.transform(tempAngularVelocity);
                tempAngularVelocity.DE(moment);
                // tempAngularVelocity is now the half-timestep angular velocity
                xWork.ME(tempAngularVelocity);
//                System.out.println("omegaDiff "+i+" "+xWork.squared()/tempAngularVelocity.squared());
                if (xWork.squared()/tempAngularVelocity.squared() < omegaTolerance) {
                    rotationTensor.ME(tempAxesTensor);
                    double sumDiff = 0, sum = 0;
                    for (int j=0; j<3; j++) {
                        for (int k=0; k<3; k++) {
                            double v = rotationTensor.component(j,k);
                            sumDiff += v*v;
                            v = tempAxesTensor.component(j,k);
                            sum += v*v;
                        }
                    }
//                    System.out.println("axes diff "+sum);
                    if (sumDiff/sum < rotationTolerance) {
                        rotationTensor.E(tempAxesTensor);
                        break;
                    }
                }
                xWork.E(tempAngularVelocity);
                rotationTensor.E(tempAxesTensor);
                
                // populate omegaTensor again with our updated estimate of omega
                omegaTensor.E(0);
                omegaTensor.setComponent(0, 1, tempAngularVelocity.getX(2));
                omegaTensor.setComponent(0, 2, -tempAngularVelocity.getX(1));
                omegaTensor.setComponent(1, 0, -tempAngularVelocity.getX(2));
                omegaTensor.setComponent(1, 2, tempAngularVelocity.getX(0));
                omegaTensor.setComponent(2, 0, tempAngularVelocity.getX(1));
                omegaTensor.setComponent(2, 1, -tempAngularVelocity.getX(0));
            }

            tempAxesTensor.E(omegaTensor);
            tempAxesTensor.TE(rotationTensor);
            tempAxesTensor.TE(timeStep);
            tempAxesTensor.PE(axesTensor);
            
            xWork.setX(0, tempAxesTensor.component(0,0));
            xWork.setX(1, tempAxesTensor.component(0,1));
            xWork.setX(2, tempAxesTensor.component(0,2));
            yWork.setX(0, tempAxesTensor.component(1,0));
            yWork.setX(1, tempAxesTensor.component(1,1));
            yWork.setX(2, tempAxesTensor.component(1,2));
            orientation.setDirections(xWork, yWork);

//            System.out.println("o "+orientation.getDirection()+" "+orientation.getSecondaryDirection());
            calcer.setOrientation(molecule, orientation);
            //advance linear velocity to half timestep
            orientedMolecule.getVelocity().PEa1Tv1(0.5*timeStep/mass, agent.force);
            
            //advance position to full timestep
            Vector transVec = ((AtomActionTranslateBy)translator.getAtomAction()).getTranslationVector();
            transVec.Ea1Tv1(timeStep, orientedMolecule.getVelocity());
            orientedMolecule.getPosition().PE(transVec);
            translator.actionPerformed(molecule);
        }

//        if(forceSum instanceof PotentialCalculationForcePressureSum){
//            pressureTensor.E(((PotentialCalculationForcePressureSum)forceSum).getPressureTensor());
//        }
        if (isothermal) {
            doThermostatInternal();
        }
        
        torqueSum.reset();
        //Compute forces and torques on each atom at t+dt
        potentialMaster.calculate(box, allAtoms, torqueSum);
        
        for (int iMolecule = 0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = moleculeList.get(iMolecule);
            IAtomList children = molecule.getChildList();
            OrientationCalc calcer = (OrientationCalc)typeAgentManager.getAgent(molecule.getType());
            if (calcer == null) {
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
                        System.out.println("second "+a+" v="+velocity+", f="+ forces.getAgent(a));
                    }
                    velocity.PEa1Tv1(0.5*timeStep* a.getType().rm(), forces.getAgent(a));  //p += f(new)*dt/2
                    currentKineticEnergy += velocity.squared()* a.getType().getMass();
                }
                // skip the rotational stuff
                continue;
            }
            
            IMoleculeOrientedKinetic orientedMolecule = (IMoleculeOrientedKinetic)molecule;
            IntegratorRigidIterative.MoleculeAgent agent = (IntegratorRigidIterative.MoleculeAgent) moleculeAgentManager.getAgent(molecule);
            Vector moment = ((ISpeciesOriented)molecule.getType()).getMomentOfInertia();
            double mass = ((ISpeciesOriented)molecule.getType()).getMass();
            //calc torque and linear force
            for (int i = 0; i<children.size(); i++) {
                IAtom atom = children.get(i);
                Vector atomForce = forces.getAgent(atom);
                if (atomForce.isZero()) {
                    continue;
                }

                agent.force.PE(atomForce);

                xWork.Ev1Mv2(atom.getPosition(), orientedMolecule.getPosition());
                xWork.XE(atomForce);
                agent.torque.PE(xWork);
            }
            
            //advance linear velocity to full timestep
            orientedMolecule.getVelocity().PEa1Tv1(0.5*timeStep/mass, agent.force);

            //advance momentum to full timestep
            Vector angularVelocity = orientedMolecule.getAngularVelocity();

            // we actually stored the half-timestep angular momentum in this field...
            // advance to full timestep
            angularVelocity.PEa1Tv1(0.5*timeStep, agent.torque);
            // we actually had this before, but we didn't save it...
            // transform to body-fixed, multiply by moment of inertia, transform back
            // ==> angular momentum
            rotationTensor.setOrientation((IOrientationFull3D)orientedMolecule.getOrientation());
            rotationTensor.transform(angularVelocity);
            angularVelocity.DE(moment);
            xWork.E(angularVelocity);
            rotationTensor.invert();
            rotationTensor.transform(angularVelocity);

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
        if (nMolecules == 0) {
            return;
        }
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
            double mass = ((ISpeciesOriented)((IMolecule)orientedMolecule).getType()).getMass();
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
                    totalMass += a.getType().getMass();
                    KE += velocity.squared() * a.getType().getMass();
                    D += 3;
                }
                continue;
            }
            IMoleculeOrientedKinetic orientedMolecule = (IMoleculeOrientedKinetic)molecule;
            Vector velocity = orientedMolecule.getVelocity();
            velocity.ME(momentum);
            KE += velocity.squared() * ((ISpeciesOriented)molecule.getType()).getMass();

            Vector moment = ((ISpeciesOriented)molecule.getType()).getMomentOfInertia();

            tempAngularVelocity.E(orientedMolecule.getAngularVelocity());
            rotationTensor.setOrientation((IOrientationFull3D)orientedMolecule.getOrientation());
            rotationTensor.transform(tempAngularVelocity);
            tempAngularVelocity.TE(moment);

            tempAngularVelocity.TE(tempAngularVelocity);
            tempAngularVelocity.DE(moment);
            KE += tempAngularVelocity.getX(0) + tempAngularVelocity.getX(1)+ tempAngularVelocity.getX(2);
            D += 6;
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
                double mass = ((ISpeciesOriented)m.getType()).getMass();
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
                double mass = ((ISpeciesOriented)m.getType()).getMass();
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
                    double mass = ((ISpeciesOriented)m.getType()).getMass();
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
            double mass = ((ISpeciesOriented)molecule.getType()).getMass();
            Vector moment = ((ISpeciesOriented)molecule.getType()).getMomentOfInertia();
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
            angularVelocity.DE(moment);
            
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
            double mass = ((ISpeciesOriented)((IMolecule)atom).getType()).getMass();
            Vector moment = ((ISpeciesOriented)((IMolecule)atom).getType()).getMomentOfInertia();
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

            IntegratorRigidIterative.MoleculeAgent agent = (IntegratorRigidIterative.MoleculeAgent) moleculeAgentManager.getAgent(molecule);
            MoleculeOrientedDynamic orientedMolecule = (MoleculeOrientedDynamic)molecule;
            //calc angular velocities
            IAtomList children = molecule.getChildList();
            for (int i = 0; i<children.size(); i++) {
                IAtom atom = children.get(i);
                Vector force = forces.getAgent(atom);
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

    public final Object makeAgent(IMolecule a) {
        return new IntegratorRigidIterative.MoleculeAgent(space);
    }
    
    public void releaseAgent(Object agent, IMolecule molecule) {}

    public Object makeAgent(ISpecies type) {
        return null;
    }

    public void releaseAgent(Object agent, ISpecies type) {}
    
    public static class BoxImposePbcMolecule implements IAction {
        public BoxImposePbcMolecule(Box box, Space space) {
            this.box = box;
            translateBy = new AtomActionTranslateBy(space);
            translator = new MoleculeChildAtomAction(translateBy);
            positionDefinition = new MoleculePositionGeometricCenter(space);
        }
        public void actionPerformed() {
            IMoleculeList molecules = box.getMoleculeList();
            Boundary boundary = box.getBoundary();
            for (int i = 0; i<molecules.size(); i++) {
                IMolecule molecule = molecules.get(i);
                if (molecule instanceof MoleculeOrientedDynamic) {
                    Vector position = ((MoleculeOrientedDynamic)molecule).getPosition();
                    Vector shift = boundary.centralImage(position);
                    if (!shift.isZero()) {
                        translateBy.setTranslationVector(shift);
                        translator.actionPerformed(molecule);
                    }
                    position.PE(shift);
                }
                else {
                    Vector position = positionDefinition.position(molecule);
                    Vector shift = boundary.centralImage(position);
                    if (!shift.isZero()) {
                        translateBy.setTranslationVector(shift);
                        translator.actionPerformed(molecule);
                    }
                }
            }
        }
        protected Box box;
        protected final AtomActionTranslateBy translateBy;
        protected final MoleculeChildAtomAction translator;
        protected final IMoleculePositionDefinition positionDefinition;
    }
}
