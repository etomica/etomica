package etomica.integrator;

import java.awt.Color;
import java.io.Serializable;

import etomica.EtomicaInfo;
import etomica.action.AtomActionTranslateBy;
import etomica.action.AtomGroupAction;
import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtom;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPositionCOM;
import etomica.atom.AtomSetSinglet;
import etomica.atom.AtomTypeAgentManager;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.AtomTypeMoleculeOriented;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomOrientedKinetic;
import etomica.atom.MoleculeOrientedDynamic;
import etomica.atom.OrientationCalc;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationLattice;
import etomica.data.meter.MeterKineticEnergyRigid;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet.MyAgent;
import etomica.lattice.LatticeCubicFcc;
import etomica.models.water.DipoleSourceWater;
import etomica.models.water.OrientationCalcWater;
import etomica.models.water.P2WaterSPCSoft;
import etomica.models.water.SpeciesWater3POriented;
import etomica.potential.P2MoleculeSoftTruncatedSwitched;
import etomica.potential.P2ReactionFieldDipole;
import etomica.potential.PotentialCalculationTorqueSum;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space3d.IOrientationFull3D;
import etomica.space3d.IVector3D;
import etomica.space3d.OrientationFull3D;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Space3D;
import etomica.units.Electron;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.Debug;
import etomica.util.Function;

public class IntegratorRigidIterative extends IntegratorMD implements AgentSource, AtomTypeAgentManager.AgentSource {

    private static final long serialVersionUID = 2L;
    protected PotentialCalculationTorqueSum torqueSum;
    private final IteratorDirective allAtoms;
    protected final Tensor pressureTensor;
    protected final Tensor workTensor;
    protected final RotationTensor3D rotationTensor;
    protected final RotationTensor3D tempRotationTensor;
    protected final IVector3D xWork;
    protected final AtomTypeAgentManager typeAgentManager;
    protected final IVector tempAngularVelocity;
    protected final AtomPositionCOM atomPositionCOM;
    protected final AtomActionTranslateBy translateBy;
    protected final AtomGroupAction translator;
    protected final OrientationFull3D tempOrientation;
    protected IAtomPositioned orientAtom;
    public int printInterval;
    protected int maxIterations;
    protected double omegaTolerance;
    protected final boolean storeAngularMomentum = false;
    
    protected AtomAgentManager agentManager;

    public IntegratorRigidIterative(ISimulation sim, PotentialMaster potentialMaster, Space _space) {
        this(sim, potentialMaster, 0.05, 1.0, _space);
    }
    
    public IntegratorRigidIterative(ISimulation sim, PotentialMaster potentialMaster,
            double timeStep, double temperature, Space _space) {
        super(potentialMaster,sim.getRandom(),timeStep,temperature, _space);
        // if you're motivated to throw away information earlier, you can use 
        // PotentialCalculationForceSum instead.
        torqueSum = new PotentialCalculationTorqueSum();
        allAtoms = new IteratorDirective();
        // allAtoms is used only for the force calculation, which has no LRC
        // but we're also calculating the pressure tensor, which does have LRC.
        // things deal with this OK.
        allAtoms.setIncludeLrc(true);
        pressureTensor = potentialMaster.getSpace().makeTensor();
        workTensor = potentialMaster.getSpace().makeTensor();
        rotationTensor = (RotationTensor3D)potentialMaster.getSpace().makeRotationTensor();
        tempRotationTensor = (RotationTensor3D)potentialMaster.getSpace().makeRotationTensor();
        xWork = (IVector3D)potentialMaster.getSpace().makeVector();
        typeAgentManager = new AtomTypeAgentManager(this, sim.getSpeciesManager(), sim.getEventManager(), true);
        tempAngularVelocity = potentialMaster.getSpace().makeVector();
        tempOrientation = new OrientationFull3D(potentialMaster.getSpace());
        atomPositionCOM = new AtomPositionCOM(potentialMaster.getSpace());
        translateBy = new AtomActionTranslateBy(potentialMaster.getSpace());
        translator = new AtomGroupAction(translateBy);
        maxIterations = 20;
        omegaTolerance = 1.e-25;
        meterKE = new MeterKineticEnergyRigid(potentialMaster.getSpace(), sim);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Molecular dynamics using velocity Verlet integration algorithm");
        return info;
    }
    
    public void setBox(IBox p) {
        if (box != null) {
            // allow agentManager to de-register itself as a BoxListener
            agentManager.dispose();
        }
        super.setBox(p);
        agentManager = new AtomAgentManager(this,p);
        torqueSum.setAgentManager(agentManager);
        ((MeterKineticEnergyRigid)meterKE).setBox(p);
    }

    public void setOrientationCalc(AtomTypeMolecule moleculeType, OrientationCalc calcer) {
        typeAgentManager.setAgent(moleculeType, calcer);
//        System.out.println("moment: "+typeAgent.moment);
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
            IAtomSet pair = Debug.getAtoms(box);
            if (pair != null) {
                IVector dr = space.makeVector();
                dr.Ev1Mv2(((IAtomPositioned)pair.getAtom(1)).getPosition(), ((IAtomPositioned)pair.getAtom(0)).getPosition());
                System.out.println(pair+" dr "+dr);
            }
        }
        IAtomSet moleculeList = box.getMoleculeList();
        int nMolecules = moleculeList.getAtomCount();
        for (int iMolecule = 0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = (IMolecule)moleculeList.getAtom(iMolecule);
            IAtomSet children = molecule.getChildList();
            OrientationCalc calcer = (OrientationCalc)typeAgentManager.getAgent(molecule.getType());
            if (calcer == null) {
                for (int iLeaf=0; iLeaf<children.getAtomCount(); iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic)children.getAtom(iLeaf);
                    MyAgent agent = (MyAgent)agentManager.getAgent(a);
                    IVector r = a.getPosition();
                    IVector v = a.getVelocity();
                    if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                        System.out.println("first "+a+" r="+r+", v="+v+", f="+agent.force);
                    }
                    v.PEa1Tv1(0.5*timeStep*((AtomTypeLeaf)a.getType()).rm(),agent.force);  // p += f(old)*dt/2
                    r.PEa1Tv1(timeStep,v);         // r += p*dt/m
                }
                continue;
            }
            
            MoleculeAgent agent = (MoleculeAgent)agentManager.getAgent(molecule);
            IAtomOrientedKinetic orientedMolecule = (IAtomOrientedKinetic)molecule;
            IOrientationFull3D orientation = (IOrientationFull3D)orientedMolecule.getOrientation();
            IVector moment = ((AtomTypeMoleculeOriented)molecule.getType()).getMomentOfInertia();

            // use the angular velocity field to store angular momentum during the time step  :(
            IVector angularMomentum = orientedMolecule.getAngularVelocity();
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

            double dtheta = Math.sqrt(tempAngularVelocity.squared());
            tempAngularVelocity.TE(1/dtheta);
            angularMomentum.PEa1Tv1(0.5*timeStep, agent.torque);

            numRigid++;
            for (int i = 0; i<maxIterations; i++) {
                iterationsTotal++;
                tempOrientation.E(orientation);
                // estimate orientation at half timestep, we need this to get the moment of inertia at half timestep
                tempOrientation.rotateBy(0.5*dtheta*timeStep, tempAngularVelocity);

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
                double omegaError = xWork.squared();
                xWork.E(tempAngularVelocity);
                dtheta = Math.sqrt(tempAngularVelocity.squared());
                tempAngularVelocity.TE(1.0/dtheta);

                if (omegaError < omegaTolerance) {
                    break;
                }
                if (i == maxIterations-1) {
                    System.out.println("omegaError still "+omegaError+" after "+maxIterations+" iterations");
//                    throw new RuntimeException("omegaError still "+omegaError+" after "+maxIterations+" iterations");
                }
            }

            orientation.rotateBy(timeStep*dtheta, tempAngularVelocity);
                
//            System.out.println("o "+orientation.getDirection()+" "+orientation.getSecondaryDirection());
            calcer.setOrientation(molecule, orientation);
            //advance linear velocity to half timestep
            orientedMolecule.getVelocity().PEa1Tv1(0.5*timeStep/((AtomTypeMoleculeOriented)molecule.getType()).getMass(), agent.force);
            
            //advance position to full timestep
            IVector transVec = ((AtomActionTranslateBy)translator.getAction()).getTranslationVector();
            transVec.Ea1Tv1(timeStep, orientedMolecule.getVelocity());
            orientedMolecule.getPosition().PE(transVec);
            translator.actionPerformed(molecule);
        }

//        if(forceSum instanceof PotentialCalculationForcePressureSum){
//            pressureTensor.E(((PotentialCalculationForcePressureSum)forceSum).getPressureTensor());
//        }
        
        torqueSum.reset();
        //Compute forces and torques on each atom at t+dt
        potential.calculate(box, allAtoms, torqueSum);
        
        for (int iMolecule = 0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = (IMolecule)moleculeList.getAtom(iMolecule);
            IAtomSet children = molecule.getChildList();
            if ((OrientationCalc)typeAgentManager.getAgent(molecule.getType()) == null) {
                // unimolecular or at least not rigid
                //Finish integration step
                for (int iLeaf=0; iLeaf<children.getAtomCount(); iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic)children.getAtom(iLeaf);
//                    System.out.println("force: "+((MyAgent)a.ia).force.toString());
                    IVector velocity = a.getVelocity();
                    workTensor.Ev1v2(velocity,velocity);
                    workTensor.TE(((AtomTypeLeaf)a.getType()).getMass());
                    pressureTensor.PE(workTensor);
                    if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                        System.out.println("second "+a+" v="+velocity+", f="+((MyAgent)agentManager.getAgent(a)).force);
                    }
                    velocity.PEa1Tv1(0.5*timeStep*((AtomTypeLeaf)a.getType()).rm(),((MyAgent)agentManager.getAgent(a)).force);  //p += f(new)*dt/2
                    currentKineticEnergy += velocity.squared()*((AtomTypeLeaf)a.getType()).getMass();
                }
                // skip the rotational stuff
                continue;
            }
            
            IAtomOrientedKinetic orientedMolecule = (IAtomOrientedKinetic)molecule;
            MoleculeAgent agent = (MoleculeAgent)agentManager.getAgent(molecule);
            //calc torque and linear force
            for (int i=0; i<children.getAtomCount(); i++) {
                IAtomPositioned atom = (IAtomPositioned)children.getAtom(i);
                IVector3D atomForce = (IVector3D)((MyAgent)agentManager.getAgent(atom)).force;
                if (atomForce.isZero()) {
                    continue;
                }

                agent.force.PE(atomForce);

                xWork.Ev1Mv2(atom.getPosition(), orientedMolecule.getPosition());
                xWork.XE(atomForce);
                agent.torque.PE(xWork);
            }
            
            //advance linear velocity to full timestep
            double mass = ((AtomTypeMoleculeOriented)molecule.getType()).getMass();
            orientedMolecule.getVelocity().PEa1Tv1(0.5*timeStep/mass, agent.force);

            //advance momentum to full timestep
            IVector angularVelocity = orientedMolecule.getAngularVelocity();
            IVector moment = ((AtomTypeMoleculeOriented)molecule.getType()).getMomentOfInertia();

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
            xWork.TE(moment);
            currentKineticEnergy += xWork.x(0) + xWork.x(1)+ xWork.x(2);
        }
        pressureTensor.TE(1/box.getBoundary().volume());

        currentKineticEnergy *= 0.5;
        if (printInterval > 0 && stepCount%printInterval == 0) {
            double PE = meterPE.getDataAsScalar();
            int moleculeCount = box.getMoleculeList().getAtomCount();
            double fac = Joule.UNIT.fromSim(1.0/moleculeCount)*Constants.AVOGADRO;
            System.out.println(currentTime+" "+(iterationsTotal/(double)numRigid)+" "+Kelvin.UNIT.fromSim(currentKineticEnergy/moleculeCount/3)+" "
                              +fac*currentKineticEnergy+" "+fac*PE+" "+fac*(PE+currentKineticEnergy));
        }

        if (isothermal) {
            doThermostat();
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
        IAtomSet moleculeList = box.getMoleculeList();
        int nMolecules = moleculeList.getAtomCount();
        int D = 0;
        for (int iMolecule = 0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = (IMolecule)moleculeList.getAtom(iMolecule);
            IAtomSet children = molecule.getChildList();
            if (typeAgentManager.getAgent(molecule.getType()) == null) {
                // unimolecular or at least not rigid
                //Finish integration step
                for (int iLeaf=0; iLeaf<children.getAtomCount(); iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic)children.getAtom(iLeaf);
//                    System.out.println("force: "+((MyAgent)a.ia).force.toString());
                    IVector velocity = a.getVelocity();
                    KE += velocity.squared() * ((AtomTypeLeaf)a.getType()).getMass();
                    D += 3;
                }
                continue;
            }
            IAtomOrientedKinetic orientedMolecule = (IAtomOrientedKinetic)molecule;
            KE += orientedMolecule.getVelocity().squared() * ((AtomTypeMoleculeOriented)molecule.getType()).getMass();

            IVector moment = ((AtomTypeMoleculeOriented)molecule.getType()).getMomentOfInertia();

            tempAngularVelocity.E(orientedMolecule.getAngularVelocity());
            rotationTensor.setOrientation((IOrientationFull3D)orientedMolecule.getOrientation());
            rotationTensor.transform(tempAngularVelocity);
            tempAngularVelocity.DE(moment);

            tempAngularVelocity.TE(tempAngularVelocity);
            tempAngularVelocity.TE(moment);
            KE += tempAngularVelocity.x(0) + tempAngularVelocity.x(1)+ tempAngularVelocity.x(2);
            D += 6;
        }
        double scale = Math.sqrt(temperature*D/KE);
        currentKineticEnergy *= scale*scale;
        for (int iMolecule = 0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = (IMolecule)moleculeList.getAtom(iMolecule);
            IAtomSet children = molecule.getChildList();
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

            IVector angularVelocity = orientedMolecule.getAngularVelocity();
            rotationTensor.setOrientation((IOrientationFull3D)orientedMolecule.getOrientation());
            rotationTensor.transform(angularVelocity);
            angularVelocity.TE(scale);
            rotationTensor.invert();
            rotationTensor.transform(angularVelocity);
        }
    }

    public void randomizeMomenta() {
        IAtomSet moleculeList = box.getMoleculeList();
        int nMolecules = moleculeList.getAtomCount();
//        System.out.println("rerandomize");
        for (int iMolecule = 0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = (IMolecule)moleculeList.getAtom(iMolecule);
            OrientationCalc calcer = (OrientationCalc)typeAgentManager.getAgent(molecule.getType());
            if (calcer == null) {
                IAtomSet children = molecule.getChildList();
                for (int i=0; i<children.getAtomCount(); i++) {
                    super.randomizeMomentum((IAtomKinetic)children.getAtom(i));
                }
                continue;
            }

            MoleculeOrientedDynamic orientedMolecule = (MoleculeOrientedDynamic)molecule;
            IVector velocity = orientedMolecule.getVelocity();
            IVector moment = ((AtomTypeMoleculeOriented)molecule.getType()).getMomentOfInertia();

            double mass = ((AtomTypeMoleculeOriented)molecule.getType()).getMass();
            int D = velocity.getD();
            for(int i=0; i<D; i++) {
                velocity.setX(i,random.nextGaussian());
            }
            velocity.TE(Math.sqrt(temperature/mass));
    
            IVector angularVelocity = orientedMolecule.getAngularVelocity();
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
            
            calcer.calcOrientation(orientedMolecule, orientedMolecule.getOrientation());
            rotationTensor.setOrientation(orientedMolecule.getOrientation());
            // body-fixed to space-fixed, so invert
            rotationTensor.invert();
            // transform to space-fixed angular momentum
            rotationTensor.transform(angularVelocity);
        }
        
    }
    
    public void randomizeMomentum(IAtomKinetic atom) {
        if (atom instanceof AtomLeaf) {
            super.randomizeMomentum(atom);
        }
        else {
//            System.out.println("rerandomize 1");
            OrientationCalc calcer = (OrientationCalc)typeAgentManager.getAgent(atom.getType());
            if (calcer == null) {
                super.randomizeMomentum(atom);
                return;
            }

            MoleculeOrientedDynamic orientedMolecule = (MoleculeOrientedDynamic)atom;
            IVector velocity = orientedMolecule.getVelocity();
            double mass = ((AtomTypeMoleculeOriented)orientedMolecule.getType()).getMass();
            int D = velocity.getD();
            for(int i=0; i<D; i++) {
                velocity.setX(i,random.nextGaussian());
            }
            velocity.TE(Math.sqrt(temperature/mass));

            IVector angularVelocity = orientedMolecule.getAngularVelocity();
            IVector moment = ((AtomTypeMoleculeOriented)orientedMolecule.getType()).getMomentOfInertia();
            for(int i=0; i<D; i++) {
                angularVelocity.setX(i,random.nextGaussian());
            }
            tempAngularVelocity.Ea1Tv1(temperature, moment);
            tempAngularVelocity.map(new Function.Sqrt());
            tempAngularVelocity.TE(angularVelocity);
            angularVelocity.E(tempAngularVelocity);
            //angularVelocity is now the correct body-fixed angular momentum
            angularVelocity.DE(moment);
            
            calcer.calcOrientation(orientedMolecule, orientedMolecule.getOrientation());
            rotationTensor.setOrientation(orientedMolecule.getOrientation());
            // body-fixed to space-fixed, so invert
            rotationTensor.invert();
            // transform to space-fixed angular momentum
            rotationTensor.transform(angularVelocity);
        }
    }
    
    public void reset() throws ConfigurationOverlapException{
        if(!initialized) return;
        
        super.reset();
        if (Debug.ON && Debug.DEBUG_NOW) {
            IAtomSet pair = Debug.getAtoms(box);
            if (pair != null) {
                IVector dr = space.makeVector();
                dr.Ev1Mv2(((IAtomPositioned)pair.getAtom(1)).getPosition(), ((IAtomPositioned)pair.getAtom(0)).getPosition());
                System.out.println(pair+" dr "+dr);
            }
        }
        
        IAtomSet moleculeList = box.getMoleculeList();
        int nMolecules = moleculeList.getAtomCount();

        for (int iMolecule=0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = (IMolecule)moleculeList.getAtom(iMolecule);
            OrientationCalc calcer = (OrientationCalc)typeAgentManager.getAgent(molecule.getType());
            if (calcer == null) {
                continue;
            }

            MoleculeOrientedDynamic orientedMolecule = (MoleculeOrientedDynamic)molecule;
            calcer.calcOrientation(molecule, orientedMolecule.getOrientation());
            orientedMolecule.getPosition().E(atomPositionCOM.position(molecule));
        }

        torqueSum.reset();
        //Compute forces on each atom at t+dt
        potential.calculate(box, allAtoms, torqueSum);

        for (int iMolecule=0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = (IMolecule)moleculeList.getAtom(iMolecule);
            if (typeAgentManager.getAgent(molecule.getType()) == null) {
                continue;
            }

            MoleculeAgent agent = (MoleculeAgent)agentManager.getAgent(molecule);
            MoleculeOrientedDynamic orientedMolecule = (MoleculeOrientedDynamic)molecule;
            //calc angular velocities
            IAtomSet children = molecule.getChildList();
            for (int i=0; i<children.getAtomCount(); i++) {
                IAtomPositioned atom = (IAtomPositioned)children.getAtom(i);
                IVector3D force = (IVector3D)((MyAgent)agentManager.getAgent(atom)).force;
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
    
    public Class getAgentClass() {
        return Object.class;
    }

    public final Object makeAgent(IAtom a) {
        if (a instanceof IMolecule) {
            return new MoleculeAgent(potential.getSpace());
        }
        return new MyAgent(potential.getSpace());
    }
    
    public void releaseAgent(Object agent, IAtom atom) {}
            
    public static class MoleculeAgent implements Integrator.Torquable, Integrator.Forcible, Serializable {  //need public so to use with instanceof
        private static final long serialVersionUID = 1L;
        public final IVector torque;
        public final IVector force;

        public MoleculeAgent(Space space) {
            torque = space.makeVector();
            force = space.makeVector();
        }
        
        public IVector torque() {return torque;}
        public IVector force() {return force;}
    }

    public Class getTypeAgentClass() {
        return OrientationCalc.class;
    }
    
    public Object makeAgent(IAtomType type) {
        return null;
    }

    public void releaseAgent(Object agent, IAtomType type) {}
    
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        boolean periodic = true;
        Box box;
        if (periodic) {
            box = new Box(new BoundaryRectangularPeriodic(sim.getSpace(), sim.getRandom(), 10), space);
        }
        else {
            box = new Box(new BoundaryRectangularNonperiodic(sim.getSpace(), sim.getRandom()), space);
        }
        sim.addBox(box);
        SpeciesWater3POriented species = new SpeciesWater3POriented(sim.getSpace(), true);
        sim.getSpeciesManager().addSpecies(species);
        int numMolecules = 256;
        box.setNMolecules(species, numMolecules);
        box.setDensity(1/18.0*Constants.AVOGADRO/1E24);
        if (true) {
            new ConfigurationLattice(new LatticeCubicFcc(), space).initializeCoordinates(box);
        }
        else {
            new ConfigurationFile("water"+numMolecules+(periodic ? "pbc":"")).initializeCoordinates(box);
        }

        PotentialMaster potentialMaster = new PotentialMaster(sim.getSpace());
        double timeInterval = 0.001;
        int maxIterations = 20;
        IntegratorRigidIterative integrator = new IntegratorRigidIterative(sim, potentialMaster, timeInterval, 1, space);
        integrator.printInterval = 10;
        integrator.setMaxIterations(maxIterations);
        integrator.setBox(box);
        OrientationCalcWater calcer = new OrientationCalcWater(sim.getSpace());
        integrator.setOrientationCalc(species.getMoleculeType(), calcer);
        integrator.setTemperature(Kelvin.UNIT.toSim(298));
        ActivityIntegrate ai = new ActivityIntegrate(integrator);
        sim.getController().addAction(ai);

        P2WaterSPCSoft p2Water = new P2WaterSPCSoft(sim.getSpace());
        
        if (periodic) {
            BoxImposePbc pbc = new BoxImposePbc(box, space);
            pbc.setApplyToMolecules(true);
            integrator.addIntervalAction(pbc);

            double boxlength = box.getBoundary().getDimensions().x(0);

            DipoleSourceWater dipoleSource = new DipoleSourceWater(sim.getSpace());
            dipoleSource.setDipoleStrength(2*Electron.UNIT.toSim(0.41)*Math.cos(109.5/2.0*Math.PI/180));
            P2ReactionFieldDipole pNRF = new P2ReactionFieldDipole(sim.getSpace());
            pNRF.setDipoleSource(dipoleSource);
            pNRF.setRange(boxlength*0.5);
            pNRF.setDielectric(78.4);

            potentialMaster.addPotential(new P2MoleculeSoftTruncatedSwitched(pNRF, boxlength*0.49, space), new ISpecies[]{species, species});
            potentialMaster.lrcMaster().addPotential(pNRF.makeP0(), new PotentialMaster.AtomIterator0(), null);

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
            SimulationGraphic graphic = new SimulationGraphic(sim, "Rigid", 1, space);
            ((ColorSchemeByType)graphic.getDisplayBox(box).getColorScheme()).setColor(species.getHydrogenType(), Color.WHITE);
            ((ColorSchemeByType)graphic.getDisplayBox(box).getColorScheme()).setColor(species.getOxygenType(), Color.RED);
            graphic.makeAndDisplayFrame();
        }
    }
}
