/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.atom.*;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.math.function.Function;
import etomica.molecule.CenterOfMass;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.OrientationCalc;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.space3d.IOrientation3D;
import etomica.space3d.RotationTensor3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesAgentManager;
import etomica.species.SpeciesManager;
import etomica.units.dimensions.Energy;
import etomica.util.Debug;
import etomica.util.random.IRandom;

public class IntegratorRigidIterative extends IntegratorMD implements SpeciesAgentManager.AgentSource<OrientationCalc> {

    protected final RotationTensor3D rotationTensor;
    protected final SpeciesAgentManager<OrientationCalc> typeAgentManager;
    protected final Vector tempAngularVelocity;
    public int printInterval;
    protected int maxIterations;
    protected double omegaTolerance;

    public IntegratorRigidIterative(SpeciesManager sm, IRandom random, PotentialCompute potentialMaster, Box box) {
        this(sm, random, potentialMaster, 0.05, 1.0, box);
    }

    public IntegratorRigidIterative(SpeciesManager sm, IRandom random, PotentialCompute potentialMaster,
                                    double timeStep, double temperature, Box box) {
        super(potentialMaster, random, timeStep, temperature, box);
        rotationTensor = (RotationTensor3D) space.makeRotationTensor();
        typeAgentManager = new SpeciesAgentManager<>(this, sm);
        tempAngularVelocity = space.makeVector();
        maxIterations = 40;
        omegaTolerance = 1.e-30;
        meterKE = new DataSourceScalar("KE", Energy.DIMENSION) {
            @Override
            public double getDataAsScalar() {
                return currentKineticEnergy;
            }
        };
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
        Vector[] torques = potentialCompute.getTorques();
        Vector[] forces = potentialCompute.getForces();
        for (int iMolecule = 0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = moleculeList.get(iMolecule);
            IAtomList children = molecule.getChildList();
            OrientationCalc calcer = typeAgentManager.getAgent(molecule.getType());
            if (calcer == null) {
                for (int iLeaf = 0; iLeaf<children.size(); iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic)children.get(iLeaf);
                    Vector force = forces[a.getLeafIndex()];
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

            Vector com = CenterOfMass.position(box, molecule);

            Vector torque = box.getSpace().makeVector();
            Vector force = box.getSpace().makeVector();
            Vector velocity = box.getSpace().makeVector();
            double mass = 0;
            for (IAtom a : molecule.getChildList()) {
                if (torques != null && torques[a.getLeafIndex()] != null) {
                    torque.PE(torques[a.getLeafIndex()]);
                }
                Vector dr = box.getSpace().makeVector();
                dr.Ev1Mv2(a.getPosition(), com);
                box.getBoundary().nearestImage(dr);
                dr.XE(forces[a.getLeafIndex()]);
                torque.PE(dr);
                force.PE(forces[a.getLeafIndex()]);
                mass += a.getType().getMass();
                velocity.PEa1Tv1(a.getType().getMass(), ((IAtomKinetic)a).getVelocity());
            }
            velocity.TE(1.0/mass);
//            System.out.println("rz "+com.getX(2)+"   vz "+velocity.getX(2)+"  fz "+force.getX(2));

//            IMoleculeOrientedKinetic orientedMolecule = (IMoleculeOrientedKinetic)molecule;
//            IOrientationFull3D orientation = (IOrientationFull3D)orientedMolecule.getOrientation();
            IOrientation3D orientation = (IOrientation3D) calcer.makeOrientation(space);

            calcer.calcOrientation(molecule, orientation);
            Vector angularMomentum = calcer.getAngularMomentum(molecule, com, box);

//            System.out.println("angular velocity "+angularMomentum);
            // transform to body-fixed, multiply by moment of inertia, transform back
            tempAngularVelocity.E(calcer.angularMomentumToVelocity(molecule, space, angularMomentum));

            // xWork = angular velocity from the previous iteration
            // save angular velocity to see if we've converged
            Vector xWork = box.getSpace().makeVector();
            xWork.E(tempAngularVelocity);
//            System.out.println("initial guess half-timestep angular velocity "+angularVelocity);
            angularMomentum.PEa1Tv1(0.5*timeStep, torque);

            double dtheta = Math.sqrt(tempAngularVelocity.squared());
            if (dtheta > 0) {
                tempAngularVelocity.TE(1/dtheta);
            }

            numRigid++;
            for (int i = 0; i<maxIterations; i++) {
                iterationsTotal++;
                IOrientation3D tempOrientation = (IOrientation3D) calcer.makeOrientation(space);
                tempOrientation.E(orientation);
                if (dtheta > 0) {
                    // estimate orientation at half timestep, we need this to get the moment of inertia at half timestep
                    tempOrientation.rotateBy(0.5*dtheta*timeStep, tempAngularVelocity);
                }

                Vector omega = calcer.angularMomentumToVelocity(tempOrientation, molecule, space, angularMomentum);

                xWork.ME(omega);
                double omegaError = xWork.squared() / omega.squared();
                xWork.E(omega);
                dtheta = Math.sqrt(omega.squared());
                if (dtheta > 0) {
                    omega.TE(1.0/dtheta);
                }

                if (omegaError < omegaTolerance) {
                    break;
                }
                if (i == maxIterations-1) {
                    System.err.println("omegaError still "+omegaError+" after "+maxIterations+" iterations");
//                    throw new RuntimeException("omegaError still "+omegaError+" after "+maxIterations+" iterations");
                }
                tempAngularVelocity.E(omega);
            }

            if (dtheta > 0) {
                orientation.rotateBy(timeStep*dtheta, tempAngularVelocity);
            }

//            System.out.println("o "+orientation.getDirection()+" "+orientation.getSecondaryDirection());
            calcer.setOrientation(molecule, box, orientation);

            calcer.setAngularMomentum(molecule, com, box, angularMomentum);
            //advance linear velocity to half timestep
            Vector dv = box.getSpace().makeVector();
            dv.Ea1Tv1(0.5*timeStep/mass, force);
            Vector dr = box.getSpace().makeVector();
            dr.Ea1Tv1(timeStep, velocity);
            dr.PEa1Tv1(timeStep, dv);
            for (IAtom a : molecule.getChildList()) {
                ((IAtomKinetic)a).getVelocity().PE(dv);
                a.getPosition().PE(dr);
            }
        }
        eventManager.forcePrecomputed();

        //Compute forces and torques on each atom at t+dt
        double PE = potentialCompute.computeAll(true);
        torques = potentialCompute.getTorques();
        forces = potentialCompute.getForces();

        eventManager.forceComputed();

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
                    if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(new AtomSetSinglet(a))) {
                        System.out.println("second "+a+" v="+velocity+", f="+forces[a.getLeafIndex()]);
                    }
                    velocity.PEa1Tv1(0.5 * timeStep * a.getType().rm(), forces[a.getLeafIndex()]);  //p += f(new)*dt/2
                    currentKineticEnergy += velocity.squared()* a.getType().getMass();
                }
                // skip the rotational stuff
                continue;
            }

            OrientationCalc calcer = typeAgentManager.getAgent(molecule.getType());
            Vector com = CenterOfMass.position(box, molecule);
            Vector torque = box.getSpace().makeVector();
            Vector force = box.getSpace().makeVector();
            double mass = 0;

            //calc torque and linear force
            for (IAtom a : molecule.getChildList()) {
                if (torques != null && torques[a.getLeafIndex()] != null) {
                    torque.PE(torques[a.getLeafIndex()]);
                }
                Vector dr = box.getSpace().makeVector();
                dr.Ev1Mv2(a.getPosition(), com);
                box.getBoundary().nearestImage(dr);
                dr.XE(forces[a.getLeafIndex()]);
                torque.PE(dr);
                force.PE(forces[a.getLeafIndex()]);
                mass += a.getType().getMass();
            }

            //advance linear velocity to full timestep
            Vector dv = box.getSpace().makeVector();
            dv.Ea1Tv1(0.5*timeStep/mass, force);
            for (IAtom a : molecule.getChildList()) {
                ((IAtomKinetic)a).getVelocity().PE(dv);
            }

            //advance momentum to full timestep
            IOrientation3D orientation = (IOrientation3D) calcer.makeOrientation(space);
            calcer.calcOrientation(molecule, orientation);
            Vector angularMomentum = calcer.getAngularMomentum(molecule, com, box);

            // we actually stored the half-timestep angular momentum in this field...
            // advance to full timestep

            angularMomentum.PEa1Tv1(0.5*timeStep, torque);

            calcer.setAngularMomentum(molecule, com, box, angularMomentum);

            Vector velocity = box.getSpace().makeVector();
            for (IAtom a : molecule.getChildList()) {
                velocity.PEa1Tv1(a.getType().getMass(), ((IAtomKinetic)a).getVelocity());
            }
            velocity.TE(1.0/mass);

            currentKineticEnergy += mass * velocity.squared();
            Vector omega = calcer.angularMomentumToVelocity(molecule, space, angularMomentum);
            omega = calcer.spaceToBody(molecule, space, omega);
            omega.TE(omega);
            currentKineticEnergy += calcer.getMomentOfInertia(molecule).dot(omega);
        }

        currentKineticEnergy *= 0.5;
        if (printInterval > 0 && stepCount%printInterval == 0) {
            int moleculeCount = box.getMoleculeList().size();
            double fac = 1; //Joule.UNIT.fromSim(1.0/moleculeCount)* Constants.AVOGADRO;
            System.out.println(currentTime+" "+(iterationsTotal/(double)numRigid)+" "+(currentKineticEnergy/moleculeCount/5)*2+" "
                              +fac*currentKineticEnergy+" "+fac*PE+" "+fac*(PE+currentKineticEnergy));
        }

        if (isothermal) {
            doThermostatInternal();
        }
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
            for (int iLeaf = 0; iLeaf<children.size(); iLeaf++) {
                IAtomKinetic a = (IAtomKinetic)children.get(iLeaf);
                double mass = a.getType().getMass();
                momentum.PEa1Tv1(mass, a.getVelocity());
                totalMass += mass;
            }
        }
        momentum.TE(1.0/totalMass);
        for (int iMolecule = 0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = moleculeList.get(iMolecule);
            IAtomList children = molecule.getChildList();
            OrientationCalc calcer = typeAgentManager.getAgent(molecule.getType());
            if (calcer == null) {
                D += children.size() * 3;
                for (int iLeaf = 0; iLeaf<children.size(); iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic)children.get(iLeaf);
                    Vector velocity = a.getVelocity();
                    KE += velocity.squared() * a.getType().getMass();
                }
            }
            else {
                for (int iLeaf = 0; iLeaf<children.size(); iLeaf++) {
                    IAtomKinetic a = (IAtomKinetic)children.get(iLeaf);
                    KE += a.getVelocity().squared() * a.getType().getMass();
                }

                D += 3 + calcer.getDOF(space);
            }
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
            for (int iLeaf = 0; iLeaf < children.size(); iLeaf++) {
                IAtomKinetic a = (IAtomKinetic) children.get(iLeaf);
                a.getVelocity().TE(scale);
                if (a instanceof IAtomOrientedKinetic) {
                    ((IAtomOrientedKinetic)a).getAngularVelocity().TE(scale);
                }
            }
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
                IAtomList children = m.getChildList();
                for (int i = 0; i < children.size(); i++) {
                    IAtomKinetic a = (IAtomKinetic)children.get(i);
                    momentum.PEa1Tv1(a.getType().getMass(), a.getVelocity());
                    totalMass += a.getType().getMass();
                }
            }
            if (totalMass == 0) return;
            momentum.TE(1.0 / totalMass);
            //momentum is now net velocity
            //set net momentum to 0
            for (int iLeaf = 0; iLeaf < nMolecules; iLeaf++) {
                IMolecule m = moleculeList.get(iLeaf);
                IAtomList children = m.getChildList();
                for (int i = 0; i < children.size(); i++) {
                    IAtomKinetic a = (IAtomKinetic)children.get(i);
                    double rm = a.getType().rm();
                    if (rm != 0 && rm != Double.POSITIVE_INFINITY) {
                        a.getVelocity().ME(momentum);
                    }
                }
            }
            if (Debug.ON) {
                momentum.E(0);
                for (int iLeaf = 0; iLeaf < nMolecules; iLeaf++) {
                    IMolecule m = moleculeList.get(iLeaf);
                    IAtomList children = m.getChildList();
                    for (int i = 0; i < children.size(); i++) {
                        IAtomKinetic a = (IAtomKinetic)children.get(i);
                        momentum.PEa1Tv1(a.getType().getMass(), a.getVelocity());
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
        currentKineticEnergy = 0;
        for (int iMolecule = 0; iMolecule<nMolecules; iMolecule++) {
            IMolecule molecule = moleculeList.get(iMolecule);
            OrientationCalc calcer = typeAgentManager.getAgent(molecule.getType());
            if (calcer == null) {
                IAtomList children = molecule.getChildList();
                for (int i = 0; i<children.size(); i++) {
                    IAtomKinetic a = (IAtomKinetic) children.get(i);
                    super.randomizeMomentum(a);
                    currentKineticEnergy += 0.5*a.getType().getMass()*a.getVelocity().squared();
                }
                continue;
            }

            Vector com = CenterOfMass.position(box, molecule);
            double mass = 0;
            for (IAtom a : molecule.getChildList()) {
                Vector dr = box.getSpace().makeVector();
                dr.Ev1Mv2(a.getPosition(), com);
                box.getBoundary().nearestImage(dr);
                mass += a.getType().getMass();
            }

            int D = box.getSpace().getD();
            Vector velocity = box.getSpace().makeVector();
            for(int i=0; i<D; i++) {
                velocity.setX(i,random.nextGaussian());
            }
            velocity.TE(Math.sqrt(temperature/mass));
            for (IAtom a : molecule.getChildList()) {
                ((IAtomKinetic)a).getVelocity().E(velocity);
            }
            currentKineticEnergy += 0.5*mass*velocity.squared();

            Vector angularVelocity = box.getSpace().makeVector();
            for(int i=0; i<D; i++) {
                angularVelocity.setX(i,random.nextGaussian());
            }

            tempAngularVelocity.Ea1Tv1(temperature, calcer.getMomentOfInertia(molecule));
            tempAngularVelocity.map(new Function.Sqrt());
            tempAngularVelocity.TE(angularVelocity);
            Vector tmp = box.getSpace().makeVector();
            tmp.E(calcer.angularMomentumToVelocity(calcer.makeOrientation(space), molecule, space, tempAngularVelocity));
            tmp.TE(tmp);
            currentKineticEnergy += 0.5*calcer.getMomentOfInertia(molecule).dot(tmp);
            //tempAngularVelocity is now the correct body-fixed angular momentum
            Vector L = calcer.bodyToSpace(molecule, space, tempAngularVelocity);
            calcer.setAngularMomentum(molecule, com, box, L);
        }
    }

    public void reset() {
        super.reset();

        if (Debug.ON && Debug.DEBUG_NOW) {
            IAtomList pair = Debug.getAtoms(box);
            if (pair != null) {
                Vector dr = space.makeVector();
                dr.Ev1Mv2(pair.get(1).getPosition(), pair.get(0).getPosition());
                System.out.println(pair + " dr " + dr);
            }
        }

        eventManager.forcePrecomputed();

        currentPotentialEnergy = potentialCompute.computeAll(true);

        eventManager.forceComputed();

        if (printInterval == 1) {
            double PE = currentPotentialEnergy;
            int moleculeCount = box.getMoleculeList().size();
            double fac = 1; //Joule.UNIT.fromSim(1.0/moleculeCount)*Constants.AVOGADRO;
            System.out.println(currentTime+" "+0+" "+(currentKineticEnergy/moleculeCount/5)*2+" "
                    +fac*currentKineticEnergy+" "+fac*PE+" "+fac*(PE+currentKineticEnergy));
        }

    }

//--------------------------------------------------------------

    public OrientationCalc makeAgent(ISpecies type) {
        return null;
    }

    public void releaseAgent(OrientationCalc agent, ISpecies type) {
    }

//    public static void main(String[] args) {
//        Space space = Space3D.getInstance();
//        Simulation sim = new Simulation(space);
//        SpeciesGeneral species = SpeciesWater3P.create(true, true);
//        sim.addSpecies(species);
//        boolean periodic = true;
//        Box box;
//        if (periodic) {
//            box = new Box(new BoundaryRectangularPeriodic(sim.getSpace(), 10), space);
//        }
//        else {
//            box = new Box(new BoundaryRectangularNonperiodic(sim.getSpace()), space);
//        }
//        sim.addBox(box);
//        int numMolecules = 256;
//        box.setNMolecules(species, numMolecules);
//        box.setDensity(1/18.0*Constants.AVOGADRO/1E24);
//        if (true) {
//            new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
//        }
//        else {
//            new ConfigurationFile("water"+numMolecules+(periodic ? "pbc":"")).initializeCoordinates(box);
//        }
//
//        PotentialMasterFasterer potentialMaster = new PotentialMasterFasterer(sim.getSpeciesManager(), box, BondingInfo.noBonding());
//        double timeInterval = 0.001;
//        int maxIterations = 20;
//        IntegratorRigidIterativeFasterer integrator = new IntegratorRigidIterativeFasterer(sim.getSpeciesManager(), sim.getRandom(), potentialMaster, timeInterval, 1, box);
//        integrator.printInterval = 10;
//        integrator.setMaxIterations(maxIterations);
//        OrientationCalcWater3P calcer = new OrientationCalcWater3P(sim.getSpace());
//        integrator.setOrientationCalc(species, calcer);
//        integrator.setTemperature(Kelvin.UNIT.toSim(298));
//
//        P2WaterSPCSoft p2Water = new P2WaterSPCSoft(sim.getSpace());
//
//        if (periodic) {
//            BoxImposePbc pbc = new BoxImposePbc(box, space);
//            pbc.setApplyToMolecules(true);
//            integrator.getEventManager().addListener(new IntegratorListenerAction(pbc));
//
//            double boxlength = box.getBoundary().getBoxSize().getX(0);
//
//            DipoleSourceMolecularWater dipoleSource = new DipoleSourceMolecularWater(sim.getSpace());
//            dipoleSource.setDipoleStrength(2*Electron.UNIT.toSim(0.41)*Math.cos(109.5/2.0*Math.PI/180));
//            IMoleculePositionDefinition positionDefinition = new MoleculePositionCOM(space) ;
//            P2ReactionFieldDipole pNRF = new P2ReactionFieldDipole(sim.getSpace(),positionDefinition);
//            pNRF.setDipoleSource(dipoleSource);
//            pNRF.setRange(boxlength*0.5);
//            pNRF.setDielectric(78.4);
//
//            potentialMaster.addPotential(new P2MoleculeSoftTruncatedSwitched(pNRF, boxlength*0.49, space), new ISpecies[]{species, species});
//            potentialMaster.lrcMaster().addPotential(pNRF.makeP0());
//
//            potentialMaster.addPotential(new P2MoleculeSoftTruncatedSwitched(p2Water, boxlength*0.49, space), new ISpecies[]{species,species});
//        }
//        else {
//            potentialMaster.addPotential(p2Water, new ISpecies[]{species,species});
//        }
//
//        if (true) {
//            sim.getController().runActivityBlocking(new ActivityIntegrate(integrator, 1000000));
//        }
//        else {
//            sim.getController().setSleepPeriod(2);
//            sim.getController().addActivity(new ActivityIntegrate(integrator));
//            SimulationGraphic graphic = new SimulationGraphic(sim, "Rigid", 1);
//            ((ColorSchemeByType) graphic.getDisplayBox(box).getColorScheme()).setColor(species.getTypeByName("H"), Color.WHITE);
//            ((ColorSchemeByType) graphic.getDisplayBox(box).getColorScheme()).setColor(species.getTypeByName("O"), Color.RED);
//            graphic.makeAndDisplayFrame();
//        }
//    }
}
