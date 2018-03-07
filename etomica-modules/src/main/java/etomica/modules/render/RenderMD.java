/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.render;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.*;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeCubicFcc;
import etomica.integrator.IntegratorListenerAction;
import etomica.modules.render.ParseObj.BondInfo;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2PenetrableSquareWell;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.random.RandomNumberGenerator;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * 
 * Three-dimensional hard-sphere molecular dynamics simulation, using
 * neighbor listing.  
 * <p>
 * Developed as a prototype and example for the construction of a basic simulation.
 *
 * @author David Kofke and Andrew Schultz
 *
 */
public class RenderMD extends Simulation {

    //the following fields are made accessible for convenience to permit simple
    //mutation of the default behavior

    private static final long serialVersionUID = 1L;
    /**
     * The Box holding the atoms. 
     */
    public final Box box;
    /**
     * The Integrator performing the dynamics.
     */
    public final IntegratorHard integrator;
    /**
     * The single hard-sphere species.
     */
    public final SpeciesSpheresMono species;

    public final P2PenetrableCar potentialBonded;
    
    public final PotentialMasterList potentialMaster;
    
    public final ParseObj parser;
    public ActivityIntegrate activityIntegrate;
    public CriterionCar criterion;

    
    /**
     * Sole public constructor, makes a simulation using a 3D space.
     */
    public RenderMD(Space _space) {
        this(_space, new RenderMDParam());
    }
    
    public RenderMD(Space _space, RenderMDParam params) {

        // invoke the superclass constructor
        // "true" is indicating to the superclass that this is a dynamic simulation
        // the PotentialMaster is selected such as to implement neighbor listing
        super(_space);
        setRandom(new RandomNumberGenerator(2));

        potentialMaster = new PotentialMasterList(this, 1, space);

        parser = new ParseObj(params.file);

        int numAtoms = parser.nAtoms;
        box = this.makeBox();
        integrator = new IntegratorHard(this, potentialMaster, box);
        integrator.setTemperature(params.temperature);
        integrator.setIsothermal(true);
        integrator.setTimeStep(params.timeStep);

        activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(0);
        getController().addAction(activityIntegrate);

        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);
        box.setNMolecules(species, numAtoms);

        potentialBonded = new P2PenetrableCar(space, parser, box);

        potentialBonded.setEpsilonCore(params.epsilonCore);
        potentialBonded.setEpsilon(params.epsilon);
        potentialBonded.setLambda(params.lambda);

        AtomType leafType = species.getLeafType();

        potentialMaster.addPotential(potentialBonded, new AtomType[]{leafType, leafType});

        IAtomList leafList = box.getLeafList();
        for (int iLeaf = 0; iLeaf < numAtoms; iLeaf++) {
            IAtom a = leafList.get(iLeaf);
            Vector pos = a.getPosition();
            pos.E(parser.vertices.get(iLeaf));
        }

        criterion = new CriterionCar(parser, box);
        potentialMaster.setCriterion(potentialBonded, criterion);
        integrator.setThermostat(params.thermostatType);

//        int bondSum = 0;
//        int bondMax = 0;
//        int bondMin = 0;
//        for (int iLeaf=0; iLeaf<numAtoms; iLeaf++) {
//            IAtom a = leafList.getAtom(iLeaf);
//            int bondCount = potential.getBondedList(a).size();
//            bondSum += bondCount;
//            bondMax = Math.max(bondMax, bondCount);
//            bondMin = Math.min(bondMin, bondCount);
//        }
//        System.out.println("Avg, max, min bonds per atom: "+((float)bondSum/numAtoms)+" "+bondMax+" "+bondMin);
        BoxImposePbc imposepbc = new BoxImposePbc(space);
        imposepbc.setBox(box);
        integrator.getEventManager().addListener(new IntegratorListenerAction(imposepbc));

        double vNew = box.getMoleculeList().size() / params.density;
        double scale = Math.pow(vNew / box.getBoundary().volume(), 1.0 / 3.0);

        if (params.initCar) {
            Vector3D dimVector = new Vector3D();
            dimVector.E(box.getBoundary().getBoxSize());
            dimVector.TE(scale);
            box.getBoundary().setBoxSize(dimVector);
        } else {
            BoxInflate inflater = new BoxInflate(box, space);
            inflater.setTargetDensity(params.density);
            inflater.actionPerformed();
            new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        }
        potentialMaster.setRange(box.getBoundary().getBoxSize().getX(0) * 0.49);

        // find neighbors now.  don't try to update later (neighbors never change)
        potentialMaster.getNeighborManager(box).reset();
    }

    public static RenderMDParam getParameters() {
        return new RenderMDParam();
    }
    
    public static class CriterionCar implements NeighborCriterion {

        protected final Map<IAtom,Set<IAtom>> bondedSet;
        public CriterionCar(ParseObj parser, Box box) {
            bondedSet = new HashMap<IAtom,Set<IAtom>>();
            IAtomList leafList = box.getLeafList();
            for (int i = 0; i<leafList.size(); i++) {
                bondedSet.put(leafList.get(i), new HashSet<IAtom>());
            }
            int nBonds = parser.bondList.size();
            for(int i=0; i<nBonds; i++) {
                BondInfo bond = parser.bondList.get(i);
                IAtom atom0 = leafList.get(bond.i0);
                IAtom atom1 = leafList.get(bond.i1);
                bondedSet.get(atom0).add(atom1);
                bondedSet.get(atom1).add(atom0);
            }
        }


        public boolean accept(IAtomList pair) {
            return bondedSet.get(pair.get(0)).contains(pair.get(1));
        }

        public boolean needUpdate(IAtom atom) {return false;}

        public void setBox(Box box) {}

        public boolean unsafe() {return false;}

        public void reset(IAtom atom) {}
    }

    public static class P2PenetrableCar extends P2PenetrableSquareWell {
        protected final Map<IAtomList, Double> bondMap;
        private final double epsMult = 1.0;

        public P2PenetrableCar(Space space, ParseObj parser, Box box) {
            super(space);
            bondMap = new HashMap<IAtomList, Double>();
            IAtomList leafList = box.getLeafList();
            int nBonds = parser.bondList.size();
            for(int i=0; i<nBonds; i++) {
                BondInfo bond = parser.bondList.get(i);
                IAtom atom0 = leafList.get(bond.i0);
                IAtom atom1 = leafList.get(bond.i1);
                bondMap.put(new AtomPair(atom0, atom1), bond.bondLengthSquared*0.9999);
            }
        }

        public void bump(IAtomList pair, double falseTime) {
            IAtomKinetic atom0 = (IAtomKinetic)pair.get(0);
            IAtomKinetic atom1 = (IAtomKinetic)pair.get(1);
            double v2old = atom1.getVelocity().Mv1Squared(atom0.getVelocity());

            setCoreDiameterSquared(bondMap.get(pair));
            super.bump(pair, falseTime);

            double v2new = atom1.getVelocity().Mv1Squared(atom0.getVelocity());
            if(v2new > 1.01*v2old) {
//                atom0.getVelocity().TE(0.1);
//                atom1.getVelocity().TE(0.1);
            }

        }

        public double collisionTime(IAtomList pair, double falseTime) {
            setCoreDiameterSquared(bondMap.get(pair));
            return super.collisionTime(pair, falseTime);
        }

        public double getEpsilon() {
            return super.getEpsilon()/epsMult;
        }

        public void setEpsilon(double eps) {
            super.setEpsilon(epsMult * eps);
        }

    }

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class RenderMDParam extends ParameterBase {
        public String file = "/Users/kofke/Documents/workspace/car self-assembly/mustang.txt";
        public double epsilonCore = 0.0;
        public double lambda = 1.01;
        public double epsilon = 1;
        public double temperature = 0.1;
        public double density = 8.0;
        public double timeStep = 0.02;
        public boolean initCar = true;
        public boolean drawBonds = true;
        public ThermostatType thermostatType = ThermostatType.ANDERSEN;
        
    }
}
