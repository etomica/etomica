/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.*;
import etomica.box.Box;
import etomica.config.ConformationLinear;
import etomica.data.DataSourceScalar;
import etomica.graphics.*;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorLangevin;
import etomica.integrator.IntegratorMD;
import etomica.potential.P1Anharmonic234;
import etomica.potential.P1AnharmonicTIA;
import etomica.potential.P2Harmonic;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeField;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space2d.Space2D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.units.dimensions.Length;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;

public class SimQuantumAOPIMDGraphic extends Simulation {

    public PotentialComputeField pcP1, pcP2, pcP1EnTIA;
    public final IntegratorMD integrator;
    public final PotentialMasterBonding pmBonding;
    public final Box box;
    public final PotentialComputeAggregate pmAgg;
    public P1Anharmonic234 p1ah, p2ah;
    public P1AnharmonicTIA p1ahUeff;
    public double betaN;
    public double k2_kin;
    public int nBeads;
    public MCMoveHOReal2 moveStageSimple, moveStageEC;
    public int dim;

    public SimQuantumAOPIMDGraphic(Space space, MoveChoice coordType, double mass, double timeStep, double gammaLangevin, int nBeads, double temperature, double k3, double k4, double omega, boolean isTIA, double hbar, boolean isExactA, boolean isCayleyA) {
        super(space);

        if (isExactA && isCayleyA) {
            System.out.println(" Can not have both isExactA and isCayleyA to be true!");
            System.exit(0);
        }
        SpeciesGeneral species = new SpeciesBuilder(space)
                .setDynamic(true)
                .addCount(AtomType.simple("A", mass / nBeads), nBeads)
                .withConformation(new ConformationLinear(space, 0))
                .build();
        addSpecies(species);
        SpeciesGeneral speciesLattice = null;
        if (space.D() > 1) {
            speciesLattice = new SpeciesBuilder(space)
                    .setDynamic(true)
                    .addAtom(AtomType.simple("L", Double.POSITIVE_INFINITY), space.makeVector())
                    .build();
            addSpecies(speciesLattice);
        }
        this.nBeads = nBeads;
        box = this.makeBox(new BoundaryRectangularNonperiodic(space));
        box.setNMolecules(species, 1);
        if (speciesLattice!=null) box.setNMolecules(speciesLattice, 1);
        //pm2 that uses the full PI potential, for data collection
        //spring P2 part (x_i-x_{i+1})^2
        pmBonding = new PotentialMasterBonding(getSpeciesManager(), box);
        this.dim = space.D();
        double beta = 1.0/temperature;
        betaN = beta/nBeads;
        double omegaN = nBeads/(hbar*beta);
        double omegaN2 = omegaN*omegaN;
        double omega2 = omega*omega;

        k2_kin = nBeads == 1 ? 0 : mass*omegaN2/nBeads;

        P2Harmonic p2Bond = new P2Harmonic(k2_kin, 0);
        List<int[]> pairs = new ArrayList<>();
        for (int i=0; i<nBeads; i++) {
            int[] p = new int[]{i,(i+1)%nBeads};
            pairs.add(p);
        }
        pmBonding.setBondingPotentialPair(species, p2Bond, pairs);

        pcP1 = new PotentialComputeField(getSpeciesManager(), box);
        pcP2 = new PotentialComputeField(getSpeciesManager(), box);

        if (isTIA){
            double facUeff = 1.0;
            p1ahUeff = new P1AnharmonicTIA(space, 1, k4, nBeads, mass*omegaN2/nBeads, facUeff);
            pcP1.setFieldPotential(species.getLeafType(), p1ahUeff);
        } else {
            p1ah = new P1Anharmonic234(space, omega2/nBeads, k3/nBeads, k4/nBeads);
            pcP1.setFieldPotential(species.getLeafType(), p1ah);
            p2ah = new P1Anharmonic234(space, 0,k3/nBeads, k4/nBeads);
            pcP2.setFieldPotential(species.getLeafType(), p2ah);
        }

        PotentialComputeAggregate.localStorageDefault = true;
        pmAgg = new PotentialComputeAggregate(pmBonding, pcP1);

        double facEn = 3.0;
        P1AnharmonicTIA p1ahEn = new P1AnharmonicTIA(space, 1, k4, nBeads, mass*omegaN2/nBeads, facEn);
        pcP1EnTIA = new PotentialComputeField(getSpeciesManager(), box);
        pcP1EnTIA.setFieldPotential(species.getLeafType(), p1ahEn);

        if (coordType == MoveChoice.Real) {
            integrator = new IntegratorLangevin(pmAgg, random, timeStep, temperature, box, gammaLangevin);
        } else if (coordType == MoveChoice.NM) {
            MCMoveHO move = new MCMoveHO(space, pmAgg, random, temperature, 0, box, hbar);
            if (isExactA) {
                integrator = new IntegratorLangevinPINMExactA(pcP1, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            } else if (isCayleyA) {
                integrator = new IntegratorLangevinPINMCayleyA(pcP1, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            } else {
                integrator = new IntegratorLangevinPINM(pmAgg, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            }
        } else if (coordType == MoveChoice.NMEC) {
            MCMoveHO move = new MCMoveHO(space, pmAgg, random, temperature, omega2, box, hbar);
            if (isExactA) {
                integrator = new IntegratorLangevinPINMExactA(pcP2, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            } else if (isCayleyA) {
                integrator = new IntegratorLangevinPINMCayleyA(pcP2, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            } else {
                integrator = new IntegratorLangevinPINM(pmAgg, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            }
        } else if (coordType == MoveChoice.Stage) {
            MCMoveHOReal2 move = new MCMoveHOReal2(space, pmAgg, random, temperature, 0, box, hbar);
            if (isExactA) {
                integrator = new IntegratorLangevinPIExactA(pcP1, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            } else if (isCayleyA) {
                integrator = new IntegratorLangevinPICayleyA(pcP1, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            } else {
                integrator = new IntegratorLangevinPI(pmAgg, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            }
        } else { //StageEC -- default
            MCMoveHOReal2 move = new MCMoveHOReal2(space, pmAgg, random, temperature, omega2, box, hbar);
            if (isExactA) {
                integrator = new IntegratorLangevinPIExactA(pcP2, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            } else if (isCayleyA) {
                integrator = new IntegratorLangevinPICayleyA(pcP2, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            } else {
                integrator = new IntegratorLangevinPI(pmAgg, random, timeStep, temperature, box, gammaLangevin, move, hbar, omega2);
            }
        }

        moveStageSimple = new MCMoveHOReal2(space, pmAgg, random, temperature, 0, box, hbar);
        moveStageEC = new MCMoveHOReal2(space, pmAgg, random, temperature, omega2, box, hbar);
        integrator.setThermostatNoDrift(false);
        integrator.setIsothermal(!true);
    }

    public Integrator getIntegrator() {
        return integrator;
    }

    /**
     * This draws a line for a bond in 2D
     */
    public static class MyBond implements Drawable {
        public final Vector v1, v2;
        public final Color c;
        public final Boundary b;
        public MyBond(Vector v1, Vector v2, Color c, Boundary b) {
            this.v1 = v1;
            this.v2 = v2;
            this.c = c;
            this.b = b;
        }
        public void draw(Graphics g, int[] origin, double toPixels) {
            int x1 = origin[0] + (int)(0.5*toPixels*b.getBoxSize().getX(0)) + (int)(toPixels*v1.getX(0));
            int y1 = origin[1] + (int)(0.5*toPixels*b.getBoxSize().getX(1)) + (int)(toPixels*v1.getX(1));
            int x2 = origin[0] + (int)(0.5*toPixels*b.getBoxSize().getX(0)) + (int)(toPixels*v2.getX(0));
            int y2 = origin[1] + (int)(0.5*toPixels*b.getBoxSize().getX(1)) + (int)(toPixels*v2.getX(1));
            g.setColor(c);
            g.drawLine(x1,y1,x2,y2);
        }
    }


    public static void main(String[] args) {
        final long startTime = System.currentTimeMillis();
        OctaneParams params = new OctaneParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            // custom parameters
            params.hbar = 1.0;
            params.steps = 100000;
            params.temperature = 1;
            params.omega = 1;
            params.k3 = 0.;
            params.k4 = 0.;
            params.onlyCentroid = false;
            params.gammaLangevin = params.omega;

//            params.coordType = MoveChoice.Real;
//            params.coordType = MoveChoice.NM;
//            params.coordType = MoveChoice.NMEC;
//            params.coordType = MoveChoice.Stage;
            params.coordType = MoveChoice.StageEC;
            params.timeStep = 0.00001;
            params.nBeads = 3;
            params.isGraphic = !false;


        }

        int nShifts = params.nShifts;
        double mass = params.mass;
        double temperature = params.temperature;
        double hbar = params.hbar;
        double omega = params.omega;
        double k3 = params.k3;
        double k4 = params.k4;
        double gammaLangevin = params.gammaLangevin;
        boolean isGraphic = params.isGraphic;
        boolean isExactA = params.isExactA;
        boolean isCayley = params.isCayleyA;
        long steps = params.steps;
        long stepsEq = steps/10;
        boolean isTIA = params.isTIA;
        boolean zerok0 = params.zerok0;
        boolean onlyCentroid = params.onlyCentroid;
        MoveChoice coordType = params.coordType;
        double omega2 = omega*omega;

        double x = 1/temperature*hbar*omega;
        int nBeads = params.nBeads;
        if (nBeads == -1){
            nBeads = (int) (20*x);
        }

        double omegaN = nBeads*temperature/hbar;
        double timeStep = params.timeStep;

        if (isTIA){
            omega2 = omega2*(1.0 + omega2/12.0/(omegaN*omegaN));
        }
//        if (zerok0) omega2 = 0;

        final SimQuantumAOPIMDGraphic sim = new SimQuantumAOPIMDGraphic(Space2D.getInstance(), coordType, mass, timeStep, gammaLangevin, nBeads, temperature, k3, k4, omega, isTIA, hbar, isExactA, isCayley);
        sim.integrator.reset();

        System.out.println(" PIMD-" + coordType);
        System.out.println(" mass: " + mass);
        System.out.println(" T: " + temperature);
        System.out.println(" hbar: " + hbar);
        System.out.println(" w: " + Math.sqrt(omega2));
        System.out.println(" wn: " + omegaN  + " , w/sqrt(n): " + Math.sqrt(omega2)/Math.sqrt(nBeads));
        System.out.println(" x = beta*hbar*w = " + hbar*omega/temperature);
        System.out.println(" nBeads: " + nBeads);
        System.out.println(" nShifts: "+ nShifts);
        System.out.println(" steps: " +  steps + " stepsEq: " + stepsEq);
        System.out.println(" timestep: " + timeStep);
        System.out.println(" k3: " + k3);
        System.out.println(" k4: " + k4);
        System.out.println(" isTIA: " + isTIA);
        System.out.println(" gammaLangevin: " + gammaLangevin);
        System.out.println(" onlyCentroid: " + onlyCentroid);
        System.out.println(" isExactA: " + isExactA);
        System.out.println(" isCayley: " + isCayley);

        System.out.println("\n Quantum Harmonic Oscillator Theory");
        System.out.println(" ====================================");
        double alpha = 1 + 0.5*Math.pow(hbar*sim.betaN*omega,2)+0.5*hbar* sim.betaN*omega*Math.sqrt(4+Math.pow(hbar* sim.betaN*omega,2));
        double alpha2 = alpha*alpha;
        double hbar2 = hbar*hbar;
        double dAlphaDBeta = 2.0/temperature*hbar2*omega2/nBeads/nBeads*alpha2/(alpha2-1);
        double dAlphadT = -1.0/temperature/temperature*dAlphaDBeta;
        double EnQ = sim.space.D()*(hbar2*omega2)*sim.betaN*alpha/(alpha*alpha-1)*(Math.pow(alpha,nBeads)+1)/(Math.pow(alpha,nBeads)-1);
        double numerator = 1 + alpha2 - Math.pow(alpha,2*nBeads)*(alpha2+1)-2*nBeads*(alpha2-1)*Math.pow(alpha,nBeads);
        double denominator = (alpha2-1)*(alpha2-1)*(Math.pow(alpha,nBeads)-1)*(Math.pow(alpha,nBeads)-1);
        double CvnQ = sim.space.D()*hbar2*omega2*sim.betaN*dAlphadT*numerator/denominator-1/temperature/temperature*EnQ*temperature;
        double EnQinf = sim.space.D()*hbar*omega*(0.5 + 1/(Math.exp(nBeads*sim.betaN*hbar*omega)-1.0));
        double CvnQinf = sim.space.D()*Math.pow(1.0/temperature*hbar*omega/2/Math.sinh(1.0/temperature*hbar*omega/2), 2);
        double EnC = sim.space.D()*temperature;
        double CvnC = sim.space.D();
        System.out.println(" En_ho_c: " + EnC);
        System.out.println(" En_ho_q: " + EnQ);
        System.out.println(" E_ho_q: " + EnQinf);
        System.out.println(" Cvn_ho_c: " + CvnC);
        System.out.println(" Cvn_ho_q: " + CvnQ);
        System.out.println(" Cv_ho_q: " + CvnQinf + "\n");


//        System.out.println("Real: " + nBeads/omegaN/Math.sqrt(nBeads));
//        double s2 = omega2*nBeads*nBeads/omegaN/omegaN*(1 + 1.0/12.0*(nBeads*nBeads-1.0)/nBeads);
//        System.out.println("SS: " + Math.sqrt(s2)/omega);
//        s2 = omega2;
//        System.out.println("SNM: " + Math.sqrt(s2)/omega);
//        System.out.println("EC-bassed: " + 1/omega);
//        System.exit(0);



        if (isGraphic) {
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.GRAPHIC_ONLY);
            int intervalG = 100;
            simGraphic.setPaintInterval(sim.box, intervalG);
            int finalNBeads = nBeads;
            ColorScheme colorScheme = new ColorScheme() {
                protected Color[] allColors;

                public Color getAtomColor(IAtom a) {
                    if (a.getType().getMass() == Double.POSITIVE_INFINITY) return sim.space.D() == 2 ? Color.BLACK : Color.WHITE;
                    if (allColors == null) {
                        allColors = new Color[768];
                        for (int i = 0; i < 256; i++) {
                            allColors[i] = new Color(255 - i, i, 0);
                        }
                        for (int i = 0; i < 256; i++) {
                            allColors[i + 256] = new Color(0, 255 - i, i);
                        }
                        for (int i = 0; i < 256; i++) {
                            allColors[i + 512] = new Color(i, 0, 255 - i);
                        }
                    }
                    return allColors[(768 * a.getIndex() / (finalNBeads))];
                }
            };

            DisplayBox displayBox = simGraphic.getDisplayBox(sim.box);
            displayBox.setColorScheme(colorScheme);
            ((DiameterHashByType) displayBox.getDiameterHash()).setDiameter(sim.getSpecies(0).getAtomType(0), 0.6);
            ((DiameterHashByType) displayBox.getDiameterHash()).setDiameter(sim.getSpecies(1).getAtomType(0), 0.1);

            if (sim.space.D() == 3) {
                AtomPair pair = new AtomPair();
                for (int j = 0; j < 1; j++) {
                    IAtomList beads = sim.box.getMoleculeList().get(j).getChildList();
                    for (int i = 0; i < nBeads; i++) {
                        pair.atom0 = beads.get(i);
                        int next = i + 1;
                        if (next == nBeads) next = 0;
                        pair.atom1 = beads.get(next);
                        ((DisplayBoxCanvasG3DSys) displayBox.canvas).makeBond(pair, null);
                    }
                }
                IAtomList beads = sim.box.getLeafList();
                for (int i = 0; i < nBeads; i++) {
                    pair.atom0 = beads.get(i);
                    pair.atom1 = beads.get(nBeads);
                    ((DisplayBoxCanvasG3DSys) displayBox.canvas).makeBond(pair, Color.BLUE);
                }
            }
            else if (sim.space.D() == 2) {
                for (int j = 0; j < 1; j++) {
                    IAtomList beads = sim.box.getMoleculeList().get(j).getChildList();
                    for (int i = 0; i < nBeads; i++) {
                        int next = i + 1;
                        if (next == nBeads) next = 0;
                        displayBox.addDrawable(new MyBond(beads.get(i).getPosition(), beads.get(next).getPosition(), Color.RED, sim.box.getBoundary()));
                    }
                }
                IAtomList beads = sim.box.getLeafList();
                for (int i = 0; i < nBeads; i++) {
                    displayBox.addDrawable(new MyBond(beads.get(i).getPosition(), beads.get(nBeads).getPosition(), Color.BLUE, sim.box.getBoundary()));
                }
            }

            simGraphic.makeAndDisplayFrame("PIMD - "+coordType);

            DataSourceScalar meterCOM = new DataSourceScalar("COM", Length.DIMENSION) {
                @Override
                public double getDataAsScalar() {
                    Vector COM = sim.box.getSpace().makeVector();
                    for (IAtom atom : sim.box.getLeafList()) {
                        COM.PE(atom.getPosition());
                    }
                    COM.TE(1.0/sim.box.getLeafList().size());
                    return COM.getX(0);
                }
            };

            return;
        }

        System.out.flush();
        // equilibration
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, stepsEq));

        System.out.println("\n equilibration finished");
        int interval = 1;
        int blocks = 100;
        long blockSize = steps / (interval * blocks);
        if (blockSize == 0) blockSize = 1;
        System.out.println(" numBlocks: " + blocks + " blocksize: " + blockSize + " interval: " + interval);

    }

    public enum MoveChoice {Real, NM, NMEC, Stage, StageEC};

    public static class OctaneParams extends ParameterBase {
        public double temperature = 1;
        public double hbar = 1;
        public double omega = 1;
        public double gammaLangevin = omega;
        public double k3 = 0.1;
        public double k4 = 0.01;
        public long steps = 10_000_000;
        public boolean isGraphic = false;
        public boolean isTIA = false;
        public boolean zerok0 = false;
        public boolean onlyCentroid = true;
        public double mass = 1.0;
        public MoveChoice coordType = MoveChoice.Real;
        public double timeStep = -1;
        public int nBeads = -1;
        public int nShifts = 0;
        public boolean isExactA = false;
        public boolean isCayleyA = false;
    }
}