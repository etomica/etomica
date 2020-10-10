/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.clathrates;

import etomica.atom.*;
import etomica.box.Box;
import etomica.box.storage.DoubleStorage;
import etomica.box.storage.Tokens;
import etomica.box.storage.VectorStorage;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationFileBinary;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveOrthorhombic;
import etomica.models.clathrates.MinimizationTIP4P.ChargeAgentSourceRPM;
import etomica.models.water.ConfigurationFileTIP4P;
import etomica.models.water.SpeciesWater4P;
import etomica.molecule.IMolecule;
import etomica.normalmode.LatticeSumMolecularCrystal.AtomicTensorAtomicPair;
import etomica.normalmode.NormalModesMolecular;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Space3D;
import etomica.space3d.Tensor3D;
import etomica.species.SpeciesGeneral;
import etomica.units.Calorie;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.units.Mole;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;

public class ClathrateHarmonicFE extends Simulation {
    protected static double[] initialU;
    protected Box box;
    protected PotentialMaster potentialMaster;
    protected SpeciesGeneral species;
    protected MeterPotentialEnergy meterPE;
    protected Potential2SoftSpherical potentialLJ;
    protected Potential2SoftSphericalLS potentialLJLS;
    protected EwaldSummation potentialES;

    public ClathrateHarmonicFE(Space space, int[] nC, double rCutRealES, double rCutLJ, double[] a0_sc, int numMolecule, String configFileName, boolean isIce, double kCut, boolean includeM) {
        super(space);
        species = SpeciesWater4P.create();
        addSpecies(species);
        box = this.makeBox(new BoundaryRectangularPeriodic(space, a0_sc));
        box.setNMolecules(species, numMolecule);
        DoubleStorage charges = box.getAtomStorage(Tokens.doubles(new ChargeAgentSourceRPM(species, isIce)));
        double sigma, epsilon;
        if (isIce) {
            sigma = 3.1668;
            epsilon = Kelvin.UNIT.toSim(106.1);//TIP4P/Ice
//			sigma = 3.1589; epsilon = Kelvin.UNIT.toSim(93.2);//TIP4P/2005			
        } else {//TIP4P
            double A = 600E3; // kcal A^12 / mol
            double C = 610.0; // kcal A^6 / mol
            double s6 = A / C;
            sigma = Math.pow(s6, 1.0 / 6.0);
            epsilon = Mole.UNIT.fromSim(Calorie.UNIT.toSim(C / s6 * 1000)) / 4.0;
//            sigma = 3.154; epsilon = Kelvin.UNIT.toSim(78.0); //TIP4P           
        }
        //To get Dij(LJ)
        potentialLJ = new P2LennardJones(space, sigma, epsilon);
        //To get U
        potentialLJLS = new Potential2SoftSphericalLS(space, rCutLJ, a0_sc, new P2LennardJones(space, sigma, epsilon));
        potentialES = new EwaldSummation(box, charges, space, kCut, rCutRealES);
//XXXX Potential Master
        potentialMaster = new PotentialMaster();
        potentialMaster.addPotential(potentialLJLS, new AtomType[]{species.getTypeByName("O"), species.getTypeByName("O")});
        potentialMaster.addPotential(potentialES, new AtomType[0]);
        potentialLJLS.setBox(box);

        if (includeM) {
            ConfigurationFile config = new ConfigurationFile(configFileName);////to duplicate with M point!
            ConfigurationFileBinary.replicate(config, box, nC, space);
        } else {
            ConfigurationFileTIP4P config = new ConfigurationFileTIP4P(configFileName, space, isIce);//to duplicate w.o. M point!
            ConfigurationFileBinary.replicate(config, box, nC, space);
        }

        MeterPotentialEnergy meterPotentialEnergy = new MeterPotentialEnergy(potentialMaster, box);
        AtomPair selfAtomLJ = new AtomPair(box.getLeafList().get(2), box.getLeafList().get(2));
        double selfELJ = potentialLJLS.energy(selfAtomLJ);// = 0 if nLJshells=0
        double E = Joule.UNIT.fromSim(meterPotentialEnergy.getDataAsScalar() / numMolecule + selfELJ) * 1.0E-3 * Constants.AVOGADRO;
        System.out.println(" E (kJ/mol) = " + E);

//		System.out.println("******************************************************************");
//		double s = Math.sqrt(Math.PI*rCutRealES*nKs/nC[0]/a0);
//		alpha = s/rCutRealES;
//		System.out.println(" s = Sqrt(PI*rCutRealES*nKs/L) = " + s);
//		System.out.println(" alpha = s/rCutRealES = " + alpha);
//		System.out.println("******************************************************************");
//        int nAtoms = box.getLeafList().getAtomCount();
//        double Q2=0;
//        for (int i=0; i < nAtoms; i++){//H
//            IAtom atom = box.getLeafList().getAtom(i);
//            double charge = atomAgentManager.getAgent(atom).charge;
//            Q2 += charge*charge;
//        }
//        Q2 /= numMolecule;
//        double boxSize = box.getBoundary().getBoxSize().getX(0);
//
//        System.out.println("Q/N*(s/2/alpha/L^3)^1/2 * exp(-s^2)/s/s  =  " +
//          Joule.UNIT.fromSim(Q2*Math.sqrt(s/2/alpha/boxSize/boxSize/boxSize)*Math.exp(-s*s) /s/s)  *  1.0E-3*Constants.AVOGADRO);
//        System.exit(0);


//		double gmToSim = Gram.UNIT.toSim(1);
//		double cmToSim = 1.0e8;
//		double mH2O = species.getOxygenType().getMass() + 2.0 * species.getHydrogenType().getMass();  
//		double I_H2O = (1.0220e-40*gmToSim*cmToSim*cmToSim)*(2.9376e-40*gmToSim*cmToSim*cmToSim)*(1.9187e-40*gmToSim*cmToSim*cmToSim); 
//		double I_H2O = 1.22; 
//		double[] lnSumN = new double[] {58.461308350846494,58.816210045241121,58.844810323435830,58.850744777620520
//			,58.852595489204397,58.853322163833958,58.853653353151124,58.853820899859144,58.853912442563434
//			,58.853965490875808};
//		double[] lnSumN = new double[] {56.935724982969120,57.282460042436469,57.310064109044561,57.301594015133816,57.246984160655664,57.184578618197946,57.094033481750614,57.019428683165764};
//		double N=46.0;
//		for (int t=180;t<=280;t+=10){
//		double t=210.0;
//	   for (int i=0;i<lnSumN.length;i++){
//		   N = (i+1)*(i+1)*(i+1)*46;
//			double T = t;
//			double kT = Kelvin.UNIT.toSim(t);
//			double kTkJMol = 0.00831435 * T;
//			double lnSum = 2730.726971235;
//			double lnSum = 2619.04;
//			double lnSum =  2724.8282433741542;
//			System.out.println(Constants.PLANCK_H);
//			double Atra = 3.0*(N-1)/(2.0*N) * Math.log( 2.0*Math.PI*mH2O*kT/(Constants.PLANCK_H*Constants.PLANCK_H) );
//			Atra *= -kTkJMol;
//			
//			double Arot = 3.0/2.0 * Math.log(2.0*Math.PI*kT)
//					+ 1.0/2.0 * Math.log(I_H2O)
//					- 3.0*Math.log(Constants.PLANCK_H);
//			Arot *= - kTkJMol;
//			
//			double Aconf = lnSum/N/2.0 - (6.0*N-3.0)/(2.0*N) * Math.log(2.0*Math.PI*kT);
//			double Aconf = lnSumN[i]/2.0 - (6.0*N-3.0)/(2.0*N) * Math.log(2.0*Math.PI*kT);
//			Aconf *= kTkJMol;
//			double U0 = -56.22;
//			double U0 = -54.21;
//			double Atot = Atra + Arot + Aconf + U0 - kTkJMol*Math.log(3.0/2.0) + kTkJMol * Math.log(2);// -3.0/2.0 * Math.log(N)/N; ;
//			System.out.println(T +" "+ Atot);
//			System.out.println(N +" "+ 1.0/N +" "+ (Atot));
//			}
//        System.exit(0);
    }

    public static void main(String[] args) {
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        final int[] nC = params.nC;
        int nBasis = params.nBasis;
        double kCut = params.kCut;
        double temperature = params.temperature;
        String configFile = params.configFile;
        final double[] a0 = params.a0;
        boolean waveVectorMethod = params.waveVectorMethod;
        boolean isIce = params.isIce;
        boolean includeM = params.includeM;
        double rCutRealES = params.rCutRealES;
        double rCutLJ = params.rCutLJ;

        final double[] a0_sc = new double[]{a0[0] * nC[0], a0[1] * nC[1], a0[2] * nC[2]};

        final int[] nLJShells = new int[]{(int) Math.ceil(rCutLJ / a0_sc[0] - 0.49999), (int) Math.ceil(rCutLJ / a0_sc[1] - 0.49999), (int) Math.ceil(rCutLJ / a0_sc[2] - 0.49999)};

        final double rCutLJ2 = rCutLJ * rCutLJ;

        System.out.println(" nX = " + nC[0]);
        System.out.println(" rCutLJ = " + rCutLJ);
        System.out.println(" kCut = " + kCut);


        int numMolecule = nBasis * nC[0] * nC[1] * nC[2];
        final Space space = Space3D.getInstance();
        final ClathrateHarmonicFE sim = new ClathrateHarmonicFE(space, nC, rCutRealES, rCutLJ, a0_sc, numMolecule, configFile, isIce, kCut, includeM);


        //Anonymous class enables you to make your code more concise.
        //They enable you to declare and instantiate a class at the same time! :)
        AtomicTensorAtomicPair atomicTensorAtomicPair = new AtomicTensorAtomicPair() {
            final Tensor identity = new Tensor3D(new double[][]{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}});
            Vector Lxyz = space.makeVector();
            Vector dr = space.makeVector();
            Vector drTmp = space.makeVector();
            Tensor D3ESLJ = space.makeTensor();
            Tensor tmpD3LJ = space.makeTensor();

            //dr == r0 - r1
            public Tensor atomicTensor(IAtom atom0, IAtom atom1) {
                D3ESLJ.E(0);
                //LJ
                if (atom0.getType() == sim.species.getTypeByName("O") && atom1.getType() == sim.species.getTypeByName("O")) {
                    dr.Ev1Mv2(atom0.getPosition(), atom1.getPosition());
                    sim.box.getBoundary().nearestImage(dr);

                    //LS of LJ
                    for (int nx = -nLJShells[0]; nx <= nLJShells[0]; nx++) {
                        Lxyz.setX(0, nx * a0_sc[0]);
                        for (int ny = -nLJShells[1]; ny <= nLJShells[1]; ny++) {
                            Lxyz.setX(1, ny * a0_sc[1]);
                            for (int nz = -nLJShells[2]; nz <= nLJShells[2]; nz++) {
                                Lxyz.setX(2, nz * a0_sc[2]);
                                drTmp.Ev1Pv2(dr, Lxyz);
                                double dr2Tmp = drTmp.squared();
                                if (dr2Tmp > rCutLJ2 || dr2Tmp == 0) continue;
                                tmpD3LJ.Ev1v2(drTmp, drTmp);
                                double dW = sim.potentialLJ.du(dr2Tmp);
                                double d2W = sim.potentialLJ.d2u(dr2Tmp);
                                tmpD3LJ.TE(1.0 / (dr2Tmp * dr2Tmp) * (dW - d2W));
                                tmpD3LJ.PEa1Tt1(-dW / dr2Tmp, identity);
                                D3ESLJ.PE(tmpD3LJ);
                            }
                        }
                    }
                    //ES
                } else if (atom0.getType() != sim.species.getTypeByName("O") && atom1.getType() != sim.species.getTypeByName("O")) {
                    D3ESLJ.PE(sim.potentialES.secondDerivative(atom0, atom1));
                }
                return D3ESLJ;
            }
        }; // End AtomicTensorAtomicPair


        if (false) {
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "string");
            final DisplayBox display = new DisplayBox(sim, sim.box);
            simGraphic.add(display);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box).getColorScheme()).setColor(sim.species.getTypeByName("H"), Color.GREEN);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box).getColorScheme()).setColor(sim.species.getTypeByName("O"), Color.RED);
            //Sabry
            ((DiameterHashByType) ((DisplayBox) simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species.getTypeByName("O"), 2.0);
            ((DiameterHashByType) ((DisplayBox) simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species.getTypeByName("H"), 1.0);

            ((DisplayBoxCanvasG3DSys) simGraphic.getDisplayBox(sim.box).canvas).setBackgroundColor(Color.white);
            ((DisplayBoxCanvasG3DSys) simGraphic.getDisplayBox(sim.box).canvas).setBoundaryFrameColor(Color.blue);
            simGraphic.makeAndDisplayFrame();
            if (!true) {
                return;
            }
        }

//Sabry		

        if (false) {
            VectorStorage forces = sim.box.getAtomStorage(Tokens.FORCES);
            PotentialCalculationForceSum pcForce = new PotentialCalculationForceSum(sim.box.getAtomStorage(Tokens.FORCES));

            IteratorDirective id = new IteratorDirective();
            id.includeLrc = false;
            sim.potentialMaster.calculate(sim.box, id, pcForce);

            Vector fn = sim.space.makeVector();
            Vector f4 = sim.space.makeVector();
            Vector dr = sim.space.makeVector();
            Vector torque = sim.space.makeVector();
//
//	        for (int m = 0; m < 46; m++){
//		        IMolecule iMol = sim.box.getMoleculeList().getMolecule(m);
//		        IAtomList atoms = iMol.getChildList();
//		        for (int n = 0; n < atoms.getAtomCount(); n++){
//		        	IAtom atomj = atoms.getAtom(n);
//		            fn.E(((IntegratorVelocityVerlet.MyAgent)atomAgentManager.getAgent(atomj)).force);
//		        	f4.PE(fn);
//	        		dr.Ev1Mv2(atomj.getPosition(), atoms.getAtom(2).getPosition());
//	        		sim.box.getBoundary().nearestImage(dr);
//	        		dr.XE(fn);
//	        		torque.PE(dr);
//		        }
//		    	//System.out.println(f4);
//		    	System.out.println(torque);
//		    	f4.E(0);
//		    	torque.E(0);
//	        }
//	        System.out.println("++++++");

            RotationTensor3D rTensor = new RotationTensor3D();
            Vector f4dx = sim.space.makeVector();
            Vector dx = space.makeVector();
            dx.setX(0, 0.0001);
            dx.setX(1, 0.0);
            dx.setX(2, 0.0);
            for (int m = 0; m < 46; m++) {
                IMolecule iMol = sim.box.getMoleculeList().get(m);
                IAtomList atoms = iMol.getChildList();
//		        for (int n = 0; n < atoms.getAtomCount(); n++){
//		        	IAtom atomj = atoms.getAtom(n);
//		            pcForce.reset();
//		            sim.potentialMaster.calculate(sim.box, id, pcForce);
//		            fn.E(((IntegratorVelocityVerlet.MyAgent)atomAgentManager.getAgent(atomj)).force);
//		        	System.out.println(fn);
//		        	f4.PE(fn);
//		        }
//		    	System.out.println(f4);
//		        System.out.println("================");
                //  dx of molecule m'
                int Mp = 0, M = 1;
                if (m == Mp) { //2ns molecule (m=1)
                    IMolecule mpMol = sim.box.getMoleculeList().get(m);
//			        for(int j=0;j<4;j++){
//				        atoms.getAtom(j).getPosition().PE(dx);
//			        	System.out.println(atoms.getAtom(j).getPosition());
//			        }
                    Vector axes = sim.space.makeVector();
                    axes.E(new double[]{1.0, 0.0, 0.0});
//			        axes.E(new double [] {0.0, 1.0, 0.0});
//			        axes.E(new double [] {0.0, 0.0, 1.0});
                    double dTheta = 0.001;
                    rTensor.setRotationAxis(axes, dTheta);
                    System.out.println("dTheta = " + dTheta);
                    MinimizationTIP4P.doTransform(mpMol, sim.box.getBoundary(), rTensor);
                    System.out.println("");
//	    	        for(int j=0;j<4;j++){
//			        	System.out.println(atoms.getAtom(j).getPosition());
//			        }
                }
                //measure force/torque on molecule m
                if (m == M) {//5th molecule (m=4)

                    pcForce.reset();
                    sim.potentialMaster.calculate(sim.box, id, pcForce);
                    for (int j = 0; j < 4; j++) {
                        IAtom atomj = atoms.get(j);

                        Vector fndx = space.makeVector();
                        fndx.E(forces.get(atomj));
//			            System.out.println(fndx);
                        f4dx.PE(fndx);
                        dr.Ev1Mv2(atomj.getPosition(), atoms.get(2).getPosition());
                        sim.box.getBoundary().nearestImage(dr);
                        dr.XE(fndx);
                        torque.PE(dr);
                    }
                    System.out.println("f4dx = " + f4dx);
                    System.out.println("torque = " + torque);
                    System.out.println("");
                    f4dx.E(0);
                    torque.E(0);
                    System.exit(0);
                }
            }//end for m
        }
        //end Sabry


        NormalModesMolecular normalModes;
        Primitive primitive;

        if (waveVectorMethod == true) {
            primitive = new PrimitiveOrthorhombic(space, a0[0], a0[1], a0[2]);
            normalModes = new NormalModesMolecular(sim.species, waveVectorMethod, sim.potentialMaster, sim.box, nC, primitive, nBasis, atomicTensorAtomicPair, space);
        } else {//Direct BIG D
            primitive = new PrimitiveOrthorhombic(space, a0[0], a0[1], a0[2]);
            normalModes = new NormalModesMolecular(sim.species, waveVectorMethod, sim.potentialMaster, sim.box, nC, primitive, nC[0] * nC[1] * nC[2] * nBasis, atomicTensorAtomicPair, space);
        }
        normalModes.calculateModes();
    }


    public static class SimParams extends ParameterBase {
        //    	public String configFile = "config_from_paper_HHO_shiftedL_2_sI";
        public String configFile = "finalPos";
        public double rCutLJ = 400;
        public double rCutRealES = 14.0;
        public double[] a0 = new double[]{12.03, 12.03, 12.03};//sI
        public int nBasis = 46;//sI
        //		public int nBasis = 136;//sII
//		public int nBasis = 68;//sH
        public double kCut = 2.6;
        //		public double[] a0 = new double[]{17.31, 17.31, 17.31};//sII
//		public double[] a0 = new double[]{12.21,21.15, 10.14};//sH
        public double temperature = 0.1;
        public boolean waveVectorMethod = true;
        public boolean isIce = false;
        public boolean includeM = true;
        int nX = 1;
        public int[] nC = new int[]{nX, nX, nX};
    }
}
