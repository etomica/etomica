package etomica.models.nitrogen;

import java.util.Arrays;

import etomica.action.IAction;
import etomica.action.MoleculeActionTranslateTo;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IMolecule;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtomPositionDefinition;
import etomica.atom.MoleculeAgentManager;
import etomica.atom.MoleculeAgentManager.MoleculeAgentSource;
import etomica.atom.MoleculePair;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataTensor;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorRigidIterative;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.listener.IntegratorListenerAction;
import etomica.nbr.list.molecule.BoxAgentSourceCellManagerListMolecular;
import etomica.nbr.list.molecule.NeighborListManagerSlantyMolecular;
import etomica.nbr.list.molecule.PotentialMasterListMolecular;
import etomica.normalmode.BasisBigCell;
import etomica.normalmode.MCMoveMoleculeCoupled;
import etomica.normalmode.MeterHarmonicEnergy;
import etomica.potential.PotentialCalculationTorqueSum;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.ISpace;
import etomica.space.Tensor;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Space3D;
import etomica.units.Degree;
import etomica.units.Kelvin;
import etomica.units.Pascal;
import etomica.units.Pixel;



/**
 * Simulation class for nitrogen molecules
 * beta-N2 crystal Structure
 * 
 * @author Tai Boon Tan
 *
 */
public class MinimizationBetaNitrogenModel extends Simulation{

	public MinimizationBetaNitrogenModel(ISpace space, int[] nC, double temperature, double pressure, double newScale) {
		super(space);
		this.space = space;
//		double a = 3.840259;//3.854;
//		double c = 6.263463;//6.284; 
		
		BoxAgentSourceCellManagerListMolecular boxAgentSource = new BoxAgentSourceCellManagerListMolecular(this, null, space);
	    BoxAgentManager boxAgentManager = new BoxAgentManager(boxAgentSource);
	     
		
		double ratio = 1.631;
		double density = 0.0236;
		double a = Math.pow(4.0/(Math.sqrt(3.0)*ratio*density), 1.0/3.0);
		double c = a*ratio;
		int numMolecule = nC[0]*nC[1]*nC[2]*2;
		int nCa = (int)Math.pow(numMolecule/1.999999999, 1.0/3.0);
		
		Basis basisHCP = new BasisHcp();
		Basis basis = new BasisBigCell(space, basisHCP, new int[]{nC[0], nC[1], nC[2]});
		
		species = new SpeciesN2(space);
		addSpecies(species);
		
		box = new Box(space);
		addBox(box);
		box.setNMolecules(species, numMolecule);		
		int [] nCells = new int[]{1,1,1};
		
		IVector[] boxDim = new IVector[3];
		boxDim[0] = space.makeVector(new double[]{nC[0]*a, 0, 0});
		boxDim[1] = space.makeVector(new double[]{-nC[1]*a*Math.cos(Degree.UNIT.toSim(60)), nC[1]*a*Math.sin(Degree.UNIT.toSim(60)), 0});
		boxDim[2] = space.makeVector(new double[]{0, 0, nC[2]*c});
		
		Boundary boundary = new BoundaryDeformablePeriodic(space,boxDim);
		primitive = new PrimitiveHexagonal(space, (nC[0])*a, nC[2]*c);
		
		coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDef.setIsBeta();
		coordinateDef.setOrientationVectorBeta(space);
		coordinateDef.initializeCoordinates(nCells);

		
		double[] u = new double[20];

		if(true){
			BetaPhaseLatticeParameter parameters = new BetaPhaseLatticeParameter();
			double[][] param = parameters.getParameter(density);
			
			int kParam=0;
			for (int i=0; i<param.length;i++){
				for (int j=0; j<param[0].length;j++){
					u[kParam]=param[i][j];
					kParam++;
				}	
			}
			
//			System.out.println("*************Parameters*************");
//			for (int i=0; i<u.length;i++){
//				System.out.print(u[i] +", ");
//				if((i+1)%5==0){
//					System.out.println();
//				}
//			}
			
			int numDOF = coordinateDef.getCoordinateDim();
			double[] newU = new double[numDOF];
			if(true){
				for(int j=0; j<numDOF; j+=10){
					if(j>0 && j%(nCa*10)==0){
						j+=nCa*10;
						if(j>=numDOF){
							break;
						}
					}
					for(int k=0; k<10;k++){
						newU[j+k]= u[k];
					}
				}
				
				for(int j=nCa*10; j<numDOF; j+=10){
					if(j>nCa*10 && j%(nCa*10)==0){
						j+=nCa*10;
						if(j>=numDOF){
							break;
						}
					}
					for(int k=0; k<10;k++){
						newU[j+k]= u[k+10];
					}
				}
			}

			coordinateDef.setToU(box.getMoleculeList(), newU);
			coordinateDef.initNominalU(box.getMoleculeList());
			
		}

		
		box.setBoundary(boundary);
		double rC = a*nC[0]*0.475;
		//System.out.println("Truncation Radius: " + rC);
		potential = new P2Nitrogen(space, rC);
		potential.setBox(box);
		
		
		
		potentialMaster = new PotentialMasterListMolecular(this, rC, boxAgentSource, boxAgentManager, new NeighborListManagerSlantyMolecular.NeighborListSlantyAgentSourceMolecular(rC, space), space);
	    potentialMaster.addPotential(potential, new ISpecies[]{species, species});
	    int cellRange = 6;
	    potentialMaster.setRange(rC);
	    potentialMaster.setCellRange(cellRange); 
	    potentialMaster.getNeighborManager(box).reset();
	      
	    potential.setRange(Double.POSITIVE_INFINITY);
	    int potentialCells = potentialMaster.getNbrCellManager(box).getLattice().getSize()[0];
	    if (potentialCells < cellRange*2+1) {
	    	throw new RuntimeException("oops ("+potentialCells+" < "+(cellRange*2+1)+")");
	    }
		
		MeterPotentialEnergy meterPotentialEnergy = new MeterPotentialEnergy(potentialMaster);
		meterPotentialEnergy.setBox(box);
		double latticeEnergy = 0; // meterPotentialEnergy.getDataAsScalar();
		System.out.println("Lattice Energy: "+ latticeEnergy);


//		NormalModesFromFile nm = new NormalModesFromFile("beta"+numMolecule+"_2ndDer_d"+density, 3);
//		meterHarm = new MeterHarmonicEnergy(coordinateDef, nm);
		
//		meterPE = new MeterPotentialEnergy(potentialMaster);
//		meterPE.setBox(box);
//		double[] uCoord = new double[coordinateDef.getCoordinateDim()];
//		
//		for (double i=-0.1; i<=0.101; i+=0.002){
//			u = new double[]{i, 0, 0, 0, 0,
//					         0, 0, 0, 0, 0, 
//					         0, 0, 0, 0, 0, 
//					         0, 0, 0, 0, 0};
//			//uCoord[89] = i;
//			int numDOF = coordinateDef.getCoordinateDim();
//			if(true){
//				for(int j=0; j<numDOF; j+=10){
//					if(j>0 && j%(nCa*10)==0){
//						j+=nCa*10;
//						if(j>=numDOF){
//							break;
//						}
//					}
//					for(int k=0; k<10;k++){
//						uCoord[j+k]= u[k];
//					}
//				}
//				
//				for(int j=nCa*10; j<numDOF; j+=10){
//					if(j>nCa*10 && j%(nCa*10)==0){
//						j+=nCa*10;
//						if(j>=numDOF){
//							break;
//						}
//					}
//					for(int k=0; k<10;k++){
//						uCoord[j+k]= u[k+10];
//					}
//				}
//			}
			
			
//			coordinateDef.setToU(box.getMoleculeList(), uCoord);
//			double pe = meterPE.getDataAsScalar();
//		//	double he = meterHarm.getDataAsScalar();
//			System.out.println(i+" "+ (pe-latticeEnergy)/numMolecule);
//		}
//		
//		System.exit(1);
		
		
		MCMoveMoleculeCoupled move = new MCMoveMoleculeCoupled(potentialMaster,getRandom(),space);
		move.setBox(box);
		move.setPotential(potential);
		move.setDoExcludeNonNeighbors(true);
		//move.setStepSize(Kelvin.UNIT.toSim(temperature));
		//((MCMoveStepTracker)move.getTracker()).setNoisyAdjustment(true);
		   
		MCMoveRotateMolecule3D rotate = new MCMoveRotateMolecule3D(potentialMaster, getRandom(), space);
		rotate.setBox(box);
		
		MCMoveVolume mcMoveVolume = new MCMoveVolume(this, potentialMaster, space);
		mcMoveVolume.setBox(box);
		pressure *= 1e9;
		mcMoveVolume.setPressure(Pascal.UNIT.toSim(pressure));
		
		integrator = new IntegratorMC(this, potentialMaster);
		integrator.getMoveManager().addMCMove(move);
		integrator.getMoveManager().addMCMove(rotate);
		//integrator.getMoveManager().addMCMove(mcMoveVolume);
		integrator.setBox(box);
		
		integrator.setTemperature(Kelvin.UNIT.toSim(temperature));
		
		activityIntegrate = new ActivityIntegrate(integrator);
		getController().addAction(activityIntegrate);
	}
	
	public static void main (String[] args){
//		System.out.println("pressure " + Pascal.UNIT.fromSim(31.48928359791796)/1e9);
//		System.exit(1);
		int nC0 = 8; 
		int nC1 = 8; 
		int nC2 = 8;
		double temperature =0.002; // in Unit Kelvin
		double pressure = 0.0; //in Unit GPa
		long simSteps = 100000;
		double newScale = 1.0;
		if(args.length > 1){
			simSteps = Long.parseLong(args[1]);
		}
		if(args.length > 2){
			temperature = Double.parseDouble(args[2]);
		}
		if(args.length > 3){
			pressure = Double.parseDouble(args[3]);
		}
		if(args.length > 4){
			nC0 = Integer.parseInt(args[4]);
		}
		if(args.length > 5){
			nC1 = Integer.parseInt(args[5]);
		}
		if(args.length > 6){
			nC2 = Integer.parseInt(args[6]);
		}
		if(args.length > 7){
			newScale = Double.parseDouble(args[7]);
		}
		
		int[] nC = new int []{nC0,nC1,nC2};
		int numMolecule =nC[0]*nC[1]*nC[2]*2;
		
		String filename = "betaN2_nA"+numMolecule+"_T"+temperature;
		
		if(args.length > 0){
			filename = args[0];
		} 
//		System.out.println("Running beta-N2 crystal structure simulation with " + simSteps + " steps" );
//		System.out.println("num Molecules: " + numMolecule+ " ; temperature: " + temperature
//				+"K ; pressure: "+ pressure+"GPa");
//		System.out.println("With volume scaling of " + newScale);
//		System.out.println("Output file: " + filename + "\n");

		double density = 0.0230;
		final SimulationAlphaNitrogenModel sim = new SimulationAlphaNitrogenModel(Space3D.getInstance(3), nC, temperature, density);
	    
		final MeterPotentialEnergy meterPotentialEnergy = new MeterPotentialEnergy(sim.potentialMaster);
		meterPotentialEnergy.setBox(sim.box);
		final double latticeEnergy = 0; //meterPotentialEnergy.getDataAsScalar();
//		System.out.println("Lattice Energy (per molecule): "+ latticeEnergy/numMolecule);

		PotentialCalculationTorqueSum pcForce = new PotentialCalculationTorqueSum();
        MoleculeAgentSource molAgentSource = new MoleculeAgentSource() {
            
            public void releaseAgent(Object agent, IMolecule molecule) {
            }
            
            public Object makeAgent(IMolecule mol) {
                return new IntegratorRigidIterative.MoleculeAgent(sim.space);
            }
            
            public Class getMoleculeAgentClass() {
                return IntegratorRigidIterative.MoleculeAgent.class;
            }
        };
        MoleculeAgentManager molAgentManager = new MoleculeAgentManager(sim, sim.box, molAgentSource);
        pcForce.setMoleculeAgentManager(molAgentManager);
        double[] d = new double[16];
        MoleculeActionTranslateTo translator = new MoleculeActionTranslateTo(sim.space);
        AtomPositionGeometricCenter pos = new AtomPositionGeometricCenter(sim.space);
        translator.setAtomPositionDefinition(pos);
        IVectorMutable p = sim.space.makeVector();
        int nA = 16;
        double step1 = 0;
        double[] x0 = new double[12];
        IVectorMutable[] torques = new IVectorMutable[4];
        IVectorMutable[] axes = new IVectorMutable[4];
        for (int i=0; i<4; i++) {
            axes[i] = sim.space.makeVector();
            torques[i] = sim.space.makeVector();
        }
        for (int outer = 0; outer < 30; outer++) {
            System.out.println("**** "+outer+" ****");
        double totalD = 0;
        double step = 0;
        double radianFac = 1e5;
//        System.out.println("box size "+sim.box.getBoundary().getBoxSize());
        IMolecule mol0 = sim.box.getMoleculeList().getMolecule(0);
//        new ConformationNitrogen(sim.space).initializePositions(mol0.getChildList());
//        IVectorMutable com = sim.space.makeVector(new double[]{-10,0,0});
//        translator.setDestination(com);
//        translator.actionPerformed(mol0);
        IMolecule mol1 = sim.box.getMoleculeList().getMolecule(1);
//        new ConformationNitrogen(sim.space).initializePositions(mol1.getChildList());
//        com = sim.space.makeVector(new double[]{10,0,0});
//        translator.setDestination(com);
//        translator.actionPerformed(mol1);
        RotationTensor3D rTensor = new RotationTensor3D();
//        rTensor.setAxial(1, 2);
//        doTransform(mol1, pos, rTensor);
//        double u01 = sim.potential.energy(new MoleculePair(mol0, mol1));
//        System.out.println("u01 "+u01);
        IVector torque = null; //sim.potential.gradientAndTorque(new MoleculePair(mol0, mol1))[1][0];
//        System.out.println("torque01 "+torque);
        IteratorDirective id0 = new IteratorDirective(null, mol0);
//        sim.potentialMaster.calculate(sim.box, id0, pcForce);
//        System.out.println("torque0 "+torque);
//        pcForce.reset();
        meterPotentialEnergy.setTarget(mol0);
//        double u0 = meterPotentialEnergy.getDataAsScalar();
//        System.out.println("u0 "+u0);
        rTensor.setAxial(1, 0.01);
        IVectorMutable com = sim.space.makeVector();
        com.E(pos.position(mol0));
        
//        torque = sim.potential.gradientAndTorque(new MoleculePair(mol0, mol1))[1][0];
//        System.out.println("torque01 "+torque);
//        for  (int i=0; i<2; i++) {
//            com.setX(0, com.getX(0) + 0.01);
//            translator.setDestination(com);
//            translator.actionPerformed(mol0);
//            System.out.println("u01 "+sim.potential.energy(new MoleculePair(mol0, mol1)));
//            torque = sim.potential.gradientAndTorque(new MoleculePair(mol0, mol1))[1][0];
//            System.out.println("torque01 "+torque);
//            if (i==0) {
//                DataTensor d2 = sim.potential.secondDerivative(new MoleculePair(mol0, mol1));
//                System.out.println(d2);
//            }
//        }
//        com.setX(0, com.getX(0) - 0.02);
//        translator.setDestination(com);
//        translator.actionPerformed(mol0);
//        System.out.println("\n\n\n");

        
        IVectorMutable uAxis = sim.coordinateDef.getMoleculeOrientation(mol0)[2];
        System.out.println("uAxis: " + uAxis);
//        IVectorMutable uAxis = sim.space.makeVector(new double[]{0.2, 0.4, 0.6});
      
        

//        uAxis.normalize();
//        sim.potential.uAxis = uAxis;
        
        double dx = 0.001;
        double dtheta = dx;
        double d0 = 0;
//        for (int i=0; i<3; i++) {
//            rTensor.setRotationAxis(uAxis, dtheta);
//            double u0 = 0;
//            for (int j=0; j<3; j++) {
//                double u = sim.potential.energy(new MoleculePair(mol0, mol1));
//                if  (j==0) {
//                    u0 = u;
//                }
//                else if (j==2) {
//                    double du = (u-u0)/(2*dx);
//                    System.out.println("d "+du);
//                    if (i==0) {
//                        d0 = du;
//                    }
//                    else if (i==2) {
//                        System.out.println("d2 "+(du-d0)/(2*dtheta));
//                    }
//                }
//                System.out.println(i+" "+j+" u01 "+u);
//                torque = sim.potential.gradientAndTorque(new MoleculePair(mol0, mol1))[1][0];
//                System.out.println(i+" "+j+" torque01 "+torque+" "+torque.dot(uAxis));
//                if (i==1 && j==1) {
//                    DataTensor d2 = sim.potential.secondDerivative(new MoleculePair(mol0, mol1));
//                    d2.x.TE(2*dx);
//                    System.out.println(d2);
//                }
//                doTransform(mol0, pos, rTensor);
//            }
//            rTensor.setRotationAxis(uAxis, -3*dtheta);
//            doTransform(mol0, pos, rTensor);
//
//            com.setX(0, com.getX(0) + dx);
//            translator.setDestination(com);
//            translator.actionPerformed(mol0);
//        }
        
        
//        System.exit(1);
//	        for  (int i=0; i<2; i++) {
//	            doTransform(mol0, pos, rTensor);
//	            System.out.println("u01 "+sim.potential.energy(new MoleculePair(mol0, mol1)));
//	            torque = sim.potential.gradientAndTorque(new MoleculePair(mol0, mol1))[1][0];
//	            System.out.println("torque01 "+torque);
//	            System.out.println("u0 "+meterPotentialEnergy.getDataAsScalar());
//	            sim.potentialMaster.calculate(sim.box, id0, pcForce);
//	            torque = ((IntegratorRigidIterative.MoleculeAgent)molAgentManager.getAgent(mol0)).torque;
//	            System.out.println("torque0 "+torque);
//	            pcForce.reset();
//	        }
//        System.exit(1);
		    for (int iter=0; iter<3; iter++) {
		        double[] g = new double[12];
		        double t = 0;
		        for (int i=0; i<4; i++) {
		            IMolecule iMol = sim.box.getMoleculeList().getMolecule(i<2 ? i : i+(nA-2));
		            IteratorDirective id = new IteratorDirective(null, iMol);
		            sim.potentialMaster.calculate(sim.box, id, pcForce);
		            IVector f = ((IntegratorRigidIterative.MoleculeAgent)molAgentManager.getAgent(iMol)).force;
		            for (int j=0; j<3; j++) {
		                g[i*3+j] = -f.getX(j);
		                t += g[i*3+j]*g[i*3+j];
		            }
		            torques[i].E(((IntegratorRigidIterative.MoleculeAgent)molAgentManager.getAgent(iMol)).torque);
		            System.out.println(torques[i]);
		            if (iter == 0) {
		                double t2 = torques[i].squared();
		                double sqrtT = Math.sqrt(t2);
		                t += t2/(radianFac*radianFac);
		                d[12+i] = sqrtT/radianFac;
		                axes[i].E(torques[i]);
		                axes[i].TE(1/sqrtT);
		            }
		            pcForce.reset();
		        }
		        System.out.println(Arrays.toString(g));
		        if (iter == 0) {
		            t = Math.sqrt(t);
		//                double[] t2 = new double[3];
		            for (int j=0; j<12; j++) {
		                d[j] = -g[j]/t;
		//                    t2[j%3] += d[j];
		                totalD += g[j]*d[j];
		            }
		            for (int j=12; j<16; j++) {
		                d[j] /= t;
		                totalD -= d[j]*torques[j-12].dot(axes[j-12])/radianFac;
		            }
		            System.out.println(Arrays.toString(d));
		//                for (int j=0; j<3; j++) {
		//                    t2[j] /= 4;
		//                }
		//                t = 0;
		//                for (int j=0; j<12; j++) {
		//                    d[j] -= t2[j%3];
		//                    t += d[j]*d[j];
		//                }
		//                t = Math.sqrt(t);
		//                for (int j=0; j<12; j++) {
		//                    d[j] /= t;
		//                }
		//                System.out.println("[COM] "+Arrays.toString(d));
		            System.out.println("totalD "+totalD);
		            if (step1 == 0) {
		                step = 0.0000001;
		            }
		            else {
		                step = step1/5;
		            }
		            System.out.println("step "+step);
		        }
		        else {
		            double newTotalD = 0;
		            for (int j=0; j<12; j++) {
		                newTotalD += g[j]*d[j];
		            }
		            System.out.println("totalD "+newTotalD);
		//                System.out.println(totalD + " "+step+" "+newTotalD+" "+totalD);
		            double oldStep = step; 
		            step = - newTotalD * step / (newTotalD - totalD);
		            if (iter == 1) {
		                step1 = oldStep + step;
		            }
		            totalD = newTotalD;
		            System.out.println("step "+step);
		        }
		
		        for (int i=0; i<4; i++) {
		            for (int k=i%2; k<numMolecule; k+=2) {
		                boolean isA = (k/nA)%2 == 0;
		                if ((i<2 && !isA) || (i>1 && isA)) continue;
		//                    System.out.println(i+" "+k);
		                IMolecule iMol = sim.box.getMoleculeList().getMolecule(k);
		                p.E(pos.position(iMol));
		                for (int j=0; j<3; j++) {
		                    if (x0[i*3+j] == 0) {
		                        x0[i*3+j] = p.getX(j);
		                    }
		                    p.setX(j, p.getX(j)+step*d[i*3+j]);
		                }
		                translator.setDestination(p);
		                translator.actionPerformed(iMol);
		                
		                rTensor.setRotationAxis(axes[i], d[12+i]/radianFac);
		                doTransform(iMol, pos, rTensor);
		            }
		        }
		    }
        }
        
        
        
        //DONE with minimization
        
        
        double[] xf = new double[12];
        double disp = 0;
        for (int i=0; i<4; i++) {
            IMolecule iMol = sim.box.getMoleculeList().getMolecule(i<2 ? i : i+(nA-2));
            p.E(pos.position(iMol));
            for (int j=0; j<3; j++) {
                System.out.println(x0[i*3+j]+" => "+p.getX(j)+"    "+(p.getX(j)-x0[i*3+j]));
                xf[i*3+j] = p.getX(j);
                double dx = xf[i*3+j] - x0[i*3+j];
                disp += dx*dx;
            }
        }
        disp = Math.sqrt(disp);
        System.out.println("disp "+disp);
        double newLatticeEnergy = meterPotentialEnergy.getDataAsScalar();
        System.out.println("Lattice Energy (per molecule): "+newLatticeEnergy/numMolecule);
        for (int l=0; l<201; l++) {
            for (int i=0; i<4; i++) {
                for (int k=i%2; k<numMolecule; k+=2) {
                    boolean isA = (k/nA)%2 == 0;
                    if ((i<2 && !isA) || (i>1 && isA)) continue;
                    IMolecule iMol = sim.box.getMoleculeList().getMolecule(k);
                    p.E(pos.position(iMol));
                    for (int j=0; j<3; j++) {
                        if (l==0) {
                            p.setX(j, p.getX(j)-(xf[i*3+j] - x0[i*3+j])*10);
                        }
                        else {
                            p.setX(j, p.getX(j)+(xf[i*3+j] - x0[i*3+j])/10.0);
                        }
                    }
                    translator.setDestination(p);
                    translator.actionPerformed(iMol);
                }
            }
            double u = 0;
            for (int i=0; i<4; i++) {
                for (int k=(4*nA)+i%2; k<numMolecule; k+=2) {
                    boolean isA = (k/nA)%2 == 0;
                    if ((i<2 && !isA) || (i>1 && isA)) continue;
                    IMolecule iMol = sim.box.getMoleculeList().getMolecule(k);
                    meterPotentialEnergy.setTarget(iMol);
                    u += meterPotentialEnergy.getDataAsScalar();
                    break;
                }
            }
//            meterPotentialEnergy.setTarget((IMolecule)null);
//            System.out.println((l-100)/100.0+" "+(meterPotentialEnergy.getDataAsScalar()-newLatticeEnergy)/numMolecule);
            System.out.println((l-100)/100.0+" "+(u/8-newLatticeEnergy/numMolecule));
        }            
        System.exit(1);

//		FileWriter fileWriter;
//		try{
//        	fileWriter = new FileWriter(filename+"_energy");
//        	
//        } catch (IOException e){
//        	fileWriter = null;
//        }
//        final FileWriter fileWriterEnergy = fileWriter;
//        
//        
//        IAction outputAction = new IAction(){
//        	public void actionPerformed(){
//        		double pe = sim.meterPE.getDataAsScalar();
//        		double he = sim.meterHarm.getDataAsScalar();
//        		System.out.println("(pe-latticeEnergy): " + (pe-latticeEnergy));
////        		try {
////        			fileWriterEnergy.write((pe-latticeEnergy) + " " + he +"\n");
////        			
////        		} catch (IOException e){
////        			
////        		}
//        		
//        	}
//        };
//        
//        IntegratorListenerAction outputActionListener = new IntegratorListenerAction(outputAction);
//        outputActionListener.setInterval((int)simSteps/400);
//        sim.integrator.getEventManager().addListener(outputActionListener);
		
//		MeterNormalizedCoordBeta meterCoord = new MeterNormalizedCoordBeta(sim.box, sim.coordinateDef, sim.species);
//		IntegratorListenerAction meterCoordListener = new IntegratorListenerAction(meterCoord);
//		meterCoordListener.setInterval(1000);                                      
//		sim.integrator.getEventManager().addListener(meterCoordListener);       
		
//		AccumulatorAverage energyAverage = new AccumulatorAverageFixed();
//		DataPump energyPump = new DataPump(meterPotentialEnergy, energyAverage);
//		
//		IntegratorListenerAction energyListener = new IntegratorListenerAction(energyPump);
//		energyListener.setInterval(numMolecule);
//		sim.integrator.getEventManager().addListener(energyListener);
		
//		MeterPressureMolecular meterPressure = new MeterPressureMolecular(sim.space);
//		meterPressure.setIntegrator(sim.integrator);
//						
//		AccumulatorAverage pressureAverage = new AccumulatorAverageCollapsing();
//		DataPump pressurePump = new DataPump(meterPressure, pressureAverage);
//		IntegratorListenerAction pressureListener = new IntegratorListenerAction(pressurePump);
//		pressureListener.setInterval((int)simSteps/100);
//		sim.integrator.getEventManager().addListener(pressureListener);
			
//		double staticPressure = meterPressure.getDataAsScalar();
//		System.out.println("Static Pressure (GPa): " + Pascal.UNIT.fromSim(staticPressure)/1e9);
		
		
		if(true){
			SimulationGraphic simGraphic = new SimulationGraphic(sim, sim.space, sim.getController());
		    simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(10));
		    simGraphic.makeAndDisplayFrame("Beta-Phase Nitrogen Crystal Structure");
		    
		    DiameterHashByType diameter = new DiameterHashByType(sim);
			diameter.setDiameter(sim.species.getNitrogenType(), 3.1);
			diameter.setDiameter(sim.species.getPType(), 0.0);
			
			simGraphic.getDisplayBox(sim.box).setDiameterHash(diameter);
			
			IAction output = new IAction(){

				public void actionPerformed() {
					System.out.println("energy: " + (meterPotentialEnergy.getDataAsScalar()-latticeEnergy));
					
				}
				
			};
			
			IntegratorListenerAction outListener = new IntegratorListenerAction(output);
			outListener.setInterval(500);
			sim.integrator.getEventManager().addListener(outListener);
			
			return;
		}
		
		sim.activityIntegrate.setMaxSteps(simSteps/5);
		sim.getController().actionPerformed();
		System.out.println("****System Equilibrated (20% of SimSteps)****");
		
		long startTime = System.currentTimeMillis();
		System.out.println("\nStart Time: " + startTime);
		sim.getController().reset();

		
		sim.activityIntegrate.setMaxSteps(simSteps);
		sim.getController().actionPerformed();
		
//		try{
//			fileWriterEnergy.close();
//	        
//		} catch (IOException e){
//	        	
//		}
//		double averageEnergy = ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
//		double errorEnergy = ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index);
		
//		double averagePressure = ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
//		double errorPressure = ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index);
		
//		System.out.println("Average energy (per molecule): "   + averageEnergy/numMolecule  
//				+ " ;error: " + errorEnergy/numMolecule);
//		System.out.println("Average pressure (GPa): " + Pascal.UNIT.fromSim(averagePressure)/1e9 
//				+ " ;error: " + Pascal.UNIT.fromSim(errorPressure)/1e9);
		
//		double[] u = sim.coordinateDef.calcU(sim.box.getMoleculeList());
//		
//		for (int i=0; i<u.length; i++){
//			
//			System.out.print (u[i] + ", ");
//			if(i>1 && i%20==19){
//				System.out.println("");
//			}
//		}
		
		
//		double  a = sim.box.getBoundary().getEdgeVector(0).getX(0)/nC0;
//		double  c = sim.box.getBoundary().getEdgeVector(2).getX(2)/nC2;
//		System.out.println("\na: " +a + " ;c: "+c +" ;c/a: " + (c/a));
//		double scaling = sim.box.getBoundary().getEdgeVector(0).getX(0)/(nC0*3.854);
//	    System.out.println("scaling: " + scaling);
//	    long endTime = System.currentTimeMillis();
//		System.out.println("End Time: " + endTime);
//		System.out.println("Time taken: " + (endTime - startTime));
			
//		meterCoord.writeUdistribution(filename);
			
	}
	
    protected static void doTransform(IMolecule molecule, IAtomPositionDefinition posDef, Tensor rotationTensor) {
        IAtomList childList = molecule.getChildList();
        IVector com = posDef.position(molecule);
        for (int iChild = 0; iChild<childList.getAtomCount(); iChild++) {
            IAtom a = childList.getAtom(iChild);
            IVectorMutable r = a.getPosition();
            r.ME(com);
            rotationTensor.transform(r);
            r.PE(com);
        }
    }


	protected Box box;
	protected ISpace space;
	protected PotentialMasterListMolecular potentialMaster;
	protected IntegratorMC integrator;
	protected ActivityIntegrate activityIntegrate;
	protected P2Nitrogen potential;
	protected CoordinateDefinitionNitrogen coordinateDef;
	protected Primitive primitive;
	protected SpeciesN2 species;
	protected MeterHarmonicEnergy meterHarm;
	protected MeterPotentialEnergy meterPE;
	private static final long serialVersionUID = 1L;
}
