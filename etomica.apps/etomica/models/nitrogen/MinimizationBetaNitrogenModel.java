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
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.data.meter.MeterPotentialEnergy;
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
 * Lattice-Energy Minimization for Beta-Nitrogen Structure
 * using gradient
 *  
 * @author Tai Boon Tan
 *
 */
public class MinimizationBetaNitrogenModel extends Simulation{

	public MinimizationBetaNitrogenModel(ISpace space, int[] nC, double density) {
		super(space);
		this.space = space;
		
		BoxAgentSourceCellManagerListMolecular boxAgentSource = new BoxAgentSourceCellManagerListMolecular(this, null, space);
	    BoxAgentManager boxAgentManager = new BoxAgentManager(boxAgentSource);
	     
		double ratio = 1.631;
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
			
			int numDOF = coordinateDef.getCoordinateDim();
			double[] newU = new double[numDOF];
			
//			double[] deviation = new double[]{
//					-1.898068924113261E-4, 9.86615834932536E-5, 1.1427886833814682E-4, -0.0, -0.0, 
//					1.7622553241736227E-4, -8.399450130980313E-6, 5.276642295370948E-5, 0.0, -0.0, 
//					-1.7719627196299825E-4, -6.038293356702695E-5, -1.3280775950974544E-4, -0.0, -0.0, 
//					1.9077763197827835E-4, -2.9879199817450797E-5, -3.423753185316514E-5, -1.7157581151746692E-9, 2.1003461237145532E-8
//			};
//			
//			for(int i=0; i<u.length; i++){
//				u[i] += deviation[i];
//				System.out.print(u[i]+", ");
//				if(i%5==4) System.out.println();
//			}
//			
//			System.exit(1);
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
		
	}
	
	public static void main (String[] args){

		int nC0 = 8; 
		int nC1 = 8; 
		int nC2 = 8;
		int[] nC = new int []{nC0,nC1,nC2};
		int numMolecule =nC[0]*nC[1]*nC[2]*2;

		double density = 0.0230;
		final MinimizationBetaNitrogenModel sim = new MinimizationBetaNitrogenModel(Space3D.getInstance(3), nC, density);
	    
		final MeterPotentialEnergy meterPotentialEnergy = new MeterPotentialEnergy(sim.potentialMaster);
		meterPotentialEnergy.setBox(sim.box);
		double latticeEnergy = meterPotentialEnergy.getDataAsScalar();
		
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
        for (int outer = 0; outer < 400; outer++) {
            System.out.println("**** "+outer+" ****");
	        double totalD = 0;
	        double step = 0;
	        double radianFac = 2e3; //????

	        RotationTensor3D rTensor = new RotationTensor3D();

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
//		            System.out.println(torques[i]);
		            
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
        
        double[] newU = sim.coordinateDef.calcU(sim.box.getMoleculeList());
        for(int i=0; i<4; i++){
        	int iMolec = i<2 ? i : i+(nA-2);
        	for (int k=0; k<5; k++){
        		System.out.print(newU[iMolec*5+k]+", ");
        	}
        	System.out.println();
        }
        
        
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
        System.out.println("Old Lattice Energy (per molecule): "+latticeEnergy/numMolecule);
        System.out.println("New Lattice Energy (per molecule): "+newLatticeEnergy/numMolecule);
        
        System.exit(1);
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
	protected P2Nitrogen potential;
	protected CoordinateDefinitionNitrogen coordinateDef;
	protected Primitive primitive;
	protected SpeciesN2 species;
	protected MeterHarmonicEnergy meterHarm;
	protected MeterPotentialEnergy meterPE;
	private static final long serialVersionUID = 1L;
}
