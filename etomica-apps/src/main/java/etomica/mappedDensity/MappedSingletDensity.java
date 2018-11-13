/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedDensity;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterProfileByVolume;
import etomica.data.types.DataFunction;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.math.function.FunctionDifferentiable;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.P2SoftSphericalTruncatedForceShifted;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.util.List;

/**
 * Simple Simulation for Monte Carlo simulation of Lennard-Jonesium.  Interactions are
 * tracked with cell lists.
 *
 * @author Andrew Schultz
 */

public class MappedSingletDensity extends Simulation {

    public final PotentialMasterCell potentialMaster;
    public final ActivityIntegrate activityIntegrate;
    public final IntegratorMC integrator;

    /**
     * Creates simulation with default parameters from {@link SimParams}
     */
    public MappedSingletDensity() {
        this(new SimParams());
    }

    /**
     * Creates simulation with the given parameters
     */
    public MappedSingletDensity(SimParams params) {
        super(Space3D.getInstance());

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        Box box = new Box(space);
        addBox(box);
        box.setNMolecules(species, params.numAtoms);
        BoxInflate inflater = new BoxInflate(box, space, params.density);
        inflater.actionPerformed();

        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);

        double rc = 3;
        potentialMaster = new PotentialMasterCell(this, rc, space);

        P2LennardJones p2lj = new P2LennardJones(space);
        P2SoftSphericalTruncatedForceShifted p2 = new P2SoftSphericalTruncatedForceShifted(space, p2lj, rc);
        AtomType atomType = species.getLeafType();
        potentialMaster.addPotential(p2, new AtomType[]{atomType, atomType});
//SINE IS ACTUALLY LN(2+SIN)
         if (params.field==Field.SINE || params.field==Field.PHISINEPSINESUM ) {       //AS FOR PHISINEPSINESUM EXTERNAL POT IS LN(2+SIN)-SAME AS FIELD.SINE
            P1Sine p1 = new P1Sine(space, 5, params.temperature);
            potentialMaster.addPotential(p1, new AtomType[]{atomType});
        } else if (params.field == Field.PARABOLIC || params.field == Field.UNIFORM || params.field == Field.PHIPARABOLICPSUMOFGAUSSIANS || params.field==Field.PHIARGPT5PARABOLICPSINESUM) {
            P1Parabolic p1 = new P1Parabolic(space);
            potentialMaster.addPotential(p1, new AtomType[]{atomType});
        } else if (params.field == Field.LNPARABOLIC  || params.field==Field.PHILNPARABOLICPFOURIERSUM) {
             P1Lnparabolic p1 = new P1Lnparabolic(space, params.temperature);
             potentialMaster.addPotential(p1, new AtomType[]{atomType});
         } else if (params.field == Field.EXPMINUSZSQ) {
             P1EXP p1 = new P1EXP(space);
             potentialMaster.addPotential(p1, new AtomType[]{atomType});
         }

        integrator = new IntegratorMC(this, potentialMaster, box);
        integrator.setTemperature(params.temperature);

        MCMoveAtom mcMoveAtom = new MCMoveAtom(random, potentialMaster, space);
        integrator.getMoveManager().addMCMove(mcMoveAtom);

        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        potentialMaster.setCellRange(2);
        potentialMaster.reset();
        integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());
    }

    public static void main(String[] args) {
        SimParams params = new SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
       //     params.field = Field.SINE;
      //         params.field = Field.PARABOLIC;
        //    params.field = Field.PHIPARABOLICPSUMOFGAUSSIANS;
        //       params.field = Field.PHISINEPSINESUM;
       //     params.field=Field.PHIARGPT5PARABOLICPSINESUM;
       //     params.field = Field.LNPARABOLIC;
        //    params.field = Field.PHILNPARABOLICPFOURIERSUM;
            params.field = Field.UNIFORM;
            //   params.field = Field.EXPMINUSZSQ;

            // modify parameters here for interactive testing
        }

        double aa=params.density;
        double bb=0.0;
        double cc=0.0;
        if(params.temperature==5.0 && params.density==0.125){   bb=0.04835;  cc=0.002934;}
        if(params.temperature==5.0 && params.density==0.25){  bb=0.07508;  cc=0.01126;}
        if(params.temperature==5.0 && params.density==0.375){  bb=0.08929;  cc=0.02323;}
        if(params.temperature==5.0 && params.density==0.5){  bb=0.09321;  cc=0.03947;}
        if(params.temperature==5.0 && params.density==0.625){  bb=0.08904;  cc=0.05789;}
        if(params.temperature==5.0 && params.density==0.75){  bb=0.08278;  cc=0.07673;}
        if(params.temperature==5.0 && params.density==0.875){   bb=0.06952;  cc=0.09625;}
        if(params.temperature==5.0 && params.density==1.0){   bb=0.06052;  cc=0.1058;}

        double a1=0.0;
        double a2=0.0;
        double a3=0.0;
        double a4=0.0;
        double a5=0.0;
        double b1=0.0;
        double b2=0.0;
        double b3=0.0;
        double b4=0.0;
        double b5=0.0;
        double c1=1.0;
        double c2=1.0;
        double c3=1.0;
        double c4=1.0;
        double c5=1.0;
        //for arg1
 //       if(params.temperature==5.0 && params.density==0.125){ a1	=	0.237;a2	=	0.04844;b1	=	-0.9846;b2	=	2.204;c1	=	2.073;c2	=	1.413;}
 //       if(params.temperature==5.0 && params.density==0.25){ a1	=	1.23;a2	=	-1.058;b1	=	3.183;b2	=	3.74;c1	=	3.066;c2	=	2.701;}
 //       if(params.temperature==5.0 && params.density==0.375){ a1	=	0.1553;a2	=	0.4751;b1	=	3.15;b2	=	1.314;c1	=	1.369; c2	=	2.01;}
 //       if(params.temperature==5.0 && params.density==0.5){ a1	=	0.2058;a2	=	0.5381;b1	=	3.427;b2	=	1.408;c1	=	-1.447;c2	=	2.12;}
 //       if(params.temperature==5.0 && params.density==0.625){ a1	=	2.135;a2	=	-1.462;a3	=	0.06693;b1	=	2.062;b2	=	2.088;b3	=	4.14;c1	=	1.95;c2	=	1.687;c3	=	0.4266;}
 //       if(params.temperature==5.0 && params.density==0.75){ a1	=	1.071;a2	=	-0.4528;a3	=	0.6855;a4	=	0.1572;b1	=	0.9675;b2	=	1.064;b3	=	2.958;b4	=	4.015;c1	=	1.103;c2	=	0.8714;c3	=	1.65;c4	=	0.4183;}
 //       if(params.temperature==5.0 && params.density==0.875){ a1	=	0.6406; a2	=	0.6802; a3	=	0.2837; a4	=	-0.1293; a5	=	0.7941; b1	=	3.841; b2	=	2.669; b3	=	-0.08254; b4	=	2.344; b5	=	1.085; c1	=	0.5945; c2	=	0.9332; c3	=	-0.02224; c4	=	0.3835; c5	=	1.486;}
 //       if(params.temperature==5.0 && params.density==1.0){ a1	=	2.628;a2	=	-0.4551;a3	=	-1.074;a4	=	-1.352;a5	=	-0.532;b1	=	1.964;b2	=	3.103;b3	=	1.431;b4	=	2.313;b5	=	0.565;c1	=	1.803;c2	=	0.4203;c3	=	0.6687;c4	=	0.7064;c5	=	0.6801;}

        //for arg0.8
        if(params.temperature==5.0 && params.density==0.125){ a1	=	0.2707;a2	=	-0.01455;b1	=	1.355;b2	=	-1.56;c1	=	2.13;c2	=	1.176;}
        if(params.temperature==5.0 && params.density==0.25){ a1	=	0.3203;a2	=	0.2001;b1	=	2.199;b2	=	0.6374;c1	=	1.901;c2	=	1.399;}
        if(params.temperature==5.0 && params.density==0.375){ a1	=	0.004659;a2	=	0.1673;a3	=	0.4294;b1	=	0.4074;b2	=	3.319;b3	=	1.328;c1	=	0.4601;c2	=	1.604;c3	=	2.103;}
        if(params.temperature==5.0 && params.density==0.5){ a1	=	0.1129;a2	=	-1.727;a3	=	1.778;b1	=	4.089;b2	=	0.822;b3	=	-0.1271;c1	=	1.263;c2	=	1.937;c3	=	2.624;}
        if(params.temperature==5.0 && params.density==0.625){ a1	=	0.5895;a2	=	-0.01595;a3	=	0.3503;a4	=	0.3237;a5	=	0.09432;b1	=	2.874;b2	=	1.802;b3	=	0.4146;b4	=	1.358;b5	=	4.243;c1	=	1.913;c2	=	0.337;c3	=	0.757;c4	=	0.9766;c5	=	0.4991;}
        if(params.temperature==5.0 && params.density==0.75){ a1	=	-2.324;a2	=	-0.7592;a3	=	0.8775;a4	=	0.4945;a5	=	0.1048;b1	=	5.917;b2	=	-7.58;b3	=	3.855;b4	=	1.036;b5	=	4.005;c1	=	-0.3608;c2	=	4.357;c3	=	2.436;c4	=	2.007;c5	=	0.3939;}
        if(params.temperature==5.0 && params.density==0.875){ a1	=	0.3637;a2	=	0.7177;a3	=	0.3241;a4	=	0.8154;a5	=	0.6729;b1	=	-0.2692;b2	=	1.052;b3	=	2.819;b4	=	3.741;b5	=	2.129;c1	=	0.6688;c2	=	0.8477;c3	=	0.5114;c4	=	0.8134;c5	=	0.7817;}
        if(params.temperature==5.0 && params.density==1.0){ a1	=	-2.01;a2	=	-1.269;a3	=	1.066;a4	=	3.727;a5	=	-2.288;b1	=	2.738;b2	=	1.901;b3	=	2.705;b4	=	1.631;b5	=	0.9388;c1	=	0.6194;c2	=	-0.839;c3	=	0.4916;c4	=	1.91;c5	=	-1.224;}


            //  if(params.temperature==5.0 && params.density==0.125){a1=1.0; a2=-0.4529; b1=0.809; b2=0.8364; c1=0.8498; c2=0.6394;}
      //  if(params.temperature==5.0 && params.density==0.25){ a1	=	0.2884;a2	=	0.686;b1	=	1.795;b2	=	-0.7539;c1	=	0.6091;c2	=	-1.041;}
      //  if(params.temperature==5.0 && params.density==0.375){ a1	=	0.8073;a2	=	0.3485;b1	=	0.8895;b2	=	2.079;c1	=	-1.178;c2	=	0.6238;}
      //  if(params.temperature==5.0 && params.density==0.5){ a1	=	-19.71;a2	=	0.3829;a3	=	-0.006994;a4	=	20.64;b1	=	1.069;b2	=	2.489;b3	=	-0.9805;b4	=	1.069;c1	=	0.902;c2	=	0.5318;c3	=	0.0214;c4	=	0.9197;}
      //  if(params.temperature==5.0 && params.density==0.625){ a1	=	0.6381;a2	=	0.2972;a3	=	-2.397;a4	=	2.915;b1	=	1.921;b2	=	-2.74;b3	=	0.006159;b4	=	0;c1	=	1.064;c2	=	0.4744;c3	=	2.086;c4	=	1.994;}
      //  if(params.temperature==5.0 && params.density==0.75){ a1	=	0.5239;a2	=	0.9187;a3	=	0.1986;a4	=	0.8576;b1	=	2.97;b2	=	0.6462;b3	=	1.253;b4	=	2.053;c1	=	0.5098;c2	=	0.9231;c3	=	0.4761;c4	=	0.7447;}
      //  if(params.temperature==5.0 && params.density==0.875){ a1	=	-0.3673;a2	=	1.452;a3	=	-8.242;a4	=	8.938;b1	=	-0.9993;b2	=	1.055;b3	=	2.766;b4	=	2.77;c1	=	0.5828;c2	=	1.123;c3	=	0.5688;c4	=	0.5953;}
      //  if(params.temperature==5.0 && params.density==1.0){ a1	=	0.6222;a2	=	1.174;a3	=	0.9962;a4	=	-0.283;b1	=	3.391;b2	=	2.247;b3	=	0.7525;b4	=	2.079;c1	=	0.4002;c2	=	0.9006;c3	=	1.061;c4	=	0.4084;}

//        if(params.temperature==5.0 && params.density==0.125){ a1	=	0.5703;a2	=	0.5239;b1	=	0.918;b2	=	-0.2593;c1	=	0.4308;c2	=	0.6005;}
//        if(params.temperature==5.0 && params.density==0.25){ a1	=	1.002;a2	=	0.8221;b1	=	-0.4634;b2	=	1.292;c1	=	0.5549;c2	=	0.4054;}
//        if(params.temperature==5.0 && params.density==0.375){ a1	=	0.2228;a2	=	0.7951;a3	=	1.152;b1	=	-0.001652;b2	=	1.627;b3	=	0.7323;c1	=	0.3568;c2	=	0.362;c3	=	-0.6928;}
//        if(params.temperature==5.0 && params.density==0.5) { a1	=	1.243;a2	=	0.931;a3	=	0.5515;b1	=	0.9433;b2	=	1.859;b3	=	-0.143;c1	=	0.6608;c2	=	0.3651;c3	=	0.4552; }
//        if(params.temperature==5.0 && params.density==0.625){ a1	=	1.325;a2	=	1.144;a3	=	1.352;b1	=	0.374;b2	=	2.06;b3	=	1.21;c1	=	0.4171;c2	=	0.3602;c3	=	0.4988; }
//        if(params.temperature==5.0 && params.density==0.75) { a1	=	0.6231;a2	=	-0.5084;a3	=	1.04;a4	=	1.449;b1	=	0.003771;b2	=	1.066;b3	=	2.288;b4	=	1.267;c1	=	1.021;c2	=	0.3388;c3	=	0.3343;c4	=	-0.7168; }
//        if(params.temperature==5.0 && params.density==0.875) { a1	=	1.315;a2	=	31.42;a3	=	1.473;a4	=	-36.17;b1	=	2.429;b2	=	0.0000438;b3	=	1.616;b4	=	0.2096;c1	=	0.3507;c2	=	0.6199;c3	=	0.4191;c4	=	0.5117; }
//ADD AT DENSITY 1

//FOURIER SUM
        double w=0.0;
        double aa0=0.0;
        double aa1=0.0;
        double aa2=0.0;
        double aa3=0.0;
        double aa4=0.0;
        double aa5=0.0;
        double aa6=0.0;
        double aa7=0.0;
        double aa8=0.0;
        double bb1=0.0;
        double bb2=0.0;
        double bb3=0.0;
        double bb4=0.0;
        double bb5=0.0;
        double bb6=0.0;
        double bb7=0.0;
        double bb8=0.0;
        if(params.temperature==5.0 && params.density==0.125) { aa0	=	0.124700000;aa1	=	-0.004748000;bb1	=	-0.000047860;aa2	=	-0.003962000;bb2	=	0.000066840;aa3	=	-0.002992000;bb3	=	0.000004400;aa4	=	-0.002137000;bb4	=	0.000073050;aa5	=	-0.001254000;bb5	=	-0.000033180;aa6	=	-0.000778400;bb6	=	-0.000027450;aa7	=	-0.000433600;bb7	=	0.000026760;aa8	=	-0.000246500;bb8	=	0.000043330;w	=	0.448500000; }
        if(params.temperature==5.0 && params.density==0.25) { aa0	=	0.250100000;aa1	=	-0.007282000;bb1	=	-0.000022140;aa2	=	-0.006333000;bb2	=	-0.000147300;aa3	=	-0.004818000;bb3	=	0.000255200;aa4	=	-0.002671000;bb4	=	-0.000231700;aa5	=	-0.001823000;bb5	=	-0.000345900;aa6	=	-0.000950700;bb6	=	-0.000041880;aa7	=	-0.000389300;bb7	=	0.000185600;aa8	=	-0.000304900;bb8	=	-0.000273400;w	=	0.491400000; }
        if(params.temperature==5.0 && params.density==0.375) { aa0	=	0.375000000;aa1	=	-0.008432000;bb1	=	0.001360000;aa2	=	-0.006944000;bb2	=	-0.000718100;aa3	=	-0.005082000;bb3	=	0.000411600;aa4	=	-0.002238000;bb4	=	-0.000439300;aa5	=	-0.001487000;bb5	=	0.000502500;aa6	=	-0.000701700;bb6	=	0.000325200;aa7	=	0.000003404;bb7	=	0.000097250;aa8	=	0.000231700;bb8	=	0.000142900;w	=	0.565900000; }
        if(params.temperature==5.0 && params.density==0.5) { aa0	=	0.500100000;aa1	=	-0.007932000;bb1	=	0.001553000;aa2	=	-0.004792000;bb2	=	0.000292000;aa3	=	-0.003789000;bb3	=	-0.000185600;aa4	=	-0.002862000;bb4	=	0.000376800;aa5	=	-0.001774000;bb5	=	-0.000787800;aa6	=	-0.000795400;bb6	=	-0.000007208;aa7	=	-0.000139200;bb7	=	0.000006755;aa8	=	0.000238400;bb8	=	0.000078790;w	=	0.618300000; }
        if(params.temperature==5.0 && params.density==0.625) { aa0	=	0.625000000;aa1	=	-0.006388000;bb1	=	0.001690000;aa2	=	-0.005215000;bb2	=	0.000550800;aa3	=	-0.002219000;bb3	=	0.000579800;aa4	=	-0.003402000;bb4	=	0.000770600;aa5	=	-0.001442000;bb5	=	-0.000147100;aa6	=	-0.001351000;bb6	=	-0.000733100;aa7	=	0.000171400;bb7	=	-0.000578300;aa8	=	0.000080490;bb8	=	0.000154300;w	=	0.658500000; }
        if(params.temperature==5.0 && params.density==0.75) { aa0	=	0.750000000;aa1	=	-0.007309000;bb1	=	-0.003973000;aa2	=	-0.006042000;bb2	=	-0.001989000;aa3	=	-0.000348100;bb3	=	-0.002074000;aa4	=	-0.002198000;bb4	=	-0.000535500;aa5	=	-0.000049120;bb5	=	-0.000008381;aa6	=	0.001013000;bb6	=	0.001290000;aa7	=	-0.000416200;bb7	=	-0.000802300;w	=	0.741200000; }
        if(params.temperature==5.0 && params.density==0.875) { aa0	=	0.874800000;aa1	=	-0.005123000;bb1	=	-0.013320000;aa2	=	-0.006292000;bb2	=	0.000337900;aa3	=	-0.001784000;bb3	=	0.000150300;aa4	=	-0.001172000;bb4	=	0.001691000;aa5	=	0.000588800;bb5	=	0.000881900;aa6	=	-0.001800000;bb6	=	-0.000843700;aa7	=	0.002573000;bb7	=	-0.000501800;aa8	=	0.000766900;bb8	=	0.000112700;w	=	0.711500000; }
        if(params.temperature==5.0 && params.density==1.0) { aa0	=	1.0;aa1	=	-0.016120000;bb1	=	0.004948000;aa2	=	0.000768600;bb2	=	0.005620000;aa3	=	-0.007674000;bb3	=	0.002280000;aa4	=	-0.001069000;bb4	=	-0.002247000;aa5	=	0.000506000;bb5	=	-0.001660000;aa6	=	0.002331000;bb6	=	-0.000510100;aa7	=	0.001027000;bb7	=	-0.003227000;aa8	=	-0.002271000;bb8	=	0.003299000;w	=	0.776200000; }

        double sa1=0.0;
        double sa2=0.0;
        double sa3=0.0;
        double sa4=0.0;
        double sa5=0.0;
        double sa6=0.0;
        double sa7=0.0;
        double sa8=0.0;
        double sb1=1.0;
        double sb2=1.0;
        double sb3=1.0;
        double sb4=1.0;
        double sb5=1.0;
        double sb6=1.0;
        double sb7=1.0;
        double sb8=1.0;
        double sc1=1.0;
        double sc2=1.0;
        double sc3=1.0;
        double sc4=1.0;
        double sc5=1.0;
        double sc6=1.0;
        double sc7=1.0;
        double sc8=1.0;
        if(params.temperature==5.0 && params.density==0.125) { sa1	=	0.1755;sb1	=	0.6441;sc1	=	1.567;sa2	=	0.1092;sb2	=	2.182;sc2	=	1.571;sa3	=	0.01281;sb3	=	3.772;sc3	=	1.566;sa4	=	0.004448;sb4	=	6.216;sc4	=	-1.567;}
        if(params.temperature==5.0 && params.density==0.25) { sa1	=	0.3907;sb1	=	0.8563;sc1	=	1.574;sa2	=	0.04736;sb2	=	2.826;sc2	=	1.588;sa3	=	-0.006766;sb3	=	5.223;sc3	=	-1.596;sa4	=	-0.0261;sb4	=	4.332;sc4	=	1.555; }
        if(params.temperature==5.0 && params.density==0.375) { sa1	=	0.5497;sb1	=	0.8236;sc1	=	1.569;sa2	=	0.03292;sb2	=	1.549;sc2	=	-1.562;sa3	=	-0.01415;sb3	=	5.14;sc3	=	1.479;sa4	=	-0.01071;sb4	=	5.587;sc4	=	-1.624; }
        if(params.temperature==5.0 && params.density==0.5) { sa1	=	0.7866;sb1	=	0.9087;sc1	=	1.573;sa2	=	0.2075;sb2	=	1.814;sc2	=	-1.578;sa3	=	0.02799;sb3	=	3.635;sc3	=	1.606;sa4	=	0.009727;sb4	=	5.452;sc4	=	-1.621;sa5	=	0.006593;sb5	=	7.264;sc5	=	1.624;sa6	=	0.003764;sb6	=	9.083;sc6	=	-1.199;sa7	=	0.001517;sb7	=	18.17;sc7	=	-1.278;sa8	=	0.00151;sb8	=	10.9;sc8	=	0.9269; }
        if(params.temperature==5.0 && params.density==0.625) { sa1	=	1.036;sb1	=	0.9926;sc1	=	1.573;sa2	=	0.4123;sb2	=	1.964;sc2	=	-1.57;sa3	=	0.1168;sb3	=	3.786;sc3	=	1.562;sa4	=	0.4995;sb4	=	5.521;sc4	=	-1.574;sa5	=	0.4583;sb5	=	5.615;sc5	=	1.566;sa6	=	0.004006;sb6	=	9.054;sc6	=	-1.055;sa7	=	0.002954;sb7	=	16.57;sc7	=	1.23;sa8	=	0.003006;sb8	=	17.94;sc8	=	-1.403; }
        if(params.temperature==5.0 && params.density==0.75) { sa1	=	1.286;sb1	=	0.9526;sc1	=	1.573;sa2	=	0.5388;sb2	=	1.71;sc2	=	-1.569;sa3	=	0.6026;sb3	=	4.409;sc3	=	1.586;sa4	=	0.5634;sb4	=	4.537;sc4	=	-1.551;sa5	=	0.01431;sb5	=	7.245;sc5	=	1.727;sa6	=	0.01098;sb6	=	9.394;sc6	=	-1.428;sa7	=	0.00889;sb7	=	10.82;sc7	=	1.714;sa8	=	0.005029;sb8	=	16.32;sc8	=	1.553; }
        if(params.temperature==5.0 && params.density==0.875) { sa1	=	1.577;sb1	=	0.9674;sc1	=	1.572;sa2	=	0.7353;sb2	=	1.667;sc2	=	-1.582;sa3	=	0.1001;sb3	=	3.699;sc3	=	1.63;sa4	=	0.04359;sb4	=	5.56;sc4	=	-1.623;sa5	=	0.02722;sb5	=	7.299;sc5	=	1.517;sa6	=	0.01905;sb6	=	9.064;sc6	=	-1.756;sa7	=	0.009713;sb7	=	10.88;sc7	=	1.665;sa8	=	0.006004;sb8	=	16.28;sc8	=	1.466; }
        if(params.temperature==5.0 && params.density==1.0) { sa1	=	1.568;sb1	=	0.9074;sc1	=	1.568;sa2	=	0.6369;sb2	=	1.819;sc2	=	-1.575;sa3	=	0.1242;sb3	=	3.63;sc3	=	1.569;sa4	=	0.05096;sb4	=	5.453;sc4	=	-1.598;sa5	=	0.03115;sb5	=	7.265;sc5	=	1.507;sa6	=	0.01709;sb6	=	9.083;sc6	=	-1.713;sa7	=	0.01261;sb7	=	16.35;sc7	=	1.994;sa8	=	0.0115;sb8	=	10.9;sc8	=	1.644; }


            MappedSingletDensity sim = new MappedSingletDensity(params);
        long steps = params.steps;
        int interval = 5* params.numAtoms;
        int blocks = 100;
        long blockSize = steps / (interval * blocks);

        System.out.println("Lennard-Jones Monte Carlo simulation");
        System.out.println("N: " + params.numAtoms);
        System.out.println("T: " + params.temperature);
        System.out.println("density: " + params.density);
        System.out.println("steps: " + params.steps);

        // equilibration
        long t1 = System.currentTimeMillis();
        sim.activityIntegrate.setMaxSteps(steps / 10);
        sim.activityIntegrate.actionPerformed();
        System.out.println("equilibration finished");

        // data collection
        MeterProfileByVolume densityMeter = new MeterProfileByVolume(sim.space);
        densityMeter.setBox(sim.box());
        densityMeter.setProfileDim(2);
        densityMeter.getXDataSource().setNValues(params.bins);  //conventional bins=100
        densityMeter.reset();

        MeterNMolecules meterNMolecules = new MeterNMolecules();
        densityMeter.setDataSource(meterNMolecules);
        AccumulatorAverageFixed acc = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pump = new DataPumpListener(densityMeter, acc, interval);
        sim.getIntegrator().getEventManager().addListener(pump);

        MeterProfileForceSum densityMeterForce = new MeterProfileForceSum(sim.box(), sim.potentialMaster, params.temperature);
        densityMeterForce.setProfileDim(2);
        densityMeterForce.getXDataSource().setNValues(params.bins);  //pap bins=1000
        densityMeterForce.reset();
        AccumulatorAverageFixed accForce = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pumpForce = new DataPumpListener(densityMeterForce, accForce, interval);
        sim.getIntegrator().getEventManager().addListener(pumpForce);

        FunctionDifferentiable f;
        double L = sim.box().getBoundary().getBoxSize().getX(2);
        switch (params.field) {
            case PARABOLIC:
                f = new FunctionParabolic(L, params.temperature);
                break;
            case SINE:
                f = new FunctionSine(5, L);
                break;
            case UNIFORM:
                f = new FunctionUniform(L);
                break;
            case LNPARABOLIC:
                f = new FunctionLnparabolic(L);
                break;
            case PHIPARABOLICPSUMOFGAUSSIANS:
                f = new FunctionPhiparabolicpsumofgaussians(L,a1,b1,c1,a2,b2,c2,a3,b3,c3,a4,b4,c4,a5,b5,c5);
                break;
            case PHISINEPSINESUM:
                f = new FunctionPhisinepsinesum(L,aa,bb,cc);
                break;
            case EXPMINUSZSQ:
                f = new FunctionExpminuszsq(L,aa0,aa1,aa2,aa3,aa4,aa5,aa6,aa7,aa8,bb1,bb2,bb3,bb4,bb5,bb6,bb7,bb8,w);
                break;
            case PHIARGPT5PARABOLICPSINESUM:
                f = new FunctionPhiargpt5parabolicpsinesum(L,sa1,sb1,sc1,sa2,sb2,sc2,sa3,sb3,sc3,sa4,sb4,sc4,sa5,sb5,sc5,sa6,sb6,sc6,sa7,sb7,sc7,sa8,sb8,sc8);
                break;
            default:
                throw new RuntimeException("not yet");
        }
        MeterProfileMappedAvg densityMeterMappedAvg = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
        densityMeterMappedAvg.setProfileDim(2);
        densityMeterMappedAvg.getXDataSource().setNValues(params.bins);  //mappedavg bins=1000
        densityMeterMappedAvg.reset();
        AccumulatorAverageFixed accMappedAvg = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pumpMappedAvg = new DataPumpListener(densityMeterMappedAvg, accMappedAvg, 5 * params.numAtoms);
        sim.getIntegrator().getEventManager().addListener(pumpMappedAvg);

        sim.activityIntegrate.setMaxSteps(steps);
        sim.getIntegrator().resetStepCount();
        sim.integrator.getMoveManager().setEquilibrating(false);
        sim.activityIntegrate.actionPerformed();

        long t2 = System.currentTimeMillis();
        IData data =  acc.getData(acc.AVERAGE);
        IData dataunc =  acc.getData(acc.ERROR);
        IData dataForce =  accForce.getData(accForce.AVERAGE);
        IData dataForceunc =  accForce.getData(accForce.ERROR);
        IData dataMappedAvg =  accMappedAvg.getData(accMappedAvg.AVERAGE);
        IData dataMappedAvgunc =  accMappedAvg.getData(accMappedAvg.ERROR);

       IData zdata= ((DataFunction.DataInfoFunction)((DataGroup.DataInfoGroup)acc.getDataInfo()).getSubDataInfo(0)).getXDataSource().getIndependentData(0);
       for (int i=0;i<zdata.getLength();i++){
           System.out.println(zdata.getValue(i)+" "+data.getValue(i)+" "+dataForce.getValue(i)+" "+dataMappedAvg.getValue(i)+" "+dataunc.getValue(i)+" "+dataForceunc.getValue(i)+" "+dataMappedAvgunc.getValue(i));
       }
        System.out.println("time: " + (t2 - t1) * 0.001);
    }

    public static class Graphic {
        public static void main(String[] args) {
            SimParams params = new SimParams();

            if (args.length > 0) {
                ParseArgs.doParseArgs(params, args);
            }
            else {
            //    params.field = Field.SINE;
           //         params.field = Field.PARABOLIC;
            //    params.field = Field.PHIPARABOLICPSUMOFGAUSSIANS;
           //      params.field = Field.PHISINEPSINESUM;
              params.field = Field.UNIFORM;
          //      params.field=Field.PHIARGPT5PARABOLICPSINESUM;
                //      params.field = Field.LNPARABOLIC;
         //       params.field = Field.PHILNPARABOLICPFOURIERSUM;
             //   params.field = Field.EXPMINUSZSQ;

                // modify parameters here for interactive testing
            }


            double aa=params.density;
            double bb=0.0;
            double cc=0.0;
            if(params.temperature==5.0 && params.density==0.125){   bb=0.04835;  cc=0.002934;}
            if(params.temperature==5.0 && params.density==0.25){  bb=0.07508;  cc=0.01126;}
            if(params.temperature==5.0 && params.density==0.375){  bb=0.08929;  cc=0.02323;}
            if(params.temperature==5.0 && params.density==0.5){  bb=0.09321;  cc=0.03947;}
            if(params.temperature==5.0 && params.density==0.625){  bb=0.08904;  cc=0.05789;}
            if(params.temperature==5.0 && params.density==0.75){  bb=0.08278;  cc=0.07673;}
            if(params.temperature==5.0 && params.density==0.875){   bb=0.06952;  cc=0.09625;}
            if(params.temperature==5.0 && params.density==1.0){   bb=0.06052;  cc=0.1058;}

            double a1=0.0;
            double a2=0.0;
            double a3=0.0;
            double a4=0.0;
            double a5=0.0;
            double b1=0.0;
            double b2=0.0;
            double b3=0.0;
            double b4=0.0;
            double b5=0.0;
            double c1=1.0;
            double c2=1.0;
            double c3=1.0;
            double c4=1.0;
            double c5=1.0;
            //for arg1
     //       if(params.temperature==5.0 && params.density==0.125){ a1	=	0.237;a2	=	0.04844;b1	=	-0.9846;b2	=	2.204;c1	=	2.073;c2	=	1.413;}
     //       if(params.temperature==5.0 && params.density==0.25){ a1	=	1.23;a2	=	-1.058;b1	=	3.183;b2	=	3.74;c1	=	3.066;c2	=	2.701;}
     //       if(params.temperature==5.0 && params.density==0.375){ a1	=	0.1553;a2	=	0.4751;b1	=	3.15;b2	=	1.314;c1	=	1.369; c2	=	2.01;}
     //       if(params.temperature==5.0 && params.density==0.5){ a1	=	0.2058;a2	=	0.5381;b1	=	3.427;b2	=	1.408;c1	=	-1.447;c2	=	2.12;}
     //       if(params.temperature==5.0 && params.density==0.625){ a1	=	2.135;a2	=	-1.462;a3	=	0.06693;b1	=	2.062;b2	=	2.088;b3	=	4.14;c1	=	1.95;c2	=	1.687;c3	=	0.4266;}
     //       if(params.temperature==5.0 && params.density==0.75){ a1	=	1.071;a2	=	-0.4528;a3	=	0.6855;a4	=	0.1572;b1	=	0.9675;b2	=	1.064;b3	=	2.958;b4	=	4.015;c1	=	1.103;c2	=	0.8714;c3	=	1.65;c4	=	0.4183;}
     //       if(params.temperature==5.0 && params.density==0.875){ a1	=	0.6406; a2	=	0.6802; a3	=	0.2837; a4	=	-0.1293; a5	=	0.7941; b1	=	3.841; b2	=	2.669; b3	=	-0.08254; b4	=	2.344; b5	=	1.085; c1	=	0.5945; c2	=	0.9332; c3	=	-0.02224; c4	=	0.3835; c5	=	1.486;}
     //       if(params.temperature==5.0 && params.density==1.0){ a1	=	2.628;a2	=	-0.4551;a3	=	-1.074;a4	=	-1.352;a5	=	-0.532;b1	=	1.964;b2	=	3.103;b3	=	1.431;b4	=	2.313;b5	=	0.565;c1	=	1.803;c2	=	0.4203;c3	=	0.6687;c4	=	0.7064;c5	=	0.6801;}

            //for arg0.8
            if(params.temperature==5.0 && params.density==0.125){ a1	=	0.2707;a2	=	-0.01455;b1	=	1.355;b2	=	-1.56;c1	=	2.13;c2	=	1.176;}
            if(params.temperature==5.0 && params.density==0.25){ a1	=	0.3203;a2	=	0.2001;b1	=	2.199;b2	=	0.6374;c1	=	1.901;c2	=	1.399;}
            if(params.temperature==5.0 && params.density==0.375){ a1	=	0.004659;a2	=	0.1673;a3	=	0.4294;b1	=	0.4074;b2	=	3.319;b3	=	1.328;c1	=	0.4601;c2	=	1.604;c3	=	2.103;}
            if(params.temperature==5.0 && params.density==0.5){ a1	=	0.1129;a2	=	-1.727;a3	=	1.778;b1	=	4.089;b2	=	0.822;b3	=	-0.1271;c1	=	1.263;c2	=	1.937;c3	=	2.624;}
            if(params.temperature==5.0 && params.density==0.625){ a1	=	0.5895;a2	=	-0.01595;a3	=	0.3503;a4	=	0.3237;a5	=	0.09432;b1	=	2.874;b2	=	1.802;b3	=	0.4146;b4	=	1.358;b5	=	4.243;c1	=	1.913;c2	=	0.337;c3	=	0.757;c4	=	0.9766;c5	=	0.4991;}
            if(params.temperature==5.0 && params.density==0.75){ a1	=	-2.324;a2	=	-0.7592;a3	=	0.8775;a4	=	0.4945;a5	=	0.1048;b1	=	5.917;b2	=	-7.58;b3	=	3.855;b4	=	1.036;b5	=	4.005;c1	=	-0.3608;c2	=	4.357;c3	=	2.436;c4	=	2.007;c5	=	0.3939;}
            if(params.temperature==5.0 && params.density==0.875){ a1	=	0.3637;a2	=	0.7177;a3	=	0.3241;a4	=	0.8154;a5	=	0.6729;b1	=	-0.2692;b2	=	1.052;b3	=	2.819;b4	=	3.741;b5	=	2.129;c1	=	0.6688;c2	=	0.8477;c3	=	0.5114;c4	=	0.8134;c5	=	0.7817;}
            if(params.temperature==5.0 && params.density==1.0){ a1	=	-2.01;a2	=	-1.269;a3	=	1.066;a4	=	3.727;a5	=	-2.288;b1	=	2.738;b2	=	1.901;b3	=	2.705;b4	=	1.631;b5	=	0.9388;c1	=	0.6194;c2	=	-0.839;c3	=	0.4916;c4	=	1.91;c5	=	-1.224;}


            double w=0.0;
            double aa0=0.0;
            double aa1=0.0;
            double aa2=0.0;
            double aa3=0.0;
            double aa4=0.0;
            double aa5=0.0;
            double aa6=0.0;
            double aa7=0.0;
            double aa8=0.0;
            double bb1=0.0;
            double bb2=0.0;
            double bb3=0.0;
            double bb4=0.0;
            double bb5=0.0;
            double bb6=0.0;
            double bb7=0.0;
            double bb8=0.0;
            if(params.temperature==5.0 && params.density==0.125) { aa0	=	0.124700000;aa1	=	-0.004748000;bb1	=	-0.000047860;aa2	=	-0.003962000;bb2	=	0.000066840;aa3	=	-0.002992000;bb3	=	0.000004400;aa4	=	-0.002137000;bb4	=	0.000073050;aa5	=	-0.001254000;bb5	=	-0.000033180;aa6	=	-0.000778400;bb6	=	-0.000027450;aa7	=	-0.000433600;bb7	=	0.000026760;aa8	=	-0.000246500;bb8	=	0.000043330;w	=	0.448500000; }
            if(params.temperature==5.0 && params.density==0.25) { aa0	=	0.250100000;aa1	=	-0.007282000;bb1	=	-0.000022140;aa2	=	-0.006333000;bb2	=	-0.000147300;aa3	=	-0.004818000;bb3	=	0.000255200;aa4	=	-0.002671000;bb4	=	-0.000231700;aa5	=	-0.001823000;bb5	=	-0.000345900;aa6	=	-0.000950700;bb6	=	-0.000041880;aa7	=	-0.000389300;bb7	=	0.000185600;aa8	=	-0.000304900;bb8	=	-0.000273400;w	=	0.491400000; }
            if(params.temperature==5.0 && params.density==0.375) { aa0	=	0.375000000;aa1	=	-0.008432000;bb1	=	0.001360000;aa2	=	-0.006944000;bb2	=	-0.000718100;aa3	=	-0.005082000;bb3	=	0.000411600;aa4	=	-0.002238000;bb4	=	-0.000439300;aa5	=	-0.001487000;bb5	=	0.000502500;aa6	=	-0.000701700;bb6	=	0.000325200;aa7	=	0.000003404;bb7	=	0.000097250;aa8	=	0.000231700;bb8	=	0.000142900;w	=	0.565900000; }
            if(params.temperature==5.0 && params.density==0.5) { aa0	=	0.500100000;aa1	=	-0.007932000;bb1	=	0.001553000;aa2	=	-0.004792000;bb2	=	0.000292000;aa3	=	-0.003789000;bb3	=	-0.000185600;aa4	=	-0.002862000;bb4	=	0.000376800;aa5	=	-0.001774000;bb5	=	-0.000787800;aa6	=	-0.000795400;bb6	=	-0.000007208;aa7	=	-0.000139200;bb7	=	0.000006755;aa8	=	0.000238400;bb8	=	0.000078790;w	=	0.618300000; }
            if(params.temperature==5.0 && params.density==0.625) { aa0	=	0.625000000;aa1	=	-0.006388000;bb1	=	0.001690000;aa2	=	-0.005215000;bb2	=	0.000550800;aa3	=	-0.002219000;bb3	=	0.000579800;aa4	=	-0.003402000;bb4	=	0.000770600;aa5	=	-0.001442000;bb5	=	-0.000147100;aa6	=	-0.001351000;bb6	=	-0.000733100;aa7	=	0.000171400;bb7	=	-0.000578300;aa8	=	0.000080490;bb8	=	0.000154300;w	=	0.658500000; }
            if(params.temperature==5.0 && params.density==0.75) { aa0	=	0.750000000;aa1	=	-0.007309000;bb1	=	-0.003973000;aa2	=	-0.006042000;bb2	=	-0.001989000;aa3	=	-0.000348100;bb3	=	-0.002074000;aa4	=	-0.002198000;bb4	=	-0.000535500;aa5	=	-0.000049120;bb5	=	-0.000008381;aa6	=	0.001013000;bb6	=	0.001290000;aa7	=	-0.000416200;bb7	=	-0.000802300;w	=	0.741200000; }
            if(params.temperature==5.0 && params.density==0.875) { aa0	=	0.874800000;aa1	=	-0.005123000;bb1	=	-0.013320000;aa2	=	-0.006292000;bb2	=	0.000337900;aa3	=	-0.001784000;bb3	=	0.000150300;aa4	=	-0.001172000;bb4	=	0.001691000;aa5	=	0.000588800;bb5	=	0.000881900;aa6	=	-0.001800000;bb6	=	-0.000843700;aa7	=	0.002573000;bb7	=	-0.000501800;aa8	=	0.000766900;bb8	=	0.000112700;w	=	0.711500000; }
            if(params.temperature==5.0 && params.density==1.0) { aa0	=	1.0;aa1	=	-0.016120000;bb1	=	0.004948000;aa2	=	0.000768600;bb2	=	0.005620000;aa3	=	-0.007674000;bb3	=	0.002280000;aa4	=	-0.001069000;bb4	=	-0.002247000;aa5	=	0.000506000;bb5	=	-0.001660000;aa6	=	0.002331000;bb6	=	-0.000510100;aa7	=	0.001027000;bb7	=	-0.003227000;aa8	=	-0.002271000;bb8	=	0.003299000;w	=	0.776200000; }

            double sa1=0.0;
            double sa2=0.0;
            double sa3=0.0;
            double sa4=0.0;
            double sa5=0.0;
            double sa6=0.0;
            double sa7=0.0;
            double sa8=0.0;
            double sb1=1.0;
            double sb2=1.0;
            double sb3=1.0;
            double sb4=1.0;
            double sb5=1.0;
            double sb6=1.0;
            double sb7=1.0;
            double sb8=1.0;
            double sc1=1.0;
            double sc2=1.0;
            double sc3=1.0;
            double sc4=1.0;
            double sc5=1.0;
            double sc6=1.0;
            double sc7=1.0;
            double sc8=1.0;
            if(params.temperature==5.0 && params.density==0.125) { sa1	=	0.1755;sb1	=	0.6441;sc1	=	1.567;sa2	=	0.1092;sb2	=	2.182;sc2	=	1.571;sa3	=	0.01281;sb3	=	3.772;sc3	=	1.566;sa4	=	0.004448;sb4	=	6.216;sc4	=	-1.567;}
            if(params.temperature==5.0 && params.density==0.25) { sa1	=	0.3907;sb1	=	0.8563;sc1	=	1.574;sa2	=	0.04736;sb2	=	2.826;sc2	=	1.588;sa3	=	-0.006766;sb3	=	5.223;sc3	=	-1.596;sa4	=	-0.0261;sb4	=	4.332;sc4	=	1.555; }
            if(params.temperature==5.0 && params.density==0.375) { sa1	=	0.5497;sb1	=	0.8236;sc1	=	1.569;sa2	=	0.03292;sb2	=	1.549;sc2	=	-1.562;sa3	=	-0.01415;sb3	=	5.14;sc3	=	1.479;sa4	=	-0.01071;sb4	=	5.587;sc4	=	-1.624; }
            if(params.temperature==5.0 && params.density==0.5) { sa1	=	0.7866;sb1	=	0.9087;sc1	=	1.573;sa2	=	0.2075;sb2	=	1.814;sc2	=	-1.578;sa3	=	0.02799;sb3	=	3.635;sc3	=	1.606;sa4	=	0.009727;sb4	=	5.452;sc4	=	-1.621;sa5	=	0.006593;sb5	=	7.264;sc5	=	1.624;sa6	=	0.003764;sb6	=	9.083;sc6	=	-1.199;sa7	=	0.001517;sb7	=	18.17;sc7	=	-1.278;sa8	=	0.00151;sb8	=	10.9;sc8	=	0.9269; }
            if(params.temperature==5.0 && params.density==0.625) { sa1	=	1.036;sb1	=	0.9926;sc1	=	1.573;sa2	=	0.4123;sb2	=	1.964;sc2	=	-1.57;sa3	=	0.1168;sb3	=	3.786;sc3	=	1.562;sa4	=	0.4995;sb4	=	5.521;sc4	=	-1.574;sa5	=	0.4583;sb5	=	5.615;sc5	=	1.566;sa6	=	0.004006;sb6	=	9.054;sc6	=	-1.055;sa7	=	0.002954;sb7	=	16.57;sc7	=	1.23;sa8	=	0.003006;sb8	=	17.94;sc8	=	-1.403; }
            if(params.temperature==5.0 && params.density==0.75) { sa1	=	1.286;sb1	=	0.9526;sc1	=	1.573;sa2	=	0.5388;sb2	=	1.71;sc2	=	-1.569;sa3	=	0.6026;sb3	=	4.409;sc3	=	1.586;sa4	=	0.5634;sb4	=	4.537;sc4	=	-1.551;sa5	=	0.01431;sb5	=	7.245;sc5	=	1.727;sa6	=	0.01098;sb6	=	9.394;sc6	=	-1.428;sa7	=	0.00889;sb7	=	10.82;sc7	=	1.714;sa8	=	0.005029;sb8	=	16.32;sc8	=	1.553; }
            if(params.temperature==5.0 && params.density==0.875) { sa1	=	1.577;sb1	=	0.9674;sc1	=	1.572;sa2	=	0.7353;sb2	=	1.667;sc2	=	-1.582;sa3	=	0.1001;sb3	=	3.699;sc3	=	1.63;sa4	=	0.04359;sb4	=	5.56;sc4	=	-1.623;sa5	=	0.02722;sb5	=	7.299;sc5	=	1.517;sa6	=	0.01905;sb6	=	9.064;sc6	=	-1.756;sa7	=	0.009713;sb7	=	10.88;sc7	=	1.665;sa8	=	0.006004;sb8	=	16.28;sc8	=	1.466; }
            if(params.temperature==5.0 && params.density==1.0) { sa1	=	1.568;sb1	=	0.9074;sc1	=	1.568;sa2	=	0.6369;sb2	=	1.819;sc2	=	-1.575;sa3	=	0.1242;sb3	=	3.63;sc3	=	1.569;sa4	=	0.05096;sb4	=	5.453;sc4	=	-1.598;sa5	=	0.03115;sb5	=	7.265;sc5	=	1.507;sa6	=	0.01709;sb6	=	9.083;sc6	=	-1.713;sa7	=	0.01261;sb7	=	16.35;sc7	=	1.994;sa8	=	0.0115;sb8	=	10.9;sc8	=	1.644; }


            MappedSingletDensity sim = new MappedSingletDensity(params);
            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);

            int blockSize = 100;
            MeterProfileByVolume densityMeter = new MeterProfileByVolume(sim.space);
            densityMeter.setBox(sim.box());
            densityMeter.setProfileDim(2);
            MeterNMolecules meterNMolecules = new MeterNMolecules();
            densityMeter.setDataSource(meterNMolecules);
            AccumulatorAverageFixed acc = new AccumulatorAverageFixed(5 * blockSize);
            DataPumpListener pump = new DataPumpListener(densityMeter, acc, params.numAtoms);
            sim.getIntegrator().getEventManager().addListener(pump);

            MeterProfileForceSum densityMeterForce = new MeterProfileForceSum(sim.box(), sim.potentialMaster, params.temperature);
            densityMeterForce.setProfileDim(2);
            AccumulatorAverageFixed accForce = new AccumulatorAverageFixed(blockSize);
            DataPumpListener pumpForce = new DataPumpListener(densityMeterForce, accForce, 5 * params.numAtoms);
            sim.getIntegrator().getEventManager().addListener(pumpForce);

            FunctionDifferentiable f;
            double L = sim.box().getBoundary().getBoxSize().getX(2);
            switch (params.field) {
                case PARABOLIC:
                    f = new FunctionParabolic(L, params.temperature);
                    break;
                case SINE:
                    f = new FunctionSine(5, L);
                    break;
                case UNIFORM:
                    f = new FunctionUniform(L);
                    break;
                case LNPARABOLIC:
                    f = new FunctionLnparabolic(L);
                    break;
                case PHISINEPSINESUM:
                    f = new FunctionPhisinepsinesum(L,aa,bb,cc);
                    break;
                case PHIPARABOLICPSUMOFGAUSSIANS:
                    f = new FunctionPhiparabolicpsumofgaussians(L,a1,b1,c1,a2,b2,c2,a3,b3,c3,a4,b4,c4,a5,b5,c5);
                    break;
                case EXPMINUSZSQ:
                    f = new FunctionExpminuszsq(L,aa0,aa1,aa2,aa3,aa4,aa5,aa6,aa7,aa8,bb1,bb2,bb3,bb4,bb5,bb6,bb7,bb8,w);
                    break;
                case PHIARGPT5PARABOLICPSINESUM:
                    f = new FunctionPhiargpt5parabolicpsinesum(L,sa1,sb1,sc1,sa2,sb2,sc2,sa3,sb3,sc3,sa4,sb4,sc4,sa5,sb5,sc5,sa6,sb6,sc6,sa7,sb7,sc7,sa8,sb8,sc8);
                    break;
                default:
                    throw new RuntimeException("not yet");
            }
            MeterProfileMappedAvg densityMeterMappedAvg = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterMappedAvg.setProfileDim(2);
            AccumulatorAverageFixed accMappedAvg = new AccumulatorAverageFixed(blockSize);
            DataPumpListener pumpMappedAvg = new DataPumpListener(densityMeterMappedAvg, accMappedAvg, 5 * params.numAtoms);
            sim.getIntegrator().getEventManager().addListener(pumpMappedAvg);

            MeterProfileMappedAvg densityMeterP = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterP.getXDataSource().setNValues(104);
            densityMeterP.setProfileDim(2);
            densityMeterP.setBehavior(MeterProfileMappedAvg.Behavior.P);
            MeterProfileMappedAvg densityMeterZidot0 = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterZidot0.getXDataSource().setNValues(104);
            densityMeterZidot0.setProfileDim(2);
            densityMeterZidot0.setBehavior(MeterProfileMappedAvg.Behavior.ZIDOT);
            densityMeterZidot0.setZidotZ(0);
            MeterProfileMappedAvg densityMeterZidot1 = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterZidot1.getXDataSource().setNValues(104);
            densityMeterZidot1.setProfileDim(2);
            densityMeterZidot1.setBehavior(MeterProfileMappedAvg.Behavior.ZIDOT);
            densityMeterZidot1.setZidotZ(L / 8);
            MeterProfileMappedAvg densityMeterZidot2 = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterZidot2.getXDataSource().setNValues(104);
            densityMeterZidot2.setProfileDim(2);
            densityMeterZidot2.setBehavior(MeterProfileMappedAvg.Behavior.ZIDOT);
            densityMeterZidot2.setZidotZ(L / 4);

            MeterProfileMappedAvg densityMeterDZidot0 = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterDZidot0.getXDataSource().setNValues(104);
            densityMeterDZidot0.setProfileDim(2);
            densityMeterDZidot0.setBehavior(MeterProfileMappedAvg.Behavior.DZIDOT);
            densityMeterDZidot0.setZidotZ(0);
            MeterProfileMappedAvg densityMeterDZidot1 = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterDZidot1.getXDataSource().setNValues(104);
            densityMeterDZidot1.setProfileDim(2);
            densityMeterDZidot1.setBehavior(MeterProfileMappedAvg.Behavior.DZIDOT);
            densityMeterDZidot1.setZidotZ(L / 8);
            MeterProfileMappedAvg densityMeterDZidot2 = new MeterProfileMappedAvg(sim.box(), sim.potentialMaster, params.temperature, f);
            densityMeterDZidot2.getXDataSource().setNValues(104);
            densityMeterDZidot2.setProfileDim(2);
            densityMeterDZidot2.setBehavior(MeterProfileMappedAvg.Behavior.DZIDOT);
            densityMeterDZidot2.setZidotZ(L / 4);

            DisplayPlot densityPlot = new DisplayPlot();
            acc.addDataSink(densityPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.AVERAGE});
            accForce.addDataSink(densityPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.AVERAGE});
            accMappedAvg.addDataSink(densityPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.AVERAGE});

            new DataPump(densityMeterP, densityPlot.getDataSet().makeDataSink()).actionPerformed();

            densityPlot.setLabel("density");
            graphic.add(densityPlot);

            DisplayPlot errorPlot = new DisplayPlot();
            acc.addDataSink(errorPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.ERROR});
            accForce.addDataSink(errorPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.ERROR});
            accMappedAvg.addDataSink(errorPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.ERROR});
            errorPlot.setLabel("error");
            graphic.add(errorPlot);

            DisplayPlot corPlot = new DisplayPlot();
            acc.addDataSink(corPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.BLOCK_CORRELATION});
            accForce.addDataSink(corPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.BLOCK_CORRELATION});
            accMappedAvg.addDataSink(corPlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{acc.BLOCK_CORRELATION});
            corPlot.setLabel("correlation");
            graphic.add(corPlot);

            DisplayPlot zidotPlot = new DisplayPlot();
            new DataPump(densityMeterZidot0, zidotPlot.getDataSet().makeDataSink()).actionPerformed();
            new DataPump(densityMeterZidot1, zidotPlot.getDataSet().makeDataSink()).actionPerformed();
            new DataPump(densityMeterZidot2, zidotPlot.getDataSet().makeDataSink()).actionPerformed();
            zidotPlot.setLabel("zidot");
            graphic.add(zidotPlot);

            DisplayPlot dzidotPlot = new DisplayPlot();
            new DataPump(densityMeterDZidot0, dzidotPlot.getDataSet().makeDataSink()).actionPerformed();
            new DataPump(densityMeterDZidot1, dzidotPlot.getDataSet().makeDataSink()).actionPerformed();
            new DataPump(densityMeterDZidot2, dzidotPlot.getDataSet().makeDataSink()).actionPerformed();
            dzidotPlot.setLabel("-dzidot");
            graphic.add(dzidotPlot);

            List<DataPump> pumps = graphic.getController().getDataStreamPumps();
            pumps.add(pump);
            pumps.add(pumpForce);
            pumps.add(pumpMappedAvg);

            graphic.makeAndDisplayFrame();
        }
    }

    enum Field {
        SINE, UNIFORM, PARABOLIC, PHISINEPSINESUM, PHIPARABOLICPSUMOFGAUSSIANS, LNPARABOLIC, PHILNPARABOLICPFOURIERSUM,PHIARGPT5PARABOLICPSINESUM,EXPMINUSZSQ
    }

    public static class SimParams extends ParameterBase {
        public long steps = 100000000;
        public double density = 0.125;
        public int bins = 1000;
        public double temperature = 5;
        public int numAtoms = 500;
        public Field field = Field.UNIFORM;
    }
}
