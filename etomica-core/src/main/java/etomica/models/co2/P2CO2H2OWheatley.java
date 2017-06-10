/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.co2;

import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.potential.IPotentialAtomic;
import etomica.potential.IPotentialTorque;
import etomica.simulation.Simulation;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.IOrientation3D;
import etomica.space3d.OrientationFull3D;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.BohrRadius;
import etomica.units.Hartree;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.random.RandomMersenneTwister;
import etomica.util.random.RandomNumberGeneratorUnix;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Ab initio potential for CO2-H2O by Wheatley and Harvey
 *
 * http://dx.doi.org/10.1063/1.3574345
 * 
 * The CO2 molecule is expected first
 * 
 * @author Andrew Schultz
 */
public class P2CO2H2OWheatley implements IPotentialTorque {

    protected static final int[][] iparams = new int[][]{
        {1,4,0,0,1,0},{1,4,0,0,6,0},{1,4,0,0,12,0},{1,5,0,0,1,0},{1,5,0,0,6,0}, // 5
        {1,5,0,0,12,0},{1,6,0,0,1,0},{1,6,0,0,6,0},{1,6,0,0,12,0},{2,4,0,0,1,0}, // 10
        {2,4,0,0,6,0},{2,4,0,0,12,0},{2,5,0,0,1,0},{2,5,0,0,6,0},{2,5,0,0,12,0}, // 15
        {2,6,0,0,1,0},{2,6,0,0,6,0},{2,6,0,0,12,0},{3,4,0,0,1,0},{3,4,0,0,6,0}, // 20
        {3,4,0,0,12,0},{3,5,0,0,1,0},{3,5,0,0,6,0},{3,5,0,0,12,0},{3,6,0,0,1,0}, // 25
        {3,6,0,0,6,0},{3,6,0,0,12,0},{1,4,1,5,1,5},{1,4,1,5,2,4},{1,4,1,5,3,3}, // 30
        {1,4,1,5,4,2},{1,4,1,5,5,1},{1,4,1,5,4,8},{1,4,1,5,6,6},{1,4,1,5,8,4}, // 35
        {1,4,1,5,4,4},{1,4,1,6,1,5},{1,4,1,6,2,4},{1,4,1,6,3,3},{1,4,1,6,4,2}, // 40
        {1,4,1,6,5,1},{1,4,1,6,4,8},{1,4,1,6,6,6},{1,4,1,6,8,4},{1,4,1,6,4,4}, // 45
        {1,4,2,4,1,5},{1,4,2,4,2,4},{1,4,2,4,3,3},{1,4,2,4,4,2},{1,4,2,4,5,1}, // 50
        {1,4,2,4,4,8},{1,4,2,4,6,6},{1,4,2,4,8,4},{1,4,2,4,4,4},{1,4,2,5,1,5}, // 55
        {1,4,2,5,2,4},{1,4,2,5,3,3},{1,4,2,5,4,2},{1,4,2,5,5,1},{1,4,2,5,4,8}, // 60
        {1,4,2,5,6,6},{1,4,2,5,8,4},{1,4,2,5,4,4},{1,4,2,6,1,5},{1,4,2,6,2,4}, // 65
        {1,4,2,6,3,3},{1,4,2,6,4,2},{1,4,2,6,5,1},{1,4,2,6,4,8},{1,4,2,6,6,6}, // 70
        {1,4,2,6,8,4},{1,4,2,6,4,4},{1,4,3,4,1,5},{1,4,3,4,2,4},{1,4,3,4,3,3}, // 75
        {1,4,3,4,4,2},{1,4,3,4,5,1},{1,4,3,4,4,8},{1,4,3,4,6,6},{1,4,3,4,8,4}, // 80
        {1,4,3,4,4,4},{1,4,3,5,1,5},{1,4,3,5,2,4},{1,4,3,5,3,3},{1,4,3,5,4,2}, // 85
        {1,4,3,5,5,1},{1,4,3,5,4,8},{1,4,3,5,6,6},{1,4,3,5,8,4},{1,4,3,5,4,4}, // 90
        {1,4,3,6,1,5},{1,4,3,6,2,4},{1,4,3,6,3,3},{1,4,3,6,4,2},{1,4,3,6,5,1}, // 95
        {1,4,3,6,4,8},{1,4,3,6,6,6},{1,4,3,6,8,4},{1,4,3,6,4,4},{1,5,1,6,1,5}, // 100
        {1,5,1,6,2,4},{1,5,1,6,3,3},{1,5,1,6,4,2},{1,5,1,6,5,1},{1,5,1,6,4,8}, // 105
        {1,5,1,6,6,6},{1,5,1,6,8,4},{1,5,1,6,4,4},{1,5,2,4,1,5},{1,5,2,4,2,4}, // 110
        {1,5,2,4,3,3},{1,5,2,4,4,2},{1,5,2,4,5,1},{1,5,2,4,4,8},{1,5,2,4,6,6}, // 115
        {1,5,2,4,8,4},{1,5,2,4,4,4},{1,5,2,5,1,5},{1,5,2,5,2,4},{1,5,2,5,3,3}, // 120
        {1,5,2,5,4,2},{1,5,2,5,5,1},{1,5,2,5,4,8},{1,5,2,5,6,6},{1,5,2,5,8,4}, // 125
        {1,5,2,5,4,4},{1,5,2,6,1,5},{1,5,2,6,2,4},{1,5,2,6,3,3},{1,5,2,6,4,2}, // 130
        {1,5,2,6,5,1},{1,5,2,6,4,8},{1,5,2,6,6,6},{1,5,2,6,8,4},{1,5,2,6,4,4}, // 135
        {1,5,3,4,1,5},{1,5,3,4,2,4},{1,5,3,4,3,3},{1,5,3,4,4,2},{1,5,3,4,5,1}, // 140
        {1,5,3,4,4,8},{1,5,3,4,6,6},{1,5,3,4,8,4},{1,5,3,4,4,4},{1,5,3,5,1,5}, // 145
        {1,5,3,5,2,4},{1,5,3,5,3,3},{1,5,3,5,4,2},{1,5,3,5,5,1},{1,5,3,5,4,8}, // 150
        {1,5,3,5,6,6},{1,5,3,5,8,4},{1,5,3,5,4,4},{1,5,3,6,1,5},{1,5,3,6,2,4}, // 155
        {1,5,3,6,3,3},{1,5,3,6,4,2},{1,5,3,6,5,1},{1,5,3,6,4,8},{1,5,3,6,6,6}, // 160
        {1,5,3,6,8,4},{1,5,3,6,4,4},{1,6,2,4,1,5},{1,6,2,4,2,4},{1,6,2,4,3,3}, // 165
        {1,6,2,4,4,2},{1,6,2,4,5,1},{1,6,2,4,4,8},{1,6,2,4,6,6},{1,6,2,4,8,4}, // 170
        {1,6,2,4,4,4},{1,6,2,5,1,5},{1,6,2,5,2,4},{1,6,2,5,3,3},{1,6,2,5,4,2}, // 175
        {1,6,2,5,5,1},{1,6,2,5,4,8},{1,6,2,5,6,6},{1,6,2,5,8,4},{1,6,2,5,4,4}, // 180
        {1,6,2,6,1,5},{1,6,2,6,2,4},{1,6,2,6,3,3},{1,6,2,6,4,2},{1,6,2,6,5,1}, // 185
        {1,6,2,6,4,8},{1,6,2,6,6,6},{1,6,2,6,8,4},{1,6,2,6,4,4},{1,6,3,4,1,5}, // 190
        {1,6,3,4,2,4},{1,6,3,4,3,3},{1,6,3,4,4,2},{1,6,3,4,5,1},{1,6,3,4,4,8}, // 195
        {1,6,3,4,6,6},{1,6,3,4,8,4},{1,6,3,4,4,4},{1,6,3,5,1,5},{1,6,3,5,2,4}, // 200
        {1,6,3,5,3,3},{1,6,3,5,4,2},{1,6,3,5,5,1},{1,6,3,5,4,8},{1,6,3,5,6,6}, // 205
        {1,6,3,5,8,4},{1,6,3,5,4,4},{1,6,3,6,1,5},{1,6,3,6,2,4},{1,6,3,6,3,3}, // 210
        {1,6,3,6,4,2},{1,6,3,6,5,1},{1,6,3,6,4,8},{1,6,3,6,6,6},{1,6,3,6,8,4}, // 215
        {1,6,3,6,4,4},{2,4,2,5,1,5},{2,4,2,5,2,4},{2,4,2,5,3,3},{2,4,2,5,4,2}, // 220
        {2,4,2,5,5,1},{2,4,2,5,4,8},{2,4,2,5,6,6},{2,4,2,5,8,4},{2,4,2,5,4,4}, // 225
        {2,4,2,6,1,5},{2,4,2,6,2,4},{2,4,2,6,3,3},{2,4,2,6,4,2},{2,4,2,6,5,1}, // 230
        {2,4,2,6,4,8},{2,4,2,6,6,6},{2,4,2,6,8,4},{2,4,2,6,4,4},{2,4,3,4,1,5}, // 235
        {2,4,3,4,2,4},{2,4,3,4,3,3},{2,4,3,4,4,2},{2,4,3,4,5,1},{2,4,3,4,4,8}, // 240
        {2,4,3,4,6,6},{2,4,3,4,8,4},{2,4,3,4,4,4},{2,4,3,5,1,5},{2,4,3,5,2,4}, // 245
        {2,4,3,5,3,3},{2,4,3,5,4,2},{2,4,3,5,5,1},{2,4,3,5,4,8},{2,4,3,5,6,6}, // 250
        {2,4,3,5,8,4},{2,4,3,5,4,4},{2,4,3,6,1,5},{2,4,3,6,2,4},{2,4,3,6,3,3}, // 255
        {2,4,3,6,4,2},{2,4,3,6,5,1},{2,4,3,6,4,8},{2,4,3,6,6,6},{2,4,3,6,8,4}, // 260
        {2,4,3,6,4,4},{2,5,2,6,1,5},{2,5,2,6,2,4},{2,5,2,6,3,3},{2,5,2,6,4,2}, // 265
        {2,5,2,6,5,1},{2,5,2,6,4,8},{2,5,2,6,6,6},{2,5,2,6,8,4},{2,5,2,6,4,4}, // 270
        {2,5,3,4,1,5},{2,5,3,4,2,4},{2,5,3,4,3,3},{2,5,3,4,4,2},{2,5,3,4,5,1}, // 275
        {2,5,3,4,4,8},{2,5,3,4,6,6},{2,5,3,4,8,4},{2,5,3,4,4,4},{2,5,3,5,1,5}, // 280
        {2,5,3,5,2,4},{2,5,3,5,3,3},{2,5,3,5,4,2},{2,5,3,5,5,1},{2,5,3,5,4,8}, // 285
        {2,5,3,5,6,6},{2,5,3,5,8,4},{2,5,3,5,4,4},{2,5,3,6,1,5},{2,5,3,6,2,4}, // 290
        {2,5,3,6,3,3},{2,5,3,6,4,2},{2,5,3,6,5,1},{2,5,3,6,4,8},{2,5,3,6,6,6}, // 295
        {2,5,3,6,8,4},{2,5,3,6,4,4},{2,6,3,4,1,5},{2,6,3,4,2,4},{2,6,3,4,3,3}, // 300
        {2,6,3,4,4,2},{2,6,3,4,5,1},{2,6,3,4,4,8},{2,6,3,4,6,6},{2,6,3,4,8,4}, // 305
        {2,6,3,4,4,4},{2,6,3,5,1,5},{2,6,3,5,2,4},{2,6,3,5,3,3},{2,6,3,5,4,2}, // 310
        {2,6,3,5,5,1},{2,6,3,5,4,8},{2,6,3,5,6,6},{2,6,3,5,8,4},{2,6,3,5,4,4}, // 315
        {2,6,3,6,1,5},{2,6,3,6,2,4},{2,6,3,6,3,3},{2,6,3,6,4,2},{2,6,3,6,5,1}, // 320
        {2,6,3,6,4,8},{2,6,3,6,6,6},{2,6,3,6,8,4},{2,6,3,6,4,4},{3,4,3,5,1,5}, // 325
        {3,4,3,5,2,4},{3,4,3,5,3,3},{3,4,3,5,4,2},{3,4,3,5,5,1},{3,4,3,5,4,8}, // 330
        {3,4,3,5,6,6},{3,4,3,5,8,4},{3,4,3,5,4,4},{3,4,3,6,1,5},{3,4,3,6,2,4}, // 335
        {3,4,3,6,3,3},{3,4,3,6,4,2},{3,4,3,6,5,1},{3,4,3,6,4,8},{3,4,3,6,6,6}, // 340
        {3,4,3,6,8,4},{3,4,3,6,4,4},{3,5,3,6,1,5},{3,5,3,6,2,4},{3,5,3,6,3,3}, // 345
        {3,5,3,6,4,2},{3,5,3,6,5,1},{3,5,3,6,4,8},{3,5,3,6,6,6},{3,5,3,6,8,4}, // 350
        {3,5,3,6,4,4}
    };
    
    // for neutrality, fparams[0] = -0.410682122652
    protected static final double[] fparams = new double[]{
        -0.410682122651,329.403386615392,-3662553.726999491453,0.205341061326,62.365152431380, // 5
        7582.858590874342,0.205341061326,62.365152431376,7582.858590873629,0.205341061326, // 10
        -2536.074260430988,669159.950024581165,-0.102670530663,-287.307644954777,9264.030266577680, // 15
        -0.102670530663,-287.307644954776,9264.030266577214,0.205341061326,-2536.074260430983, // 20
        669159.950024574529,-0.102670530663,-287.307644954776,9264.030266579422,-0.102670530663, // 25
        -287.307644954776,9264.030266557407,-486.178549321415,-40.209398150392,618.218855154491, // 30
        539.481682660762,-21.578841432931,-66833.954647235179,601923.312853231910,-919444.207744437037, // 35
        -2282.248036649052,-486.178549347973,-40.209398148894,618.218855141273,539.481682627901, // 40
        -21.578841447981,-66833.954647217222,601923.312853238662,-919444.207744430751,-2282.248036647884, // 45
        -252.721059889493,257.227861069700,100.795810524327,4.904319026998,137.933882614129, // 50
        1303997.287576365052,-9531192.092034170404,8256957.595972028561,21045.274819385319,866.004895345221, // 55
        229.170163896768,-471.311819717143,-448.459994793814,-423.382264173876,-320614.124487735797, // 60
        41750.290727830041,811319.998104245053,-5.089190663452,866.004895340688,229.170163876742, // 65
        -471.311819812918,-448.459994843892,-423.382264221343,-320614.124487735215,41750.290727794527, // 70
        811319.998104202677,-5.089190662862,-252.721059839961,257.227861112610,100.795810520979, // 75
        4.904318998592,137.933882592005,1303997.287576298695,-9531192.092034230009,8256957.595971997827, // 80
        21045.274819386177,866.004895348882,229.170163919666,-471.311819639149,-448.459994709298, // 85
        -423.382264197864,-320614.124487736437,41750.290727830215,811319.998104272061,-5.089190661581, // 90
        866.004895338573,229.170163898811,-471.311819707802,-448.459994701939,-423.382264173612, // 95
        -320614.124487748777,41750.290727834814,811319.998104226426,-5.089190663279,-198.860530662686, // 100
        -161.300976304516,-90.578689058104,-161.300976262708,-198.860530660960,29082.932127982462, // 105
        -283766.747027040808,29082.932127970453,6759.676869134118,1891.911166076532,170.699802044474, // 110
        -893.998651829788,-42.311447006776,451.634271166922,-1114184.426450481871,1184071.967609238811, // 115
        -306627.094639008690,-3592.469970187148,-99.755683683332,95.896168275822,156.428973170855, // 120
        68.359150680550,8.553182192121,-13780.947467422735,13298.882434405145,-40550.223331271460, // 125
        6089.761613821893,-123.554543385642,178.789596416152,-4.966429182017,10.920516857571, // 130
        105.682989904764,23502.916714671901,-17579.494768423119,9842.064321376052,-1205.220460076951, // 135
        1891.911166003081,170.699802031625,-893.998651813597,-42.311447021553,451.634271158530, // 140
        -1114184.426450459519,1184071.967609260464,-306627.094639009330,-3592.469970186777,-99.755683640776, // 145
        95.896168252015,156.428973169444,68.359150682547,8.553182170931,-13780.947467422959, // 150
        13298.882434396301,-40550.223331271576,6089.761613822790,-123.554543391370,178.789596404035, // 155
        -4.966429172690,10.920516846645,105.682989919725,23502.916714697800,-17579.494768412282, // 160
        9842.064321369966,-1205.220460076330,1891.911166078827,170.699802102709,-893.998651841541, // 165
        -42.311446977059,451.634271223439,-1114184.426450496772,1184071.967609263957,-306627.094638999843, // 170
        -3592.469970187306,-123.554543408066,178.789596403869,-4.966429223922,10.920516831703, // 175
        105.682989870520,23502.916714667903,-17579.494768428114,9842.064321378561,-1205.220460076847, // 180
        -99.755683663055,95.896168287128,156.428973177453,68.359150699453,8.553182181344, // 185
        -13780.947467425345,13298.882434394194,-40550.223331268622,6089.761613822568,1891.911166099821, // 190
        170.699802058577,-893.998651842820,-42.311447012376,451.634271187350,-1114184.426450514235, // 195
        1184071.967609287472,-306627.094639016665,-3592.469970187275,-123.554543363411,178.789596431993, // 200
        -4.966429146548,10.920516869876,105.682989946030,23502.916714663577,-17579.494768417906, // 205
        9842.064321387326,-1205.220460075945,-99.755683673250,95.896168275532,156.428973184732, // 210
        68.359150698674,8.553182173811,-13780.947467403739,13298.882434379908,-40550.223331265654, // 215
        6089.761613822919,997.125392259654,34.278607591125,-1659.697562007138,-465.811045939354, // 220
        2557.589133156398,-341146.783994329046,1014895.417913309531,-1033569.314793183818,1082.835693003179, // 225
        997.125392252379,34.278607576690,-1659.697562043677,-465.811045977848,2557.589133167420, // 230
        -341146.783994335448,1014895.417913318845,-1033569.314793183585,1082.835693002077,-261.764725492745, // 235
        654.136758553088,899.376868277612,654.136758543920,-261.764725553206,3138181.608432486653, // 240
        -3304854.139011175837,3138181.608432290610,-60845.039450381722,-1054.203464840779,-586.590549894482, // 245
        2094.096147115303,-425.984272418326,-2221.016134285837,1177256.517866900656,-4133498.991013298742, // 250
        2829292.445748528931,9507.515153724649,-1054.203464846250,-586.590549924458,2094.096147102502, // 255
        -425.984272484660,-2221.016134252017,1177256.517866872251,-4133498.991013232619,2829292.445748589933, // 260
        9507.515153724049,37.460473480616,50.848606046284,-87.940069220654,50.848606040769, // 265
        37.460473472647,-23351.337732167813,85829.985222433563,-23351.337732171723,-832.278410381126, // 270
        -2221.016134274729,-425.984272349422,2094.096147089235,-586.590549870994,-1054.203464876866, // 275
        2829292.445748496801,-4133498.991013328079,1177256.517866898561,9507.515153723660,-195.582073871740, // 280
        -45.765618540328,7.424357575607,-45.765618543371,-195.582073879490,381065.938593674393, // 285
        -3435123.759196050465,381065.938593679806,-27576.834771619055,-22.371371159362,240.419544182708, // 290
        -498.557468053889,240.419544166827,-22.371371197146,-35285.396494612971,65709.212732101834, // 295
        -35285.396494621120,818.759432015389,-2221.016134231167,-425.984272297278,2094.096147149725, // 300
        -586.590549850817,-1054.203464863569,2829292.445748568047,-4133498.991013324820,1177256.517866891110, // 305
        9507.515153725724,-22.371371120028,240.419544197652,-498.557468063896,240.419544178865, // 310
        -22.371371187494,-35285.396494620851,65709.212732116270,-35285.396494620756,818.759432014876, // 315
        -195.582073865387,-45.765618546576,7.424357609330,-45.765618571478,-195.582073864197, // 320
        381065.938593644067,-3435123.759196044877,381065.938593676721,-27576.834771620062,997.125392287735, // 325
        34.278607602479,-1659.697562103051,-465.811045938373,2557.589133247131,-341146.783994365716, // 330
        1014895.417913247016,-1033569.314793192898,1082.835693003887,997.125392258501,34.278607578982, // 335
        -1659.697562005123,-465.811045940191,2557.589133109340,-341146.783994313271,1014895.417913287645, // 340
        -1033569.314793178812,1082.835693003840,37.460473479493,50.848606052746,-87.940069223158, // 345
        50.848606041344,37.460473525859,-23351.337732129618,85829.985222404430,-23351.337732163920, // 350
        -832.278410380238};
    static {
        for (int i=0; i<fparams.length; i++) {
            fparams[i] = Hartree.UNIT.toSim(fparams[i]);
        }
    }
    
    protected static final double sitesCO2L = 2.2114;
    protected static final double sitesOH = 1.12169, sitesHH = 1.45365;
    protected final Vector[][] sitePos;
    protected final Space space;
    protected final Tensor rot0, rot1;
    protected final Vector or01, or02, orH2O3;
    protected final Vector[] allOr0, allOr1;
    protected final Vector[][] gradientAndTorque;
    protected final Vector drij, torque;
    protected final double mass1;
    
    public P2CO2H2OWheatley(Space space) {
        this.space = space;
        sitePos = new Vector[2][3];
        for (int i=0; i<2; i++) {
            for (int j=0; j<3; j++) {
                sitePos[i][j] = space.makeVector();
            }
        }
        or01 = space.makeVector();
        or02 = space.makeVector();
        orH2O3 = space.makeVector();
        allOr0 = new Vector[]{null, or01, or02};
        allOr1 = new Vector[]{null, null, orH2O3};
        rot0 = space.makeTensor();
        rot1 = space.makeTensor();
        gradientAndTorque = new Vector[2][2];
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                gradientAndTorque[i][j] = space.makeVector();
            }
        }
        drij = space.makeVector();
        torque = space.makeVector();
        mass1 = Oxygen.INSTANCE.getMass() + 2*Hydrogen.INSTANCE.getMass();
    }

    protected void getPerp(Vector or, Vector perp1, Vector perp2) {
        int max = 0;
        if (Math.abs(or.getX(1)) > Math.abs(or.getX(0))) max=1;
        if (Math.abs(or.getX(2)) > Math.abs(or.getX(max))) max=2;
        int min = 0;
        if (Math.abs(or.getX(1)) < Math.abs(or.getX(0))) min=1;
        if (Math.abs(or.getX(2)) < Math.abs(or.getX(min))) min=2;
        perp1.E(0);
        perp1.setX(min, or.getX(max));
        perp1.setX(max, -or.getX(min));
        perp1.normalize();
        perp1.PEa1Tv1(-perp1.dot(or), or);
        perp1.normalize();
        perp2.E(or);
        perp2.XE(perp1);
    }
    
    public double energy(IAtomList atoms) {
        IAtomOriented atom0 = (IAtomOriented)atoms.getAtom(0);
        IAtomOriented atom1 = (IAtomOriented)atoms.getAtom(1);
        Vector cm0 = atom0.getPosition();
        Vector cm1 = atom1.getPosition();
//        System.out.println("cm1 "+cm1);
        double bohrConv = BohrRadius.UNIT.fromSim(1);
//        System.out.println(bohrConv);
        for (int i=0; i<3; i++) {
            sitePos[0][i].Ea1Tv1(bohrConv, cm0);
            sitePos[1][i].Ea1Tv1(bohrConv, cm1);
        }
        
        // everything in BohrRadii after this
        
        //CO2
        Vector orCO2 = atom0.getOrientation().getDirection();
        sitePos[0][1].PEa1Tv1(+sitesCO2L, orCO2);
        sitePos[0][2].PEa1Tv1(-sitesCO2L, orCO2);
        Vector orH2O1 = atom1.getOrientation().getDirection();
        Vector orH2O2 = ((OrientationFull3D)atom1.getOrientation()).getSecondaryDirection();
        double cmx = 2*Hydrogen.INSTANCE.getMass()*sitesOH/mass1;
        sitePos[1][0].PEa1Tv1(-cmx, orH2O1);
        sitePos[1][1].PEa1Tv1(sitesOH-cmx, orH2O1);
        sitePos[1][2].PEa1Tv1(sitesOH-cmx, orH2O1);
        sitePos[1][1].PEa1Tv1(+sitesHH, orH2O2);
        sitePos[1][2].PEa1Tv1(-sitesHH, orH2O2);
        CO = Math.sqrt(sitePos[0][0].Mv1Squared(sitePos[1][0]));
        if (CO < coreCO && !debugging) return Double.POSITIVE_INFINITY;
        boolean checkme = false;
        if (CO < checkCoreCO) checkme = true;
//        System.out.println("CO2 C "+sitePos[0][0]);
//        System.out.println("CO2 O1 "+sitePos[0][1]);
//        System.out.println("CO2 O2 "+sitePos[0][2]);
//        System.out.println("H2O O "+sitePos[1][0]);
//        System.out.println("H2O H1 "+sitePos[1][1]);
//        System.out.println("H2O H2 "+sitePos[1][2]);
//        System.exit(1);
        double energy = 0;
//        System.out.println("*****");
        minR = Double.POSITIVE_INFINITY;
        for (int i=0; i<iparams.length; i++) {
            int ia1 = iparams[i][0]-1;
            int ib1 = iparams[i][1]-3-1;
            int ia2 = iparams[i][2]-1;
            int ib2 = iparams[i][3]-3-1;
            int ir1 = iparams[i][4];
            int ir2 = iparams[i][5];
            double c1 = fparams[i];
//            System.out.println("ia1 "+ia1+" "+sitePos[0][ia1]);
            double r21 = sitePos[1][ib1].Mv1Squared(sitePos[0][ia1]);
            double rr1 = Math.sqrt(r21);
            if (rr1 < minR) minR = rr1;
            if (rr1 < core && !debugging) {
                return Double.POSITIVE_INFINITY;
            }
            if (rr1 < checkCore) {
                checkme = true;
            }
//            System.out.println(rr1);
            double rval1 = 1.0/Math.pow(rr1, ir1);
            if (ia2>-1) {
//                System.out.println(ia1+" "+ib1+" "+ia2+" "+ib2+" "+ir1+" "+rr1+" "+drij);
                double r22 = sitePos[1][ib2].Mv1Squared(sitePos[0][ia2]);
                double rr2 = Math.sqrt(r22);
                if (rr2 < minR) minR = rr2;
                double rval2 = 1.0/Math.pow(rr2, ir2);
                energy += rval1*rval2*c1;
                if (rr2 < core && !debugging) {
                    return Double.POSITIVE_INFINITY;
                }
                if (rr2 < checkCore) {
                    checkme = true;
                }
            }
            else {
                drij.Ev1Mv2(sitePos[1][ib1], sitePos[0][ia1]);
//                System.out.println(sitePos[0][ia1]);
//                System.out.println(sitePos[1][ib1]);
//                System.out.println(ia1+" "+ib1+" "+ir1+" "+rr1+" "+drij);
                energy += rval1*c1;
            }
//            System.out.println((i+1)+" "+Hartree.UNIT.fromSim(energy));
        }
        if (false && checkme && !debugging && energy < 1000) {
            Vector mine = space.makeVector();
            debugging = true;
            mine.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
            System.out.println("** close attraction in CO2-H2O "+energy +" "+minR+" "+Math.sqrt(mine.squared())*bohrConv);
            for (int i=80; i<150; i++) {
                atom1.getPosition().E(atom0.getPosition());
                atom1.getPosition().PEa1Tv1(i*0.01, mine);
                double uu = energy(atoms);
                System.out.print(String.format("%6.3f  %6.3f  %6.3f  %g\n", minR, CO, 0.01*i*Math.sqrt(mine.squared())*bohrConv, uu));
            }
            System.exit(1);
        }
        else if (!debugging && energy < -10000) System.out.println("-- "+energy);
//        System.exit(1);
        return energy;
    }
    protected double minR, CO;
    protected static final double core = 3.1, coreCO = 4;
    protected static final double checkCore = 3.2, checkCoreCO = 4.2;
    protected boolean debugging = false;
    
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public void setBox(Box box) {
        
    }

    public int nBody() {
        return 2;
    }

    public Vector[][] gradientAndTorque(IAtomList atoms) {
        IAtomOriented atom0 = (IAtomOriented)atoms.getAtom(0);
        IAtomOriented atom1 = (IAtomOriented)atoms.getAtom(1);
        Vector cm0 = atom0.getPosition();
        Vector cm1 = atom1.getPosition();
        for (int i=0; i<3; i++) {
            sitePos[0][i].E(cm0);
            sitePos[1][i].E(cm1);
        }
        //CO2
        Vector orCO2 = atom0.getOrientation().getDirection();
        sitePos[0][1].PEa1Tv1(+sitesCO2L, orCO2);
        sitePos[0][2].PEa1Tv1(-sitesCO2L, orCO2);
        Vector orH2O1 = atom1.getOrientation().getDirection();
        Vector orH2O2 = ((OrientationFull3D)atom1.getOrientation()).getSecondaryDirection();
        double cmx = 2*Hydrogen.INSTANCE.getMass()*sitesOH/mass1;
        sitePos[1][0].PEa1Tv1(-cmx, orH2O1);
        sitePos[1][1].PEa1Tv1(+sitesOH-cmx, orH2O1);
        sitePos[1][2].PEa1Tv1(+sitesOH-cmx, orH2O1);
        sitePos[1][1].PEa1Tv1(+sitesHH, orH2O2);
        sitePos[1][2].PEa1Tv1(-sitesHH, orH2O2);
        
        getPerp(orCO2, or01, or02);
        allOr0[0] = orCO2;
        rot0.E(allOr0);
        rot0.invert();

        orH2O3.E(orH2O1);
        orH2O3.XE(orH2O2);
        allOr1[0] = orH2O1;
        allOr1[1] = orH2O2;
        rot1.E(allOr1);
        rot1.invert();
        gradientAndTorque[0][0].E(0);
        gradientAndTorque[0][1].E(0);
        gradientAndTorque[1][0].E(0);
        gradientAndTorque[1][1].E(0);
        
        for (int i=0; i<iparams.length; i++) {
            int ia1 = iparams[i][0]-1;
            int ib1 = iparams[i][1]-3-1;
            int ia2 = iparams[i][2]-1;
            int ib2 = iparams[i][3]-3-1;
            int ir1 = iparams[i][4];
            int ir2 = iparams[i][5];
            double c1 = fparams[i];
            
            drij.Ev1Mv2(sitePos[1][ib1], sitePos[0][ia1]);
            double r21 = sitePos[1][ib1].Mv1Squared(sitePos[0][ia1]);
            double rr1 = Math.sqrt(r21);
            double rval1 = 1.0/Math.pow(rr1, ir1);
            if (ia2>-1) {
                double r22 = sitePos[1][ib2].Mv1Squared(sitePos[0][ia2]);
                double rr2 = Math.sqrt(r22);
                double rval2 = 1.0/Math.pow(rr2, ir2);
            }
            else {
                double rdu = -ir1*c1*rval1;
//                System.out.println(ia1+" "+ib1+" "+ir1+" "+rr1+" "+rdu/rr1+" "+drij);

                drij.TE(rdu/r21);
                gradientAndTorque[0][1].PE(drij);
//                System.out.println("grad 1 "+gradientAndTorque[0][1]);
                gradientAndTorque[0][0].ME(drij);
                torque.Ev1Mv2(sitePos[1][ib1], cm1);
                torque.XE(drij);
                gradientAndTorque[1][1].ME(torque);
                torque.Ev1Mv2(sitePos[0][ia1], cm0);
                torque.XE(drij);
                gradientAndTorque[1][0].PE(torque);
            }
        }
        return gradientAndTorque;
    }
//    protected final int istart=3, istop=4;

    /**
     * Creates a Feynman-Hibbs implementation of a semi-classical potential based on this model for the given temperature.
     *
     */
    public P2CO2H2OSC makeSemiclassical(double temperature) {
        return new P2CO2H2OSC(temperature);
    }

    /**
     * Creates a Takahashi-Imada implementation of a semi-classical potential based on this model for the given temperature.
     *
     */
    public P2CO2H2OSCTI makeSemiclassicalTI(double temperature) {
        return new P2CO2H2OSCTI(temperature);
    }

    public class P2CO2H2OSC implements IPotentialAtomic {

        protected final Tensor tt0Tensor, rr0Tensor, rr1Tensor;
        protected final Tensor ijTensor, rTensor0, rTensor1, identity;
        protected final Tensor ijRTensor;
        protected final double mass0;
        protected final double moment0;
        protected final Vector moment1;
        protected final double fac;
        protected final Vector d2, r2a, r2ab, r2b, ga, gb;
        protected final Vector drijRot;
        protected final Vector bdrij;

        public P2CO2H2OSC(double temperature) {
            ijTensor = space.makeTensor();
            identity = space.makeTensor();
            tt0Tensor = space.makeTensor();
            rr0Tensor = space.makeTensor();
            rr1Tensor = space.makeTensor();
            rTensor0 = space.makeTensor();
            rTensor1 = space.makeTensor();
            ijRTensor = space.makeTensor();
            identity.E(new double[][]{{1,0,0},{0,1,0},{0,0,1}});
            drijRot = space.makeVector();
            mass0 = Carbon.INSTANCE.getMass() + 2*Oxygen.INSTANCE.getMass();
            moment0 = 2*Oxygen.INSTANCE.getMass()*sitesCO2L*sitesCO2L;

            moment1 = space.makeVector();
            double cm1x = 2*Hydrogen.INSTANCE.getMass()*sitesOH/mass1;
            moment1.setX(0, 2*Hydrogen.INSTANCE.getMass()*(sitesHH*sitesHH));
            moment1.setX(1, 2*Hydrogen.INSTANCE.getMass()*(sitesOH-cm1x)*(sitesOH-cm1x) + Oxygen.INSTANCE.getMass()*cm1x*cm1x);
            double cmH2 = (sitesOH-cm1x)*(sitesOH-cm1x) + (sitesHH*sitesHH);
            moment1.setX(2, 2*Hydrogen.INSTANCE.getMass()*cmH2 + Oxygen.INSTANCE.getMass()*cm1x*cm1x);
            d2 = space.makeVector();
            
            double hbar = Constants.PLANCK_H/(2*Math.PI);
            fac = hbar*hbar/(24/2)/temperature;
            bdrij = space.makeVector();
            r2a = space.makeVector();
            r2ab = space.makeVector();
            r2b = space.makeVector();
            ga = space.makeVector();
            gb = space.makeVector();
        }
        
        public double getRange() {
            return Double.POSITIVE_INFINITY;
        }

        public void setBox(Box box) {
            
        }

        public int nBody() {
            return 2;
        }
        public double[][] d2tot = new double[2][6];
        public double energy(IAtomList atoms) {
            IAtomOriented atom0 = (IAtomOriented)atoms.getAtom(0);
            IAtomOriented atom1 = (IAtomOriented)atoms.getAtom(1);
            Vector cm0 = atom0.getPosition();
            Vector cm1 = atom1.getPosition();
            double bohrConv = BohrRadius.UNIT.fromSim(1);
            for (int i=0; i<3; i++) {
                sitePos[0][i].Ea1Tv1(bohrConv, cm0);
                sitePos[1][i].Ea1Tv1(bohrConv, cm1);
            }

            //CO2
            Vector orCO2 = atom0.getOrientation().getDirection();
            sitePos[0][1].PEa1Tv1(+sitesCO2L, orCO2);
            sitePos[0][2].PEa1Tv1(-sitesCO2L, orCO2);
            Vector orH2O1 = atom1.getOrientation().getDirection();
            Vector orH2O2 = ((OrientationFull3D)atom1.getOrientation()).getSecondaryDirection();
            double cmx = 2*Hydrogen.INSTANCE.getMass()*sitesOH/mass1;
            sitePos[1][0].PEa1Tv1(-cmx, orH2O1);
            sitePos[1][1].PEa1Tv1(+sitesOH-cmx, orH2O1);
            sitePos[1][2].PEa1Tv1(+sitesOH-cmx, orH2O1);
            sitePos[1][1].PEa1Tv1(+sitesHH, orH2O2);
            sitePos[1][2].PEa1Tv1(-sitesHH, orH2O2);
            
            getPerp(orCO2, or01, or02);
            allOr0[0] = orCO2;
            rot0.E(allOr0);
            rot0.invert();

            orH2O3.E(orH2O1);
            orH2O3.XE(orH2O2);
            allOr1[0] = orH2O1;
            allOr1[1] = orH2O2;
            rot1.E(allOr1);
            rot1.invert();
            
            double d2tsum = 0, d2rsum = 0;
            double d2x = 0, d2y = 0, d2z = 0;
            double d2tx = 0, d2ty = 0, d2tz = 0;
            for (int i=0; i<2; i++) {
                for (int j=0; j<6; j++) {
                    d2tot[i][j] = 0;
                }
            }
            
            CO = Math.sqrt(sitePos[0][0].Mv1Squared(sitePos[1][0]));
            if (CO < 4 && !debugging) return Double.POSITIVE_INFINITY;
            double energy = 0;
            for (int i=0; i<iparams.length; i++) {
                int ia1 = iparams[i][0]-1;
                int ib1 = iparams[i][1]-3-1;
                int ia2 = iparams[i][2]-1;
                int ib2 = iparams[i][3]-3-1;
                int ir1 = iparams[i][4];
                int ir2 = iparams[i][5];
                double c1 = fparams[i];
                
                drij.Ev1Mv2(sitePos[1][ib1], sitePos[0][ia1]);
                double r21 = drij.squared();
                double rr1 = Math.sqrt(r21);
                double rval1 = 1.0/Math.pow(rr1, ir1);
                if (rr1 < core) {
                    return Double.POSITIVE_INFINITY;
                }
                if (ia2>-1) {
//                    System.out.println(ia1+" "+ib1+" "+ia2+" "+ib2+" "+ir1+" "+ir2);
                    bdrij.Ev1Mv2(sitePos[1][ib2], sitePos[0][ia2]);
                    double r22 = bdrij.squared();
                    double rr2 = Math.sqrt(r22);
                    if (rr2 < core) {
                        return Double.POSITIVE_INFINITY;
                    }
                    double rval2 = 1.0/Math.pow(rr2, ir2);
                    energy += rval1*rval2*c1;
                    
                    double adu = -c1*ir1*rval1*rval2/r21;
                    double bdu = -c1*ir2*rval1*rval2/r22;
                    double adu2 = -adu*(ir1+2)/r21;
                    double abdu2 = -adu*ir2/r22;
                    double bdu2 = -bdu*(ir2+2)/r22;
                    
                    drijRot.E(drij);
                    rot0.transform(drijRot);
                    ga.Ea1Tv1(adu, drijRot);

                    r2a.E(drijRot);
                    r2a.TE(drijRot);
                    r2ab.E(drijRot);
                    drijRot.E(bdrij);
                    rot0.transform(drijRot);
                    gb.Ea1Tv1(bdu, drijRot);
                    r2ab.TE(drijRot);
                    r2b.E(drijRot);
                    r2b.TE(drijRot);

                    double dxa = 0;
                    if (ia1>0) {
                        dxa = (2*ia1-3)*sitesCO2L;
                    }
                    double dxb = 0;
                    if (ia2>0) {
                        dxb = (2*ia2-3)*sitesCO2L;
                    }
                    
                    double d2a = adu + r2a.getX(0)*adu2;
                    double d2ab = r2ab.getX(0)*abdu2;
                    double d2b = bdu + r2b.getX(0)*bdu2;

                    double d2x0 = d2a + 2*d2ab + d2b;

                    d2a = adu + r2a.getX(1)*adu2;
                    d2ab = r2ab.getX(1)*abdu2;
                    d2b = bdu + r2b.getX(1)*bdu2;
                    
                    double d2y0 = d2a + 2*d2ab + d2b;
                    double d2thetaz0 = d2a*dxa*dxa + 2*d2ab*dxa*dxb + d2b*dxb*dxb - ga.getX(0)*dxa - gb.getX(0)*dxb;

                    d2a = adu + r2a.getX(2)*adu2;
                    d2ab = r2ab.getX(2)*abdu2;
                    d2b = bdu + r2b.getX(2)*bdu2;
                    
                    double d2z0 = d2a + 2*d2ab + d2b;
                    double d2thetay0 = d2a*dxa*dxa + 2*d2ab*dxa*dxb + d2b*dxb*dxb - ga.getX(0)*dxa - gb.getX(0)*dxb;
                    d2tot[0][0] += d2x0;
                    d2tot[0][1] += d2y0;
                    d2tot[0][2] += d2z0;
                    d2tot[0][4] += d2thetay0;
                    d2tot[0][5] += d2thetaz0;

                    
                    drijRot.E(drij);
                    rot1.transform(drijRot);
                    ga.Ea1Tv1(adu, drijRot);
                    double d2xyaa = drijRot.getX(0)*drijRot.getX(1)*adu2;
                    double d2xyab = drijRot.getX(0)*abdu2;
                    double d2xyba = drijRot.getX(1)*abdu2;
                    
                    r2a.E(drijRot);
                    r2a.TE(drijRot);
                    r2ab.E(drijRot);
                    drijRot.E(bdrij);
                    rot1.transform(drijRot);
                    gb.Ea1Tv1(bdu, drijRot);
                    d2xyab *= drijRot.getX(1);
                    d2xyba *= drijRot.getX(0);
                    double d2xybb = drijRot.getX(0)*drijRot.getX(1)*bdu2;
                    r2ab.TE(drijRot);
                    r2b.E(drijRot);
                    r2b.TE(drijRot);
                    

                    dxa = dxb = 0;
                    double dya = 0, dyb = 0;
                    if (ib1>0) {
                        // hydrogen
                        dxa = sitesOH - cmx;
                        dya = sitesHH;
                        if (ib1==2) dya = -dya;
                    }
                    else {
                        // oxygen
                        dxa = -cmx;
                    }
                    if (ib2>0) {
                        // hydrogen
                        dxb = sitesOH - cmx;
                        dyb = sitesHH;
                        if (ib2==2) dyb = -dyb;
                    }
                    else {
                        // oxygen
                        dxb = -cmx;
                    }

                    d2a = adu + r2a.getX(0)*adu2;
                    d2ab = r2ab.getX(0)*abdu2;
                    d2b = bdu + r2b.getX(0)*bdu2;
                    
                    double d2x1 = d2a + 2*d2ab + d2b;
//                    System.out.println(energy+" "+(ga.getX(0)+gb.getX(0))+" "+d2x1);
                    double d2thetaz1 = d2a*dya*dya + 2*d2ab*dya*dyb + d2b*dyb*dyb - ga.getX(0)*dxa - gb.getX(0)*dxb;

                    d2a = adu + r2a.getX(1)*adu2;
                    d2ab = r2ab.getX(1)*abdu2;
                    d2b = bdu + r2b.getX(1)*bdu2;
                    
                    double d2y1 = d2a + 2*d2ab + d2b;
                    d2thetaz1 += d2a*dxa*dxa + 2*d2ab*dxa*dxb + d2b*dxb*dxb - ga.getX(1)*dya - gb.getX(1)*dyb;
                    
                    d2thetaz1 -= 2*d2xyaa*dxa*dya + 2*d2xyab*dxb*dya + 2*d2xyba*dxa*dyb + 2*d2xybb*dxb*dyb;

                    d2a = adu + r2a.getX(2)*adu2;
                    d2ab = r2ab.getX(2)*abdu2;
                    d2b = bdu + r2b.getX(2)*bdu2;
                    double d2thetax1 = d2a*dya*dya + 2*d2ab*dya*dyb + d2b*dyb*dyb - ga.getX(1)*dya - gb.getX(1)*dyb;
                    double d2thetay1 = d2a*dxa*dxa + 2*d2ab*dxa*dxb + d2b*dxb*dxb - ga.getX(0)*dxa - gb.getX(0)*dxb;
                    
                    double d2z1 = d2a + 2*d2ab + d2b;

                    d2tot[1][0] += d2x1;
                    d2tot[1][1] += d2y1;
                    d2tot[1][2] += d2z1;
                    d2tot[1][3] += d2thetax1;
                    d2tot[1][4] += d2thetay1;
                    d2tot[1][5] += d2thetaz1;

                    d2tsum += 0.5*(d2x0 + d2y0 + d2z0)/mass0;
                    d2tsum += 0.5*(d2x1 + d2y1 + d2z1)/mass1;
                    d2rsum += 0.5*(d2thetay0 + d2thetaz0)/moment0;
                    d2rsum += 0.5*(d2thetax1/moment1.getX(0) + d2thetay1/moment1.getX(1) + d2thetaz1/moment1.getX(2));
                }
                else {
                    energy += rval1*c1;
//                    double x = sitePos[1][ib1].getX(0);
//                    energy += x*x;
                    double rdu = -ir1*c1*rval1;
                    double r2d2u = ir1*(ir1+1)*c1*rval1;
//                    System.out.println(rr1+" "+rval1*c1+" "+rdu/rr1+" "+r2d2u/(rr1*rr1));
                    drijRot.E(drij);
                    rot0.transform(drijRot);

                    d2.E(drijRot);
                    d2.TE(d2);
                    double d2fac = (rdu - r2d2u)/(r21*r21);
                    d2.TE(-d2fac);
                    d2.PE(rdu/r21);
                    d2tsum += 0.5*(d2.getX(0) + d2.getX(1) + d2.getX(2))/mass0;
                    d2tot[0][0] += d2.getX(0);
                    d2tot[0][1] += d2.getX(1);
                    d2tot[0][2] += d2.getX(2);
                    
                    if (ia1>0) {
                        // oxygen
                        // d2[1] = d2u/dthetaz
                        // d2[2] = d2u/dthetay
                        double dx = (2*ia1-3)*sitesCO2L;

                        d2.TE(dx*dx);
                        double gdx = drijRot.getX(0)*rdu/r21*dx;
                        double d2thetay = d2.getX(2) - gdx;
                        double d2thetaz = d2.getX(1) - gdx;
                        
                        // d2ur0[0] = 0
                        d2rsum += 0.5*(d2thetay + d2thetaz)/moment0;
                        d2tot[0][4] += d2thetay;
                        d2tot[0][5] += d2thetaz;
                    }
                    
                    drijRot.E(drij);
                    rot1.transform(drijRot);

                    d2.E(drijRot);
                    d2.TE(d2);
                    d2.TE(-d2fac);
                    d2.PE(rdu/r21);
                    double d2xy = -d2fac*drijRot.getX(0)*drijRot.getX(1);
                    // the sum here is the same as the sum above for CO2, except for the mass
                    // but we need water's d2 anyway for the rotational derivatives below.
                    d2tsum += 0.5*(d2.getX(0) + d2.getX(1) + d2.getX(2))/mass1;
                    d2tot[1][0] += d2.getX(0);
                    d2tot[1][1] += d2.getX(1);
                    d2tot[1][2] += d2.getX(2);
                    
                    double dx = 0, dy = 0;
                    if (ib1>0) {
                        // hydrogen
                        dx = sitesOH - cmx;
                        dy = sitesHH;
                        if (ib1==2) dy = -dy;
                    }
                    else {
                        // oxygen
                        dx = -cmx;
                    }
                    drijRot.TE(rdu/r21);
//                    System.out.println("g "+drijRot+" "+(drijRot.getX(1)*dy+drijRot.getX(0)*dx)+" "+sitePos[1][ib1]);
                    double d2thetax = d2.getX(2)*dy*dy - drijRot.getX(1)*dy;
                    double d2thetay = d2.getX(2)*dx*dx - drijRot.getX(0)*dx;
//                    System.out.println("t "+rdu/r21*(-drijRot.getX(0)*dy+drijRot.getX(1)*dx)+" "+sitePos[1][1]+" "+dx+" "+dy);
//                    System.out.println("t1 "+rdu/r21*(-drijRot.getX(0)*dy)+" "+(d2.getX(0)*dy*dy+drijRot.getX(0)*dx*rdu/r21));
//                    System.out.println("t2 "+rdu/r21*(+drijRot.getX(1)*dx)+" "+(d2.getX(1)*dx*dx+drijRot.getX(1)*dy*rdu/r21));
//                    System.out.println("t "+rdu/r21*drijRot.getX(0)+" "+(d2.getX(0)*dy));
//                    System.out.println("x "+dx+" "+(-dy));
//                    System.out.println("gxy "+drijRot.getX(0)+" "+drijRot.getX(1));
//                    System.out.println((drijRot.getX(1)*dy+drijRot.getX(0)*dx));
//                    System.out.println(-d2.getX(0)+" "+(2*d2xy)+" "+(-d2.getX(1)));
                    double d2thetaz = (d2.getX(0)*dy*dy - 2*d2xy*dx*dy + d2.getX(1)*dx*dx) - (drijRot.getX(1)*dy+drijRot.getX(0)*dx);
//                    System.out.println(d2thetaz);
//                    if (i==3) System.exit(1);
                    d2rsum += 0.5*(d2thetax/moment1.getX(0) + d2thetay/moment1.getX(1) + d2thetaz/moment1.getX(2));
                    d2tot[1][3] += d2thetax;
                    d2tot[1][4] += d2thetay;
                    d2tot[1][5] += d2thetaz;
                }
            }
            for (int i=0; i<2; i++) {
                for (int j=0; j<3; j++) {
                    d2tot[i][j] *= bohrConv*bohrConv;
                }
            }
            if (energy > 10000) {
                return Double.POSITIVE_INFINITY;
            }
            if (d2tsum == Double.POSITIVE_INFINITY || d2rsum == Double.POSITIVE_INFINITY) { return Double.POSITIVE_INFINITY; }
//            System.out.println(energy+" "+fac*d2tsum*bohrConv*bohrConv+" "+fac*d2rsum*bohrConv*bohrConv);
//            System.out.println((energy + fac*(d2tsum* + d2rsum)*bohrConv*bohrConv));
            return energy + fac*(d2tsum + d2rsum)*bohrConv*bohrConv;
        }
    }

    public class P2CO2H2OSCTI implements IPotentialAtomic {

        protected final Vector[] gi, ti;
        protected final double mass0;
        protected final double moment0;
        protected final Vector moment1;
        protected final double temperature, fac;
        protected final Vector ga, gb, ta, tb, tTmp;
        protected final Vector drijRot;
        protected final Vector bdrij;

        public P2CO2H2OSCTI(double temperature) {
            gi = new Vector[2];
            gi[0] = space.makeVector();
            gi[1] = space.makeVector();
            ti = new Vector[2];
            ti[0] = space.makeVector();
            ti[1] = space.makeVector();
            drijRot = space.makeVector();
            mass0 = Carbon.INSTANCE.getMass() + 2*Oxygen.INSTANCE.getMass();
            moment0 = 2*Oxygen.INSTANCE.getMass()*sitesCO2L*sitesCO2L;

            moment1 = space.makeVector();
            double cm1x = 2*Hydrogen.INSTANCE.getMass()*sitesOH/mass1;
            // 0: around O-H bisector
            moment1.setX(0, 2*Hydrogen.INSTANCE.getMass()*(sitesHH*sitesHH));
            // 1: around H-H axis
            moment1.setX(1, 2*Hydrogen.INSTANCE.getMass()*(sitesOH-cm1x)*(sitesOH-cm1x) + Oxygen.INSTANCE.getMass()*cm1x*cm1x);
            double cmH2 = (sitesOH-cm1x)*(sitesOH-cm1x) + (sitesHH*sitesHH);
            // 2: in-plane rotation
            moment1.setX(2, 2*Hydrogen.INSTANCE.getMass()*cmH2 + Oxygen.INSTANCE.getMass()*cm1x*cm1x);
            
            this.temperature = temperature;
            double hbar = Constants.PLANCK_H/(2*Math.PI);
            fac = hbar*hbar/(24/2)/(temperature*temperature);
            bdrij = space.makeVector();
            ga = space.makeVector();
            gb = space.makeVector();
            ta = space.makeVector();
            tb = space.makeVector();
            tTmp = space.makeVector();
        }
        
        public double getRange() {
            return Double.POSITIVE_INFINITY;
        }

        public void setBox(Box box) {
            
        }

        public int nBody() {
            return 2;
        }

        public double energy(IAtomList atoms) {
            IAtomOriented atom0 = (IAtomOriented)atoms.getAtom(0);
            IAtomOriented atom1 = (IAtomOriented)atoms.getAtom(1);
            Vector cm0 = atom0.getPosition();
            Vector cm1 = atom1.getPosition();
            double bohrConv = BohrRadius.UNIT.fromSim(1);
            for (int i=0; i<3; i++) {
                sitePos[0][i].Ea1Tv1(bohrConv, cm0);
                sitePos[1][i].Ea1Tv1(bohrConv, cm1);
            }

            //CO2
            Vector orCO2 = atom0.getOrientation().getDirection();
            sitePos[0][1].PEa1Tv1(+sitesCO2L, orCO2);
            sitePos[0][2].PEa1Tv1(-sitesCO2L, orCO2);
            Vector orH2O1 = atom1.getOrientation().getDirection();
            Vector orH2O2 = ((OrientationFull3D)atom1.getOrientation()).getSecondaryDirection();
            double cmx = 2*Hydrogen.INSTANCE.getMass()*sitesOH/mass1;
            sitePos[1][0].PEa1Tv1(-cmx, orH2O1);
            sitePos[1][1].PEa1Tv1(+sitesOH-cmx, orH2O1);
            sitePos[1][2].PEa1Tv1(+sitesOH-cmx, orH2O1);
            sitePos[1][1].PEa1Tv1(+sitesHH, orH2O2);
            sitePos[1][2].PEa1Tv1(-sitesHH, orH2O2);
            
            getPerp(orCO2, or01, or02);
            allOr0[0] = orCO2;
            rot0.E(allOr0);
            rot0.invert();

            orH2O3.E(orH2O1);
            orH2O3.XE(orH2O2);
            allOr1[0] = orH2O1;
            allOr1[1] = orH2O2;
            rot1.E(allOr1);
            rot1.invert();
            
            for (int i=0; i<2; i++) {
            	gi[i].E(0);
            	ti[i].E(0);
            }
            
            CO = Math.sqrt(sitePos[0][0].Mv1Squared(sitePos[1][0]));
            if (CO < 4 && !debugging) return Double.POSITIVE_INFINITY;
            double energy = 0;
            for (int i=0; i<iparams.length; i++) {
                int ia1 = iparams[i][0]-1;
                int ib1 = iparams[i][1]-3-1;
                int ia2 = iparams[i][2]-1;
                int ib2 = iparams[i][3]-3-1;
                int ir1 = iparams[i][4];
                int ir2 = iparams[i][5];
                double c1 = fparams[i];
                
                drij.Ev1Mv2(sitePos[1][ib1], sitePos[0][ia1]);
                double r21 = drij.squared();
                double rr1 = Math.sqrt(r21);
                double rval1 = 1.0/Math.pow(rr1, ir1);
                if (rr1 < core) {
                    return Double.POSITIVE_INFINITY;
                }
                if (ia2>-1) {
//                    System.out.println(ia1+" "+ib1+" "+ia2+" "+ib2+" "+ir1+" "+ir2);
                    bdrij.Ev1Mv2(sitePos[1][ib2], sitePos[0][ia2]);
                    double r22 = bdrij.squared();
                    double rr2 = Math.sqrt(r22);
                    if (rr2 < core) {
                        return Double.POSITIVE_INFINITY;
                    }
                    double rval2 = 1.0/Math.pow(rr2, ir2);
                    energy += rval1*rval2*c1;
                    
                    double adu = -c1*ir1*rval1*rval2/r21;
                    double bdu = -c1*ir2*rval1*rval2/r22;
                    
                    drijRot.E(drij);
                    rot0.transform(drijRot);
                    ga.Ea1Tv1(adu, drijRot);
                    gi[0].PE(ga);

                    drijRot.E(bdrij);
                    rot0.transform(drijRot);
                    gb.Ea1Tv1(bdu, drijRot);
                    gi[0].PE(gb);

                    double dxa = 0;
                    if (ia1>0) {
                        dxa = (2*ia1-3)*sitesCO2L;
                        tTmp.E(0);
                        tTmp.setX(0, dxa);
                        tTmp.XE(ga);
                        ti[0].PE(tTmp);
                    }
                    double dxb = 0;
                    if (ia2>0) {
                        dxb = (2*ia2-3)*sitesCO2L;
                        tTmp.E(0);
                        tTmp.setX(0, dxb);
                        tTmp.XE(gb);
                        ti[0].PE(tTmp);
                    }

                    drijRot.E(drij);
                    rot1.transform(drijRot);
                    ga.Ea1Tv1(adu, drijRot);
                    gi[1].PE(ga);
                    
                    drijRot.E(bdrij);
                    rot1.transform(drijRot);
                    gb.Ea1Tv1(bdu, drijRot);
                    gi[1].PE(gb);

                    dxa = dxb = 0;
                    double dya = 0, dyb = 0;
                    if (ib1>0) {
                        // hydrogen
                        dxa = sitesOH - cmx;
                        dya = sitesHH;
                        if (ib1==2) dya = -dya;
                    }
                    else {
                        // oxygen
                        dxa = -cmx;
                    }
                    if (ib2>0) {
                        // hydrogen
                        dxb = sitesOH - cmx;
                        dyb = sitesHH;
                        if (ib2==2) dyb = -dyb;
                    }
                    else {
                        // oxygen
                        dxb = -cmx;
                    }

                    tTmp.setX(0, dxa);
                    tTmp.setX(1, dya);
                    tTmp.setX(2, 0);
                    tTmp.XE(ga);
                    ti[1].PE(tTmp);

                    tTmp.setX(0, dxb);
                    tTmp.setX(1, dyb);
                    tTmp.setX(2, 0);
                    tTmp.XE(gb);
                    ti[1].PE(tTmp);
                }
                else {
                    energy += rval1*c1;
//                    double x = sitePos[1][ib1].getX(0);
//                    energy += x*x;
                    double rdu = -ir1*c1*rval1;
//                    System.out.println(rr1+" "+rval1*c1+" "+rdu/rr1+" "+r2d2u/(rr1*rr1));
                    drijRot.E(drij);
                    rot0.transform(drijRot);
                    ga.Ea1Tv1(rdu/r21, drijRot);
                    gi[0].PE(ga);
                    
                    if (ia1>0) {
                        // oxygen
                        // d2[1] = d2u/dthetaz
                        // d2[2] = d2u/dthetay
                        double dx = (2*ia1-3)*sitesCO2L;

                        tTmp.E(0);
                        tTmp.setX(0, dx);
                        tTmp.XE(ga);
                        ti[0].PE(tTmp);
                    }
                    
                    drijRot.E(drij);
                    rot1.transform(drijRot);

                    double dx = 0, dy = 0;
                    if (ib1>0) {
                        // hydrogen
                        dx = sitesOH - cmx;
                        dy = sitesHH;
                        if (ib1==2) dy = -dy;
                    }
                    else {
                        // oxygen
                        dx = -cmx;
                    }
                    drijRot.TE(rdu/r21);
                    gi[1].PE(drijRot);
                    tTmp.setX(0, dx);
                    tTmp.setX(1, dy);
                    tTmp.setX(2, 0);
                    tTmp.XE(drijRot);
                    ti[1].PE(tTmp);
                }
            }
            if (energy > 10000) {
                return Double.POSITIVE_INFINITY;
            }
            double sum = 0;
            sum += 0.5*gi[0].squared()/mass0;
            sum += 0.5*gi[1].squared()/mass1;
            sum += 0.5*ti[0].squared()/moment0;
            for (int i=0; i<3; i++) {
            	double tx = ti[1].getX(i);
            	sum += 0.5*tx*tx/moment1.getX(i);
            }
            sum *= bohrConv*bohrConv;
//            System.out.println(energy+" "+fac*sum);
            return energy + fac*sum;
        }
    }

    public double virial(IAtomList atoms) {
        return 0;
    }

    public Vector[] gradient(IAtomList atoms) {
        return null;
    }

    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return null;
    }
    
    /**
     * Reads geometry.txt file and computes energy for each configuration
     */
    public static void mainTestU(String[] args) throws IOException {
//        System.out.println(Hartree.UNIT.toSim(0.1));
//        System.exit(1);
        Space space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        SpeciesSpheresRotating speciesCO2 = new SpeciesSpheresRotating(space, new ElementSimple("CO2", Carbon.INSTANCE.getMass()+2*Oxygen.INSTANCE.getMass()));
        sim.addSpecies(speciesCO2);
        SpeciesSpheresRotating speciesH2O = new SpeciesSpheresRotating(space, new ElementSimple("H2O", Oxygen.INSTANCE.getMass()+2*Hydrogen.INSTANCE.getMass()));
        speciesH2O.setAxisSymmetric(false);
        sim.addSpecies(speciesH2O);
        Box box = new etomica.box.Box(space);
        sim.addBox(box);
        box.setNMolecules(speciesCO2, 1);
        box.setNMolecules(speciesH2O, 1);
        P2CO2H2OWheatley p2 = new P2CO2H2OWheatley(space);
        p2.setBox(box);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{100,100,100}));
        IAtomList pair = box.getLeafList();
        IAtomOriented atom0 = (IAtomOriented)pair.getAtom(0);
        IAtomOriented atom1 = (IAtomOriented)pair.getAtom(1);
        FileReader fReader = new FileReader("geometry.txt");
        BufferedReader bufReader = new BufferedReader(fReader);
        String line = null;
        Vector x = space.makeVector(new double[]{1,0,0});
        Vector y = space.makeVector(new double[]{0,1,0});
        Vector zm = space.makeVector(new double[]{0,0,-1});
        Vector z = space.makeVector(new double[]{0,0,1});
        atom0.getOrientation().setDirection(z);
        Vector axis = space.makeVector();
        Vector p = space.makeVector();
        double cmx = 2*Hydrogen.INSTANCE.getMass()*BohrRadius.UNIT.toSim(sitesOH)/p2.mass1;
        int iline = 1;
        while ((line = bufReader.readLine()) != null) {
//            if (iline < 1305) {
//                iline++;
//                continue;
//            }
            String[] fields = line.trim().split(" +");
//            double zwb = Double.parseDouble(fields[1]);
            double xw = BohrRadius.UNIT.toSim(Double.parseDouble(fields[0]));
            double zw = BohrRadius.UNIT.toSim(Double.parseDouble(fields[1]));
//            System.out.println(zwb+" => "+zw);
            double alpha = Double.parseDouble(fields[2]);
            double cosBeta = Double.parseDouble(fields[3]);
            double gamma = Double.parseDouble(fields[4]);
            ((OrientationFull3D)atom1.getOrientation()).setDirections(zm, x);
            ((OrientationFull3D)atom1.getOrientation()).rotateBy(alpha, zm);
            axis.E(atom1.getOrientation().getDirection());
            axis.XE(((OrientationFull3D)atom1.getOrientation()).getSecondaryDirection());
            ((OrientationFull3D)atom1.getOrientation()).rotateBy(-Math.acos(cosBeta), axis);
            axis.E(atom1.getOrientation().getDirection());
            ((OrientationFull3D)atom1.getOrientation()).rotateBy(gamma, axis);
            
            p.E(new double[]{xw,0,zw});
//            System.out.println("p "+p);
            p.PEa1Tv1(cmx, atom1.getOrientation().getDirection());
            atom1.getPosition().E(p);
//            ((OrientationFull3D)atom1.getOrientation()).setDirections(z,x);
//            System.out.println("H2O or1 "+((OrientationFull3D)atom1.getOrientation()).getDirection());
//            System.out.println("H2O or2 "+((OrientationFull3D)atom1.getOrientation()).getSecondaryDirection());
//            System.out.println("=> "+p);
            double u = p2.energy(pair);
            System.out.println(Hartree.UNIT.fromSim(u));
//            if (iline==1) {
//                System.out.println(BohrRadius.UNIT.fromSim(xw)+" "+BohrRadius.UNIT.fromSim(zw)+" "+alpha+" "+cosBeta+" "+gamma);
//                break;
//            }
            iline++;
        }
        bufReader.close();
    }

    /**
     * Check charges (all sums should be 0 for neutral molecules)
     */
    public static void mainQ(String[] args) {
        double[] c = P2CO2H2OWheatley.fparams;
        for (int i=0; i<3; i++) {
            double sum = 0;
            for (int j=0; j<3; j++) {
                sum += c[i*9+j*3];
            }
            System.out.println("i "+i+" "+Hartree.UNIT.fromSim(sum));
        }
        for (int j=0; j<3; j++) {
            double sum = 0;
            for (int i=0; i<3; i++) {
                sum += c[i*9+j*3];
            }
            System.out.println("j "+j+" "+Hartree.UNIT.fromSim(sum));
        }
    }
    
    /**
     * Randomly moves molcules all over, computing and checking 2nd derivative
     */
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        double temperature = Kelvin.UNIT.toSim(200);
        Simulation sim = new Simulation(space);
        SpeciesSpheresRotating speciesCO2 = new SpeciesSpheresRotating(space, new ElementSimple("CO2", Carbon.INSTANCE.getMass()+2*Oxygen.INSTANCE.getMass()));
        sim.addSpecies(speciesCO2);
        SpeciesSpheresRotating speciesH2O = new SpeciesSpheresRotating(space, new ElementSimple("H2O", Oxygen.INSTANCE.getMass()+2*Hydrogen.INSTANCE.getMass()));
        speciesH2O.setAxisSymmetric(false);
        sim.addSpecies(speciesH2O);
        Box box = new etomica.box.Box(space);
        sim.addBox(box);
        box.setNMolecules(speciesCO2, 1);
        box.setNMolecules(speciesH2O, 1);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{100,100,100}));
        IAtomList pair = box.getLeafList();
        IAtomOriented atom0 = (IAtomOriented)pair.getAtom(0);
        IAtomOriented atom1 = (IAtomOriented)pair.getAtom(1);
//        ((IAtomOriented)pair.getAtom(0)).getOrientation().setDirection(space.makeVector(new double[]{Math.cos(22.5/180.0*Math.PI), Math.sin(22.5/180.0*Math.PI),0}));
        Vector o1 = space.makeVector(new double[]{0,0,-1});
        Vector o2 = space.makeVector(new double[]{1,0,0});
//        ((OrientationFull3D)atom1.getOrientation()).setDirections(o1, o2);
        P2CO2H2OWheatley p2 = new P2CO2H2OWheatley(space);
        P2CO2H2OSC p2SC = p2.makeSemiclassical(temperature);
        System.out.println("or: "+atom0.getOrientation().getDirection()+" "+atom1.getOrientation().getDirection()+" "+((OrientationFull3D)atom1.getOrientation()).getSecondaryDirection());
        double lu = 0, lg = 0;
        double dx = 0.0001;
        Vector x = null;
        Vector y = space.makeVector(new double[]{0,1,0});
        Vector z = space.makeVector(new double[]{0,0,1});
        atom1.getPosition().setX(0, 5);
        double cmx = 2*Hydrogen.INSTANCE.getMass()*sitesOH/p2.mass1;
//        ((OrientationFull3D)atom1.getOrientation()).rotateBy(-Math.atan2(p2.sitesHH*0.5,p2.sitesOH-cmx), z);
//        ((OrientationFull3D)atom1.getOrientation()).rotateBy(Math.PI/2, y);
        Vector vdx = space.makeVector();
//        ((OrientationFull3D)atom1.getOrientation()).rotateBy(0.28, z);

        RandomMersenneTwister random = new RandomMersenneTwister(RandomNumberGeneratorUnix.getRandSeedArray());
        Vector[] xyzAxes = new Vector[]{x,y,z};
        double[] u = new double[3];
        for (int j=0; j<30; j++) {
            int imol = random.nextInt(2);
//            imol = 0;
            IAtomOriented iAtom = imol == 0 ? atom0 : atom1;
            x = iAtom.getOrientation().getDirection();
            xyzAxes[0] = x;
            if (imol==0) {
                p2.getPerp(x, y, z);
            }
            else {
                y.E(((OrientationFull3D)iAtom.getOrientation()).getSecondaryDirection());
                z.E(x);
                z.XE(y);
            }
//            System.out.println(x+" "+y+" "+z);
            boolean rot = random.nextInt(2) == 0;
//            rot = true;
            double d2 = 0;
            int xyz = random.nextInt(3);
//            xyz = 2;
            vdx.E(xyzAxes[xyz]);
            if (random.nextInt(2)==0) {
                if (rot) {
                    ((IOrientation3D)iAtom.getOrientation()).rotateBy(0.5, vdx);
                }
                else {
                    iAtom.getPosition().PEa1Tv1(0.5, vdx);
                }
                continue;
            }
            if (rot) {
//                System.out.println("vdx "+vdx);
//                if (imol==1 && xyz==0) dx *= 0.1;
                for (int i=0; i<3; i++) {
                    ((IOrientation3D)iAtom.getOrientation()).rotateBy(dx, vdx);
                    u[i] = p2.energy(pair);
//                    if (imol==1 && xyz==0) System.out.println(u[i]);
                    if (i==1) {
                        p2SC.energy(pair);
                        d2 = p2SC.d2tot[imol][3+xyz];
                    }
                }
            }
            else {
                for (int i=0; i<3; i++) {
                    iAtom.getPosition().PEa1Tv1(dx, vdx);
                    u[i] = p2.energy(pair);
                    if (i==1) {
                        p2SC.energy(pair);
                        d2 = p2SC.d2tot[imol][xyz];
                    }
                }
            }
            double d2fd = (u[0] - 2*u[1] + u[2])/(dx*dx);
            if (d2==0 && d2fd==0) continue;
            double check = (d2fd-d2)/(0.5*(d2+d2fd));
            System.out.print(String.format("%d %d %d %d %+10.4e %+10.4e %10.4e", j, imol, rot?1:0, xyz, d2, d2fd, check));
            if (Math.abs(check) > 0.01) {
                System.out.println("  oops");
                if (imol==1 || !rot) throw new RuntimeException("oops");
            }
            else {
                System.out.println();
            }
        }
    }

    public static void main2(String[] args) {
        Space space = Space3D.getInstance();
        double temperature = Kelvin.UNIT.toSim(200);
        Simulation sim = new Simulation(space);
        SpeciesSpheresRotating speciesCO2 = new SpeciesSpheresRotating(space, new ElementSimple("CO2", Carbon.INSTANCE.getMass()+2*Oxygen.INSTANCE.getMass()));
        sim.addSpecies(speciesCO2);
        SpeciesSpheresRotating speciesH2O = new SpeciesSpheresRotating(space, new ElementSimple("H2O", Oxygen.INSTANCE.getMass()+2*Hydrogen.INSTANCE.getMass()));
        speciesH2O.setAxisSymmetric(false);
        sim.addSpecies(speciesH2O);
        Box box = new etomica.box.Box(space);
        sim.addBox(box);
        box.setNMolecules(speciesCO2, 1);
        box.setNMolecules(speciesH2O, 1);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{100,100,100}));
        IAtomList pair = box.getLeafList();
        IAtomOriented atom0 = (IAtomOriented)pair.getAtom(0);
        IAtomOriented atom1 = (IAtomOriented)pair.getAtom(1);
        Vector p1 = atom1.getPosition();
        IOrientation or0 = atom0.getOrientation();
        OrientationFull3D or1 = (OrientationFull3D)atom1.getOrientation();
        P2CO2H2OWheatley p2 = new P2CO2H2OWheatley(space);
        P2CO2H2OSCTI p2SC = p2.makeSemiclassicalTI(temperature);
        P2CO2H2OSCTI p2CSC = p2.makeSemiclassicalTI(Double.POSITIVE_INFINITY);
//        ((OrientationFull3D)atom1.getOrientation()).setDirections(o1, o2);
        System.out.println("OCO: "+or0.getDirection());

        double cmx = 2*Hydrogen.INSTANCE.getMass()*sitesOH/p2.mass1;
        p1.setX(0, BohrRadius.UNIT.toSim(9+cmx));
        double uc = p2.energy(pair);
        double usc = p2SC.energy(pair);
        double ucsc = p2CSC.energy(pair);
        System.out.println("rO:  9,0,0");
        System.out.println("O-H: "+or1.getDirection());
        System.out.println("H-H: "+or1.getSecondaryDirection());
        System.out.println(String.format("uc:     %+22.15e", uc));
        System.out.println(String.format("uTI:    %+22.15e", usc));
        System.out.println(String.format("uTI-uc: %+22.15e", usc-uc));

        System.out.println();
        p1.setX(0, BohrRadius.UNIT.toSim(6));
        p1.setX(1, BohrRadius.UNIT.toSim(6));
        p1.setX(2, BohrRadius.UNIT.toSim(cmx));
        or1.setDirections(space.makeVector(new double[]{0,0,1}),space.makeVector(new double[]{1,0,0}));
        uc = p2.energy(pair);
        usc = p2SC.energy(pair);
        ucsc = p2CSC.energy(pair);
        System.out.println("rO:  6,6,0");
        System.out.println("O-H: "+or1.getDirection());
        System.out.println("H-H: "+or1.getSecondaryDirection());
        System.out.println(String.format("uc:     %+22.15e", uc));
        System.out.println(String.format("uTI:    %+22.15e", usc));
        System.out.println(String.format("uTI-uc: %+22.15e", usc-uc));

    }

    /**
     * Writes x-y grid of energies
     */
    public static void main2DGrid(String[] args) throws IOException {
        Space space = Space3D.getInstance();
        double temperature = Kelvin.UNIT.toSim(1);
        Simulation sim = new Simulation(space);
        SpeciesSpheresRotating speciesCO2 = new SpeciesSpheresRotating(space, new ElementSimple("CO2", Carbon.INSTANCE.getMass()+2*Oxygen.INSTANCE.getMass()));
        sim.addSpecies(speciesCO2);
        SpeciesSpheresRotating speciesH2O = new SpeciesSpheresRotating(space, new ElementSimple("H2O", Oxygen.INSTANCE.getMass()+2*Hydrogen.INSTANCE.getMass()));
        speciesH2O.setAxisSymmetric(false);
        sim.addSpecies(speciesH2O);
        Box box = new etomica.box.Box(space);
        sim.addBox(box);
        box.setNMolecules(speciesCO2, 1);
        box.setNMolecules(speciesH2O, 1);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{100,100,100}));
        IAtomList pair = box.getLeafList();
        IAtomOriented atom0 = (IAtomOriented)pair.getAtom(0);
        IAtomOriented atom1 = (IAtomOriented)pair.getAtom(1);
        Vector p0 = atom1.getPosition();
        Vector p1 = atom1.getPosition();
        IOrientation or0 = atom0.getOrientation();
        OrientationFull3D or1 = (OrientationFull3D)atom1.getOrientation();
        P2CO2H2OWheatley p2 = new P2CO2H2OWheatley(space);
        P2CO2H2OSC p2SCFH = p2.makeSemiclassical(temperature);
        P2CO2H2OSCTI p2SCTI = p2.makeSemiclassicalTI(temperature);
//        ((OrientationFull3D)atom1.getOrientation()).setDirections(o1, o2);
        System.out.println("OCO: "+or0.getDirection());
        System.out.println("OM: "+or1.getDirection());
        System.out.println("HH: "+or1.getSecondaryDirection());
        FileWriter fw = new FileWriter("/tmp/foo.dat");
        
        double dxy = 0.15;
        for (int ix=-70; ix<=70; ix++) {
            double x = dxy*ix;
            p0.setX(0, x);
            for (int iy=-70; iy<=70; iy++) {
                double y = dxy*iy;
                p0.setX(1, y);
                double uc = p2.energy(pair);
                double usc = p2SCFH.energy(pair);
                double uscTI = p2SCTI.energy(pair);
                fw.write(String.format("{%2.2f,%2.2f,%3.3f,%3.3f,%3.3f},\n",x,y,Kelvin.UNIT.fromSim(uc),Kelvin.UNIT.fromSim(usc-uc),Kelvin.UNIT.fromSim(uscTI-uc)));
//                fw.write(String.format("%2.2f %2.2f %3.3f %3.3f %3.3f \n",x,y,Kelvin.UNIT.fromSim(uc),Kelvin.UNIT.fromSim(usc),Kelvin.UNIT.fromSim(uscTI)));

            }
        }
        fw.close();
        
    }

}
