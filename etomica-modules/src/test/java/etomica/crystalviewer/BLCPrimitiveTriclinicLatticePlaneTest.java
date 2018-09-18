/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.crystalviewer;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.lattice.crystal.PrimitiveTriclinic;
import etomica.space.Vector;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.fail;

public class BLCPrimitiveTriclinicLatticePlaneTest {

	private String funcName = "";

	private double epsilon = 1.0E-5;;

	private final int DEFAULT_SIZE = 7;
	private final int DEFAULT_MILLER[] = {0,0,1};
	private final double DEFAULT_ALPHA = Math.PI * 2 * (90.0 / 360.0);
	private final double DEFAULT_BETA = Math.PI * 2 * (90.0 / 360.0);
	private final double DEFAULT_GAMMA = Math.PI * 2 * (90.0 / 360.0);
	private final int DEFAULT_BOX[] = {DEFAULT_SIZE, DEFAULT_SIZE, DEFAULT_SIZE};

	private LatticePlaneTestUtility lptu = null;

	@BeforeEach
	public void setUp() throws Exception {
		if (lptu == null) {
			lptu = new LatticePlaneTestUtility();			
	        lptu.createLatticeAndBox(lptu.TRICLINIC, DEFAULT_MILLER, DEFAULT_BOX);
	        lptu.setDimensions(DEFAULT_SIZE);
		}
	}

	@AfterEach
	public void tearDown() throws Exception {
	}

	private double[] makeArray(Vector v) {
	    return new double[] {v.getX(0), v.getX(1), v.getX(2)};
	}

    /*
     * Miller indices = 0, 0, 1
     * size of cell (A,B,C) = 1.0
     * cells per side = 7
     * plane = 1.0
     * alpha = 90 degrees
     * beta = 90 degrees
     * gamma = 90 degrees
     */
    @Test
	public void testStandard() {

    	int idx = 0;
    	double cubicSize = 1.0;
    	double plane = 1.0;
    	IAtomList leafList = null;

        lptu.createLatticeAndBox(lptu.TRICLINIC, DEFAULT_MILLER, DEFAULT_BOX);

        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeA(cubicSize);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeB(cubicSize);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeC(cubicSize);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleAlpha(DEFAULT_ALPHA);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleBeta(DEFAULT_BETA);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleGamma(DEFAULT_GAMMA);
        lptu.setDimensions(DEFAULT_SIZE);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());
        double spacePos = lptu.getLatticePlaneSpacePosition();

    	leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a =  leafList.get(idx);
                if(a.getPosition().getX(2) >= spacePos-epsilon &&
                   a.getPosition().getX(2) <= spacePos+epsilon) {
            	    Assertions.assertTrue(lptu.getLatticePlane().inPlane(a.getPosition()));
                }
                else {
            	    Assertions.assertFalse(lptu.getLatticePlane().inPlane(a.getPosition()));
                }
		    }
		}
        catch (org.opentest4j.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            if(a.getPosition().getX(2) >= spacePos-epsilon &&
               a.getPosition().getX(2) <= spacePos+epsilon) {
                System.out.println(funcName + " -> Atom position : " + a.getPosition() +
            			" should be in plane but is not.");
            	fail();
            }
            else {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is."); 
            	fail();
            }
        }

    } // End testStandard()

    /*
     * Miller indices = 0, 0, 1
     * size of cell (A) = 1.5
     * size of cell (B) = 1.75
     * size of cell (C) = 2.0
     * cells per side = 7
     * plane = -2.0
     * alpha = 90 degrees
     * beta = 90 degrees
     * gamma = 90 degrees
     */
    @Test
	public void testCellSizeAllDifferent() {

    	int idx = 0;
    	double cubicSizeA = 1.5;
    	double cubicSizeB = 1.75;
    	double cubicSizeC = 2.0;
    	double plane = -2.0;
    	IAtomList leafList = null;

        lptu.createLatticeAndBox(lptu.TRICLINIC, DEFAULT_MILLER, DEFAULT_BOX);

        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeA(cubicSizeA);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeB(cubicSizeB);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleAlpha(DEFAULT_ALPHA);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleBeta(DEFAULT_BETA);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleGamma(DEFAULT_GAMMA);
        lptu.setDimensions(DEFAULT_SIZE);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());
        double spacePos = lptu.getLatticePlaneSpacePosition();

    	leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a =  leafList.get(idx);

                if(a.getPosition().getX(2) >= spacePos-epsilon &&
                   a.getPosition().getX(2) <= spacePos+epsilon) {
            	    Assertions.assertTrue(lptu.getLatticePlane().inPlane(
							(a.getPosition())));
                }
                else {
            	    Assertions.assertFalse(lptu.getLatticePlane().inPlane(
							(a.getPosition())));
                }
		    }
		}
        catch (org.opentest4j.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            if(a.getPosition().getX(2) >= spacePos-epsilon &&
               a.getPosition().getX(2) <= spacePos+epsilon) {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should be in plane but is not.");
            	fail();
            }
            else {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");
            	fail();
            }
        }
    } // End testCellSizeAllDifferent()

    /*
     * Miller indices = 0, 0, 1
     * size of cell (A) = 1.25
     * size of cell (B) = 1.5
     * size of cell (C) = 1.75
     * cells per side = 7
     * plane = 3.0
     * alpha = 100.0 degrees
     * beta = 110.0 degrees
     * gamma = 120.0 degrees
     */
    @Test
	public void testCellSizeAllDifferentAnglesAllDifferent() {

    	int idx = 0;
    	double cubicSizeA = 1.25;
    	double cubicSizeB = 1.5;
    	double cubicSizeC = 1.75;
    	double plane = 3.0;
    	double alpha = Math.PI * 2.0 * (100.0 / 360.0);
    	double beta = Math.PI * 2.0 * (110.0 / 360.0);
    	double gamma = Math.PI * 2.0 * (80.0 / 360.0);
    	IAtomList leafList = null;

        lptu.createLatticeAndBox(lptu.TRICLINIC, DEFAULT_MILLER, DEFAULT_BOX);

        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeA(cubicSizeA);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeB(cubicSizeB);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleAlpha(alpha);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleBeta(beta);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleGamma(gamma);
        lptu.setDimensions(DEFAULT_SIZE);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());
        double spacePos = lptu.getLatticePlaneSpacePosition();

    	leafList = lptu.getBox().getLeafList();

       	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a =  leafList.get(idx);

                if(a.getPosition().getX(2) >= spacePos-epsilon &&
                   a.getPosition().getX(2) <= spacePos+epsilon) {
            	    Assertions.assertTrue(lptu.getLatticePlane().inPlane(a.getPosition()));
                }
                else {
            	    Assertions.assertFalse(lptu.getLatticePlane().inPlane(a.getPosition()));
                }
		    }
		}
        catch (org.opentest4j.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            if(a.getPosition().getX(2) >= spacePos-epsilon &&
               a.getPosition().getX(2) <= spacePos+epsilon) {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should be in plane but is not.");
            	fail();
            }
            else {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");
            	fail();
            }
        }

    } // End testCellSizeAllDifferentAnglesAllDifferent()

    /*
     * Miller indices = 2, 0, 2
     * size of cell (A,B,C) = 1.0
     * cells per side = 5
     * plane = 6.0
     * alpha = 100.0 degrees
     * beta = 100.0 degrees
     * gamma = 100.0 degrees

     */
    @Test
	public void testOddMillerIndices() {

    	int idx = 0;
    	double cubicSizeA = 1.1;
    	double cubicSizeB = 1.2;
    	double cubicSizeC = 1.3;
    	IAtomList leafList = null;
    	int size = 7;
    	double plane = 6.0;
    	int itemsFound = 0;
    	double alpha = Math.PI * 2.0 * (100.0 / 360.0);
    	double beta = Math.PI * 2.0 * (100.0 / 360.0);
    	double gamma = Math.PI * 2.0 * (80.0 / 360.0);
    	int[] millerIndices = new int[] { 2, 0, 2 };
        double actualPlane[][] =
           { { -1.302361332501978, -4.1135695860001675, 3.798478780989239 },
                { -1.0939835193016618, -2.9318002823855185, 3.798478780989239 },
                { -0.8856057061013454, -1.750030978770869, 3.798478780989239 },
                { -0.6772278929010289, -0.56826167515622, 3.798478780989239 },
                { -0.4688500797007116, 0.61350762845843, 3.798478780989239 },
                { -0.2604722665003951, 1.7952769320730786, 3.798478780989239 },
                { -0.05209445330007867, 2.9770462356877285, 3.798478780989239 },
                { 0.023381298465031453, -3.924149027614761, 2.53231918732616 },
                { 0.2317591116653479, -2.742379724000112, 2.53231918732616 },
                { 0.44013692486566436, -1.5606104203854627, 2.53231918732616 },
                { 0.6485147380659808, -0.37884111677081345, 2.53231918732616 },
                { 0.8568925512662968, 0.802928186843836, 2.53231918732616 },
                { 1.0652703644666133, 1.984697490458485, 2.53231918732616 },
                { 1.2736481776669297, 3.166466794073135, 2.53231918732616 },
                { 1.3491239294320403, -3.734728469229355, 1.2661595936630796 },
                { 1.5575017426323567, -2.5529591656147055, 1.2661595936630796 },
                { 1.7658795558326732, -1.3711898620000562, 1.2661595936630796 },
                { 1.9742573690329897, -0.18942055838540695, 1.2661595936630796 },
                { 2.182635182233306, 0.9923487452292425, 1.2661595936630796 },
                { 2.3910129954336226, 2.1741180488438907, 1.2661595936630796 },
                { 2.599390808633939, 3.3558873524585415, 1.2661595936630796 },
                { 2.6748665603990505, -3.5453079108439485, 0.0 },
                { 2.883244373599367, -2.363538607229299, 0.0 },
                { 3.0916221867996834, -1.1817693036146497, 0.0 },
                { 3.3, -4.440892098500626E-16, 0.0 },
                { 3.5083778132003163, 1.1817693036146482, 0.0 },
                { 3.7167556264006327, 2.363538607229297, 0.0 },
                { 3.925133439600949, 3.545307910843947, 0.0 } };

        DoubleTwoDArray dd = new DoubleTwoDArray(actualPlane);

        lptu.createLatticeAndBox(lptu.TRICLINIC, millerIndices, new int[] {size, size, size});
        
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeA(cubicSizeA);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeB(cubicSizeB);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleAlpha(alpha);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleBeta(beta);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleGamma(gamma);
        lptu.setDimensions(size);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getBox().getLeafList();

    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a = leafList.get(idx);

			    if(dd.contains(makeArray(a.getPosition())) == true) {
			    	itemsFound++;
            	    Assertions.assertTrue(lptu.getLatticePlane().inPlane(a.getPosition()));
                }
                else {
            	    Assertions.assertFalse(lptu.getLatticePlane().inPlane(a.getPosition()));
                }
		    }
		}
        catch (org.opentest4j.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            if(dd.contains(makeArray(a.getPosition()))) {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should be in plane but is not.");
            }
            else {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");
            }
        	fail();
        }

        Assertions.assertEquals(actualPlane.length, itemsFound);

    } // End testOddMillerIndices()


    /*
     * Miller indices = 1, 1, 1
     * size of cell (A) = 1.1
     * size of cell (B) = 1.2
     * size of cell (C) = 1.3
     * cells per side = 7
     * plane = 2.95
     * alpha = 97.0 degrees
     * beta = 104.0 degrees
     * gamma = 115.0 degrees
     */
    @Test
	public void testPlaneMinusFiveHundreths() {

    	int idx = 0;
    	double cubicSizeA = 1.1;
    	double cubicSizeB = 1.2;
    	double cubicSizeC = 1.3;
    	double alpha = Math.PI * 2.0 * (97.0 / 360.0);
    	double beta = Math.PI * 2.0 * (104.0 / 360.0);
    	double gamma = Math.PI * 2.0 * (115.0 / 360.0);
    	IAtomList leafList = null;
    	int size = 7;
    	double plane = 2.95;
    	int[] millerIndices = new int[] { 1, 1, 1 };

        lptu.createLatticeAndBox(lptu.TRICLINIC, millerIndices, new int[] {size, size, size});
        
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeA(cubicSizeA);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeB(cubicSizeB);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleAlpha(alpha);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleBeta(beta);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleGamma(gamma);
        lptu.setDimensions(size);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getBox().getLeafList();


    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a = leafList.get(idx);

            	Assertions.assertFalse(lptu.getLatticePlane().inPlane(a.getPosition()));
		    }
		}
        catch (org.opentest4j.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");
            fail();
        }

    } // End testPlaneMinusFiveHundreths()

    /*
     * Miller indices = 1, 1, 1
     * size of cell (A) = 1.1
     * size of cell (B) = 1.2
     * size of cell (C) = 1.3
     * cells per side = 7
     * plane = 3.0
     * alpha = 97.0 degrees
     * beta = 104.0 degrees
     * gamma = 115.0 degrees
     */
    @Test
	public void testPlane() {

    	int idx = 0;
    	double cubicSizeA = 1.1;
    	double cubicSizeB = 1.2;
    	double cubicSizeC = 1.3;
    	IAtomList leafList = null;
    	int size = 7;
    	double plane = 3.0;
    	int itemsFound = 0;
    	double alpha = Math.PI * 2.0 * (97.0 / 360.0);
    	double beta = Math.PI * 2.0 * (104.0 / 360.0);
    	double gamma = Math.PI * 2.0 * (115.0 / 360.0);
    	int[] millerIndices = new int[] { 1, 1, 1 };
        double actualPlane[][] =
          { { -5.7649211351052205, 2.2983239304216774, 3.6592047969123156 },
            { -4.157779221016382, 1.210754585977698, 3.6592047969123156 },
            { -4.350422670825653, 2.6197852980584315, 2.43946986460821 },
            { -2.5506373069275425, 0.12318524153371735, 3.6592047969123156 },
            { -2.743280756736814, 1.5322159536144522, 2.43946986460821 },
            { -2.935924206546085, 2.9412466656951857, 1.219734932304105 },
            { -0.9434953928387033, -0.9643841029102633, 3.6592047969123156 },
            { -1.1361388426479744, 0.4446466091704715, 2.43946986460821 },
            { -1.3287822924572454, 1.853677321251206, 1.219734932304105 },
            { -1.5214257422665174, 3.2627080333319407, 0.0 },
            { 0.6636465212501359, -2.051953447354243, 3.6592047969123156 },
            { 0.4710030714408646, -0.6429227352735087, 2.43946986460821 },
            { 0.27835962163159333, 0.7661079768072265, 1.219734932304105 },
            { 0.08571617182232205, 2.1751386888879614, 0.0 },
            { -0.10692727798694923, 3.584169400968695, -1.2197349323041053 },
            { 2.2707884353389747, -3.139522791798224, 3.6592047969123156 },
            { 2.078144985529703, -1.730492079717489, 2.43946986460821 },
            { 1.885501535720432, -0.32146136763675415, 1.219734932304105 },
            { 1.692858085911161, 1.0875693444439811, 0.0 },
            { 1.5002146361018895, 2.4966000565247155, -1.2197349323041053 },
            { 1.307571186292618, 3.9056307686054508, -2.4394698646082107 },
            { 3.877930349427814, -4.227092136242204, 3.6592047969123156 },
            { 3.685286899618543, -2.818061424161469, 2.43946986460821 },
            { 3.4926434498092718, -1.4090307120807342, 1.219734932304105 },
            { 3.3000000000000007, 4.440892098500626E-16, 0.0 },
            { 3.107356550190729, 1.4090307120807353, -1.2197349323041053 },
            { 2.9147131003814577, 2.8180614241614705, -2.4394698646082107 },
            { 2.7220696505721866, 4.227092136242205, -3.6592047969123156 } };
        DoubleTwoDArray dd = new DoubleTwoDArray(actualPlane);

        lptu.createLatticeAndBox(lptu.TRICLINIC, millerIndices, new int[] {size, size, size});
        
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeA(cubicSizeA);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeB(cubicSizeB);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleAlpha(alpha);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleBeta(beta);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleGamma(gamma);
        lptu.setDimensions(size);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getBox().getLeafList();


    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a = leafList.get(idx);

			    if(dd.contains(makeArray(a.getPosition())) == true) {
			    	itemsFound++;
            	    Assertions.assertTrue(lptu.getLatticePlane().inPlane(a.getPosition()));
                }
                else {
            	    Assertions.assertFalse(lptu.getLatticePlane().inPlane(a.getPosition()));
                }
		    }
		}
        catch (org.opentest4j.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            if(dd.contains(makeArray(a.getPosition()))) {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should be in plane but is not.");
            }
            else {
            	System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");
            }
         	fail();
        }

        Assertions.assertEquals(actualPlane.length, itemsFound);

    } // End testPlane()

    /*
     * Miller indices = 1, 1, 1
     * size of cell (A) = 1.1
     * size of cell (B) = 1.2
     * size of cell (C) = 1.3
     * cells per side = 7
     * plane = 3.05
     * alpha = 97.0 degrees
     * beta = 104.0 degrees
     * gamma = 115.0 degrees
     */
    @Test
	public void testPlanePlusFiveHundreths() {

    	int idx = 0;
    	double cubicSizeA = 1.1;
    	double cubicSizeB = 1.2;
    	double cubicSizeC = 1.3;
    	IAtomList leafList = null;
    	int size = 7;
    	double plane = 3.05;
    	double alpha = Math.PI * 2.0 * (97.0 / 360.0);
    	double beta = Math.PI * 2.0 * (104.0 / 360.0);
    	double gamma = Math.PI * 2.0 * (115.0 / 360.0);
    	int[] millerIndices = new int[] { 1, 1, 1 };

        lptu.createLatticeAndBox(lptu.TRICLINIC, millerIndices, new int[] {size, size, size});
        
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeA(cubicSizeA);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeB(cubicSizeB);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setSizeC(cubicSizeC);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleAlpha(alpha);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleBeta(beta);
        ((PrimitiveTriclinic)lptu.getLattice().getPrimitive()).setAngleGamma(gamma);
        lptu.setDimensions(size);
        lptu.setLatticePlanePosition(plane);

        // This needs to come after lattice changes
        lptu.getLatticePlane().setPrimitive(lptu.getLattice().getPrimitive());

    	leafList = lptu.getBox().getLeafList();


    	try {
		    for(idx = 0; idx < leafList.size(); idx++) {
			    IAtom a = leafList.get(idx);

            	Assertions.assertFalse(lptu.getLatticePlane().inPlane(a.getPosition()));
		    }
		}
        catch (org.opentest4j.AssertionFailedError e) {
		    IAtom a =  leafList.get(idx);
            System.out.println(funcName + " ->Atom position : " + a.getPosition() +
            			" should not be in plane but is.");
            fail();
        }

    } // End testPlanePlusFiveHundreths()

    public class DoubleTwoDArray {
    	private double[][] array;
    	private double epsilon = 1.0E-5;

    	public DoubleTwoDArray(double[][] in) {
    		array = in;
    	}
    	
    	public boolean contains(double[] val) {
    		boolean b = false;

    		if(array[0].length == val.length) {
    			for(int i = 0; i < array.length; i++) {
    				boolean matching = true;
    				for(int j = 0; j < array[0].length; j++) {
    					if(val[j] < array[i][j] - epsilon ||
    					   val[j] > array[i][j] + epsilon) {
    						matching = false;
    						break;
    					}
    				}
    				if(matching == true) {
    					b = true;
    					break;
    				}
    			}
    		}
            return b;
    	}
    	
    	public void setEpsilon(double e) {
    		this.epsilon = e;
    	}

    }  // End public class DoubleTwoDArray

}
