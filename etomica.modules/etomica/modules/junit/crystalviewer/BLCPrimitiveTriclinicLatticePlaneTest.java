package etomica.modules.junit.crystalviewer;

import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;
import etomica.lattice.crystal.PrimitiveTriclinic;
import etomica.space.IVector;
import junit.framework.TestCase;

public class BLCPrimitiveTriclinicLatticePlaneTest extends TestCase {

	private String funcName = "";

	private double epsilon = 1.0E-5;;

	private final int DEFAULT_SIZE = 7;
	private final int DEFAULT_MILLER[] = {0,0,1};
	private final double DEFAULT_ALPHA = Math.PI * 2 * (90.0 / 360.0);
	private final double DEFAULT_BETA = Math.PI * 2 * (90.0 / 360.0);
	private final double DEFAULT_GAMMA = Math.PI * 2 * (90.0 / 360.0);
	private final int DEFAULT_BOX[] = {DEFAULT_SIZE, DEFAULT_SIZE, DEFAULT_SIZE};

	private LatticePlaneTestUtility lptu = null;

	public BLCPrimitiveTriclinicLatticePlaneTest(String name) {
		super(name);
		funcName = name;
	}

	protected void setUp() throws Exception {
		super.setUp();
		if (lptu == null) {
			lptu = new LatticePlaneTestUtility();			
	        lptu.createLatticeAndBox(lptu.TRICLINIC, DEFAULT_MILLER, DEFAULT_BOX);
	        lptu.setDimensions(DEFAULT_SIZE);
		}
	}

	protected void tearDown() throws Exception {
		super.tearDown();
	}

	private double[] makeArray(IVector v) {
	    return new double[] {v.x(0), v.x(1), v.x(2)};
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
    public void testStandard() {

    	int idx = 0;
    	double cubicSize = 1.0;
    	double plane = 1.0;
    	AtomSet leafList = null;

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
		    for(idx = 0; idx < leafList.getAtomCount(); idx++) {
			    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);
                if(a.getPosition().x(2) >= spacePos-epsilon &&
                   a.getPosition().x(2) <= spacePos+epsilon) {
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);
            if(a.getPosition().x(2) >= spacePos-epsilon &&
               a.getPosition().x(2) <= spacePos+epsilon) {
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
    public void testCellSizeAllDifferent() {

    	int idx = 0;
    	double cubicSizeA = 1.5;
    	double cubicSizeB = 1.75;
    	double cubicSizeC = 2.0;
    	double plane = -2.0;
    	AtomSet leafList = null;

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
		    for(idx = 0; idx < leafList.getAtomCount(); idx++) {
			    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);

                if(a.getPosition().x(2) >= spacePos-epsilon &&
                   a.getPosition().x(2) <= spacePos+epsilon) {
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);
            if(a.getPosition().x(2) >= spacePos-epsilon &&
               a.getPosition().x(2) <= spacePos+epsilon) {
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
    public void testCellSizeAllDifferentAnglesAllDifferent() {

    	int idx = 0;
    	double cubicSizeA = 1.25;
    	double cubicSizeB = 1.5;
    	double cubicSizeC = 1.75;
    	double plane = 3.0;
    	double alpha = Math.PI * 2.0 * (100.0 / 360.0);
    	double beta = Math.PI * 2.0 * (110.0 / 360.0);
    	double gamma = Math.PI * 2.0 * (120.0 / 360.0);
    	AtomSet leafList = null;

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
		    for(idx = 0; idx < leafList.getAtomCount(); idx++) {
			    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);

                if(a.getPosition().x(2) >= spacePos-epsilon &&
                   a.getPosition().x(2) <= spacePos+epsilon) {
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);
            if(a.getPosition().x(2) >= spacePos-epsilon &&
               a.getPosition().x(2) <= spacePos+epsilon) {
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
    public void testOddMillerIndices() {

    	int idx = 0;
    	double cubicSizeA = 1.1;
    	double cubicSizeB = 1.2;
    	double cubicSizeC = 1.3;
    	AtomSet leafList = null;
    	int size = 7;
    	double plane = 6.0;
    	int itemsFound = 0;
    	double alpha = Math.PI * 2.0 * (100.0 / 360.0);
    	double beta = Math.PI * 2.0 * (100.0 / 360.0);
    	double gamma = Math.PI * 2.0 * (100.0 / 360.0);
    	int[] millerIndices = new int[] { 2, 0, 2 };
        double actualPlane[][] =
           { { -0.052094453300079335, -4.3523966852137566, 3.754992688589829 },
             { -0.2604722665003958, -3.1706273815991066, 3.754992688589829 },
             { -0.46885007970071224, -1.9888580779844571, 3.754992688589829 },
             { -0.6772278929010282, -0.8070887743698072, 3.754992688589829 },
             { -0.8856057061013447, 0.37468052924484274, 3.754992688589829 },
             { -1.0939835193016612, 1.5564498328594918, 3.754992688589829 },
             { -1.3023613325019776, 2.7382191364741417, 3.754992688589829 },
             { 1.2736481776669302, -4.083367093757154, 2.5033284590598863 },
             { 1.0652703644666137, -2.901597790142504, 2.5033284590598863 },
             { 0.8568925512662973, -1.7198284865278544, 2.5033284590598863 },
             { 0.6485147380659808, -0.5380591829132046, 2.5033284590598863 },
             { 0.4401369248656648, 0.6437101207014453, 2.5033284590598863 },
             { 0.23175911166534835, 1.8254794243160943, 2.5033284590598863 },
             { 0.023381298465031897, 3.0072487279307443, 2.5033284590598863 },
             { 2.599390808633939, -3.8143375023005515, 1.251664229529943 },
             { 2.3910129954336226, -2.6325681986859015, 1.251664229529943 },
             { 2.182635182233306, -1.450798895071252, 1.251664229529943 },
             { 1.9742573690329897, -0.26902959145660255, 1.251664229529943 },
             { 1.7658795558326732, 0.9127397121580478, 1.251664229529943 },
             { 1.5575017426323567, 2.094509015772697, 1.251664229529943 },
             { 1.3491239294320403, 3.2762783193873477, 1.251664229529943 },
             { 3.925133439600949, -3.545307910843949, 0.0 },
             { 3.7167556264006327, -2.3635386072292994, 0.0 },
             { 3.5083778132003163, -1.1817693036146495, 0.0 },
             { 3.3, 0.0, 0.0 },
             { 3.0916221867996834, 1.18176930361465, 0.0 },
             { 2.883244373599368, 2.3635386072292985, 0.0 },
             { 2.6748665603990514, 3.5453079108439502, 0.0 } };

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
		    for(idx = 0; idx < leafList.getAtomCount(); idx++) {
			    IAtomPositioned a = (IAtomPositioned)leafList.getAtom(idx);

			    if(dd.contains(makeArray(a.getPosition())) == true) {
			    	itemsFound++;
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);
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

        assertEquals(actualPlane.length, itemsFound);

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
    public void testPlaneMinusFiveHundreths() {

    	int idx = 0;
    	double cubicSizeA = 1.1;
    	double cubicSizeB = 1.2;
    	double cubicSizeC = 1.3;
    	double alpha = Math.PI * 2.0 * (97.0 / 360.0);
    	double beta = Math.PI * 2.0 * (104.0 / 360.0);
    	double gamma = Math.PI * 2.0 * (115.0 / 360.0);
    	AtomSet leafList = null;
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
		    for(idx = 0; idx < leafList.getAtomCount(); idx++) {
			    IAtomPositioned a = (IAtomPositioned)leafList.getAtom(idx);

            	assertFalse(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);
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
    public void testPlane() {

    	int idx = 0;
    	double cubicSizeA = 1.1;
    	double cubicSizeB = 1.2;
    	double cubicSizeC = 1.3;
    	AtomSet leafList = null;
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
		    for(idx = 0; idx < leafList.getAtomCount(); idx++) {
			    IAtomPositioned a = (IAtomPositioned)leafList.getAtom(idx);

			    if(dd.contains(makeArray(a.getPosition())) == true) {
			    	itemsFound++;
            	    assertTrue(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
                }
                else {
            	    assertFalse(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
                }
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);
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

        assertEquals(actualPlane.length, itemsFound);

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
    public void testPlanePlusFiveHundreths() {

    	int idx = 0;
    	double cubicSizeA = 1.1;
    	double cubicSizeB = 1.2;
    	double cubicSizeC = 1.3;
    	AtomSet leafList = null;
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
		    for(idx = 0; idx < leafList.getAtomCount(); idx++) {
			    IAtomPositioned a = (IAtomPositioned)leafList.getAtom(idx);

            	assertFalse(lptu.getLatticePlane().inPlane(
            	    		(etomica.space3d.Vector3D)(a.getPosition())));
		    }
		}
        catch (junit.framework.AssertionFailedError e) {
		    IAtomPositioned a = (IAtomPositioned) leafList.getAtom(idx);
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
