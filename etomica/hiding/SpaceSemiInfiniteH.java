package simulate;


public class SpaceSemiInfiniteH extends simulate.SpacePeriodicHSlit
{
    private final int fullDrawSize[] = new int[D];
    private final int specialCopyOrigin[] = new int[D];

	public SpaceSemiInfiniteH()
	{
	    super();
	}

/*    protected void resetCentralOrigin(int[] phaseSize) {
        Space.uEa1T_v1Mv2_(centralOrigin,0.5,phaseSize,drawSize);  // Maybe get this from space;  origin = 0.5*(phaseSize - spaceSize)
        centralOrigin[D-1] = phaseSize[D-1] - drawSize[D-1];  //D-1 index is "height"
        uEv1(fullDrawSize,drawSize);
        fullDrawSize[D-1] = phaseSize[D-1];
        uEv1(specialCopyOrigin,centralOrigin);
        specialCopyOrigin[D-1] = 0;
    }
  */  
    public int[] getDrawSize() {return fullDrawSize;}  // drawSize height is full height of phase, regardless of imageShells

    public int[] getCopyOrigin() {return specialCopyOrigin;}
    
    
    /** Returns coordinate shifts needed to draw all images that overflow into central image
    */
    // 0, 1, or 3 shifts may be returned
    public double[][] getOverflowShifts(double[] r, double distance) {
        int shiftX = 0;
        int shiftY = 0;
        if(r[0]-distance < 0.0) {shiftX = +1;}
        else if(r[0]+distance > dimensions[0]) {shiftX = -1;}
        
        if(shiftX == 0) {return shift0;}
        else { //shiftX != 0
            shift1[0][0] = shiftX*dimensions[0];
            shift1[0][1] = 0.0;
            return shift1;
        }
    }
}