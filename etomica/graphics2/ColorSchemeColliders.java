package etomica.graphics2;

import etomica.api.IAtom;
import etomica.integrator.IntegratorHard;

/**
 * This colorScheme acts to color differently the two atoms that are scheduled to collide next.
 * Highlight colors are specified by the colliderColor and partnerColor fields; all other
 * atoms are colored with the baseColor.  Applies only to with a hard-potential MD integrator.
 */
public class ColorSchemeColliders implements ColorScheme {
    
    public ColorSchemeColliders(IntegratorHard integrator) 
    {
        super();
        this.integrator = integrator;
        colorsProvided = false;
        setDefaultColors();
    }
    /** Set colors to the default colors in the device - red, blue and gray */
    public void setDefaultColors( )
    {
        setColors( Color.RED, Color.BLUE, Color.GRAY50 );    
    }
    /** Set arbitrary colors */
    public void setColors( Color collider, Color partner, Color others )
    {
    	colliderColor = collider;
    	partnerColor = partner;
    	defaultColor = others;
    	colorsProvided = true;
    }
    /**
     * Applies the special colors to the colliding pair while coloring all other atoms with baseColor.
     */ 
    public int atomColor(IAtom a) {
    	if ( !colorsProvided )
    		throw new RuntimeException( "Colors not yet provided to ColorSchemeColliders object - use SetColor(Device) or SetColor(Color,Color,Color) before calling atomColor()");
        IntegratorHard.Agent colliderAgent = integrator.colliderAgent();
        if(colliderAgent == null) return 2;
        else if(a == colliderAgent.atom) return 0;
        else if(a == colliderAgent.collisionPartner) return 1;
        else return 2;
    }
    /** Color applied to the downList atom of the colliding pair */
    protected Color colliderColor;
    /**Color applied to the upList atom of the colliding pair */
    protected Color partnerColor;
    /**Color applied to the atomns that do not have a colliderAgent */
    protected Color defaultColor;
    /**The integrator that has the collision information */
    IntegratorHard integrator;
    /**Indicate whether colors were provided - one of SetColors() functions must be called */
    protected boolean colorsProvided = false;
	/* (non-Javadoc)
	 * @see etomica.graphics2.ColorScheme#getNumColors()
	 */
	public int getNumColors() {
		// TODO Auto-generated method stub
		return 3;
	}

	/* (non-Javadoc)
	 * @see etomica.graphics2.ColorScheme#getColor(int)
	 */
	public Color getColor(int index) {
		if ( index==0 ) return colliderColor;
		if ( index==1 ) return partnerColor;
		return defaultColor;
	}
    
}

