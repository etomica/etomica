package simulate;

import java.awt.*;
import simulate.*;
import java.beans.*;

public abstract class View extends java.awt.Component
{
    public Display parentDisplay;
    
	public View()
	{
	}
	
    public abstract void updateView();
    
}