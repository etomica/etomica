package simulate;
import java.awt.*;
import java.beans.*;

public class SpeciesDumbbellBeanInfo extends SimpleBeanInfo {

 private final static Class beanClass = simulate.SpeciesDumbbell.class;
 
 public SpeciesDumbbellBeanInfo()
    {
    }
   	/**
	 * Gets a BeanInfo for the superclass of this bean.
	 * @return BeanInfo[] containing this bean's superclass BeanInfo
	 */

	public BeanInfo[] getAdditionalBeanInfo()
	{
		try
		{
			BeanInfo[] bi = new BeanInfo[1];
			bi[0] = Introspector.getBeanInfo(beanClass.getSuperclass());
			return bi;
		}
		catch (IntrospectionException e)
		{
			throw new Error(e.toString());
		}
	}
   
	public BeanDescriptor getBeanDescriptor()
	{
		BeanDescriptor bd = new BeanDescriptor(beanClass);
		return bd;
	}

	/**
	 * Gets an image that may be used to visually represent this bean
	 * (in the toolbar, on a form, etc).
	 * @param iconKind the type of icon desired, one of: BeanInfo.ICON_MONO_16x16,
	 * BeanInfo.ICON_COLOR_16x16, BeanInfo.ICON_MONO_32x32, or BeanInfo.ICON_COLOR_32x32.
	 * @return an image for this bean
	 * @see BeanInfo#ICON_MONO_16x16
	 * @see BeanInfo#ICON_COLOR_16x16
	 * @see BeanInfo#ICON_MONO_32x32
	 * @see BeanInfo#ICON_COLOR_32x32
	 */
	public java.awt.Image getIcon(int nIconKind)
	{
		java.awt.Image img = null;
		if (nIconKind == BeanInfo.ICON_COLOR_16x16)
				img = loadImage("Simulation_COLOR_16x16.gif");
		return img;
	}

	/**
	 * Returns descriptions of this bean's properties.
	 */

  public PropertyDescriptor[] getPropertyDescriptors()
  {
    try{
      PropertyDescriptor pdSpeciesOrigin = new PropertyDescriptor("speciesOrigin", beanClass);
      pdSpeciesOrigin.setPropertyEditorClass(DoubleArrayEditor.class);
      PropertyDescriptor pdMass = new PropertyDescriptor("mass", beanClass);
      pdMass.setPropertyEditorClass(DoubleArrayEditor.class);
      PropertyDescriptor pdDiameter = new PropertyDescriptor("diameter", beanClass);
      pdDiameter.setPropertyEditorClass(DoubleArrayEditor.class);
      PropertyDescriptor pdL = new PropertyDescriptor("l", beanClass);
      PropertyDescriptor[] aryPrpDescriptors = {pdSpeciesOrigin,pdMass,pdDiameter,pdL,pdDiameter};
      return aryPrpDescriptors;
    }
    catch (Exception e){
      return null;
    }
  }
  
  static
  { PropertyEditorManager.registerEditor(double[].class, DoubleArrayEditor.class);
  }
}


