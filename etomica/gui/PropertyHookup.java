package etomica.gui;

/**
 * This class manages hookups between properties, so that a
 * bound property change on object X turns into a property
 * set on a related property on object Y.
 * <P>
 * We do this by associating a PropertyHookup adaptor with each
 * source object that we are interested in.  As part of the adaptor
 * we keep track of which target setter methods to call when a given
 * property changes.
 */

import java.lang.reflect.*;
import java.beans.*;
import java.io.*;
import java.util.Hashtable;
import java.util.Vector;

public class PropertyHookup implements PropertyChangeListener, Serializable {

    static final long serialVersionUID = 4502052857914084293L;

    /**
     * Create a property hookup, so that a change to the named bound
     * property on the source object turns into a call on the "setter"
     * method of the given target object.
     */

    public void attach(Object source,
		       String propertyName, Method getter,
		       Object targetObject, Method setter) {

	Vector targets = (Vector) targetsByPropertyName.get(propertyName);
	if (targets == null) {
	    targets = new Vector();
	    targetsByPropertyName.put(propertyName, targets);
	}
	PropertyHookupTarget target;
	for (int i = 0; i < targets.size(); i++) {
	    target = (PropertyHookupTarget) targets.elementAt(i);
	    if (target.setter == setter && target.object == targetObject) {
		// We've already got this hookup.  Just return.
		return;
	    }
	}
	targets.addElement(new PropertyHookupTarget(targetObject,setter));

	// propagate the initial value.
	try {
	    Object args1[] = { };
	    Object value = getter.invoke(source, args1);
	    Object args2[] = { value };
	    setter.invoke(targetObject, args2);
	} catch (InvocationTargetException ex) {
	    System.err.println("Property propagation failed");
	    ex.getTargetException().printStackTrace();
        } catch (Exception ex) {
	    System.err.println("Property propagation failed");
	    ex.printStackTrace();
	}
    }

    /**
     * Version of the above for when the objects have been initialized somehow
     */

    public void attach(String propertyName,
		       Object targetObject, Method setter) {

	Vector targets = (Vector) targetsByPropertyName.get(propertyName);
	if (targets == null) {
	    targets = new Vector();
	    targetsByPropertyName.put(propertyName, targets);
	}
	PropertyHookupTarget target;
	for (int i = 0; i < targets.size(); i++) {
	    target = (PropertyHookupTarget) targets.elementAt(i);
	    if (target.setter == setter && target.object == targetObject) {
		// We've already got this hookup.  Just return.
		return;
	    }
	}
	targets.addElement(new PropertyHookupTarget(targetObject,setter));
    }

    /**
     * Constructor for a new property hookup adaptor.
     */

    public PropertyHookup(Object source) {
	this.source = source;	
	targetsByPropertyName = new Hashtable();
    }

    /**
     * This is the method that gets called when a bound property
     * changes on the source object.
     * We map the property name to a list of targets and then
     * call each of the target "setter" methods.
     */

    synchronized public void propertyChange(PropertyChangeEvent evt) {
	String propertyName = evt.getPropertyName();
	Vector targets = (Vector) targetsByPropertyName.get(propertyName);
	if (targets == null) {
	    return;
	}
	Object args[] = { evt.getNewValue() };
	for (int i = 0; i < targets.size(); i++) {
	    PropertyHookupTarget target
		= (PropertyHookupTarget)targets.elementAt(i);
	    try {
	        target.setter.invoke(target.object, args);
	    } catch (InvocationTargetException ex) {
	 	System.err.println("Property set failed");
		ex.getTargetException().printStackTrace();
	    } catch (Exception ex) {
	 	System.err.println("Unexpected Property set exception");
		ex.printStackTrace();
	    }
	}
    }

    public void vetoablePropertyChange(PropertyChangeEvent evt)
					throws PropertyVetoException {
	propertyChange(evt);
    }

    // Access to information for code generation
    // Target data
    public Hashtable getTargetsByProperty() {
	return targetsByPropertyName;
    }

    // Same but  for PropertyHookupTarget.
    // We do it here because  I don't want to make that class public with more thought.
    public Object getTargetObject(Object o) {
	PropertyHookupTarget x = (PropertyHookupTarget) o;
	return x.object;
    }
    public Method getSetterMethod(Object o) {
	PropertyHookupTarget x = (PropertyHookupTarget) o;
	return x.setter;
    }


    // Event source
    Object source;

    // Table that maps from property names to a vector of PropertyHookupTargets
    Hashtable targetsByPropertyName;
}


// Information for an event delivery target.
class PropertyHookupTarget implements Serializable {

    static final long serialVersionUID = -8352305996623495352L;

    private static int ourVersion = 1;
    Object object;
    Method setter;

    PropertyHookupTarget(Object object, Method setter) {
	this.object = object;
	this.setter = setter;
    }

    private void writeObject(ObjectOutputStream s)
 				   throws IOException {
	// Because Method objects aren't serializable, we 
        // serialize the name of the setter method.
	s.writeInt(ourVersion);
	s.writeObject(setter.toString());
	s.writeObject(object);
    }
 
    private void readObject(ObjectInputStream s)
  		  throws ClassNotFoundException, IOException {
	int version = s.readInt();
	String setterName = (String)s.readObject();
	object = s.readObject();

	// We do a rather expensive search for a setter method
	// matching the given setterName.
	setter = null;
	Method methods[] = object.getClass().getMethods();
	for (int i = 0; i < methods.length; i++) {
	    if (methods[i].toString().equals(setterName)) {
		setter = methods[i];
		break;
	    }
	}
	if (setter == null) {
	    throw new IOException("PropertyHookupTarget : no suitable setter" +
			"\n    " + setterName + 
			"\n    in class " + object.getClass());
	}
    }
}

