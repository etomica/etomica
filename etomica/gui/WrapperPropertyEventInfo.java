package simulate.gui;

import java.lang.reflect.Method;
import java.beans.EventSetDescriptor;

// Subclass of WrapperEventInfo for when the connection is really to
// a "bound property"

public class WrapperPropertyEventInfo extends WrapperEventInfo {

    public WrapperPropertyEventInfo(Object target,
				    String pname,
				    Method smethod) {
	super(target, "sunw.beanbox.PropertyHookup", "propertyChange");
	propertyName = pname;
	setterName = smethod.getName();
	setterTypes = initStringFromType(smethod.getParameterTypes());
    }

    public String getPropertyName() {
	return propertyName;
    }

    public String getSetterName() {
	return setterName;
    }

    public String getSetterTypes() {
	return setterTypes;
    }

    private String initStringFromType(Class[] klass) {
	StringBuffer buf = new StringBuffer();
	buf.append("new String[] {");
	for (int i=0; i<klass.length; i++) {
	    buf.append("\"").append(klass[i].getName()).append("\"");
	    if (i != klass.length-1) {
		buf.append(", ");
	    }
	}
	buf.append("}");
	return buf.toString();
    }

    private String setterName;
    private String setterTypes;
    private String propertyName;
}

