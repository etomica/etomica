package simulate;

import java.util.*;
public interface IntegrationIntervalListener extends java.util.EventListener
{
    public abstract void integrationIntervalAction(IntegrationIntervalEvent evt);
}