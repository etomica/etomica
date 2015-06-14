package etomica.normalmode.nptdemo.fluid;

import etomica.api.IVector;

public interface IStuff {

    public abstract void setTemperature(double newTemperature);

    public abstract IVector[] stuff();

}