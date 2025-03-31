package resub.operator;


import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;


@Description("UniformOperator but on a random number of indices")
public class UniformOperatorVariable extends Operator {
	
	final public Input<Parameter<?>> parameterInput = new Input<>("parameter", "a real or integer parameter to sample individual values for", Validate.REQUIRED, Parameter.class);
	final public Input<Double> pInput = new Input<>("p", "probability of changing any given index", 0.1);
	   
	
	Parameter<?> parameter;
	
	 double lower, upper;
    int lowerIndex, upperIndex;

	@Override
    public void initAndValidate() {
        parameter = parameterInput.get();
        if (parameter instanceof RealParameter) {
            lower = (Double) parameter.getLower();
            upper = (Double) parameter.getUpper();
        } else if (parameter instanceof IntegerParameter) {
            lowerIndex = (Integer) parameter.getLower();
            upperIndex = (Integer) parameter.getUpper();
        } else {
            throw new IllegalArgumentException("parameter should be a RealParameter or IntergerParameter, not " + parameter.getClass().getName());
        }

    }

    @Override
    public double proposal() {
    	
    	int indexToChange = Randomizer.nextInt(parameter.getDimension()); // Ensure that we always change at least 1 thing
        for (int n = 0; n < parameter.getDimension(); ++n) {
        	
        	// Change each index with probability p, ensuring that at least one thing is always changed
        	if (n == indexToChange || Randomizer.nextFloat() < pInput.get()) {
        	
	            // Do not worry about duplication, does not matter
	            int index = Randomizer.nextInt(parameter.getDimension());
	
	            if (parameter instanceof IntegerParameter) {
	                int newValue = Randomizer.nextInt(upperIndex - lowerIndex + 1) + lowerIndex; // from 0 to n-1, n must > 0,
	                ((IntegerParameter) parameter).setValue(index, newValue);
	            } else {
	                double newValue = Randomizer.nextDouble() * (upper - lower) + lower;
	                ((RealParameter) parameter).setValue(index, newValue);
	            }
            
        	}

        }

        return 0.0;
    }

	
	

}
