package resub.operator;



import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.StateNode;
import beast.base.inference.operator.kernel.KernelOperator;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;


@Description("Same as BactrianIntervalOperator but it allows the parameter range to vary")
public class VariableRangeBactrianIntervalOperator extends KernelOperator {
     final public Input<RealParameter> parameterInput = new Input<>("parameter", "the parameter to operate a random walk on.", Validate.REQUIRED);
    public final Input<Double> scaleFactorInput = new Input<>("scaleFactor", "scaling factor: larger means more bold proposals", 1.0);
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);
    final public Input<Boolean> inclusiveInput = new Input<>("inclusive", "are the upper and lower limits inclusive i.e. should limit values be accepted (default true)", true);
    final public Input<Function> lowerInput = new Input<>("lower", "lower limit of the parameter (default 0)");
    final public Input<Function> upperInput = new Input<>("upper", "upper limit of the parameter (default 1)");
    
    
    double scaleFactor;
    boolean inclusive;

    @Override
	public void initAndValidate() {
    	super.initAndValidate();
        scaleFactor = scaleFactorInput.get();
        inclusive = inclusiveInput.get();
    }

    @Override
    public double proposal() {

        RealParameter param = (RealParameter) InputUtil.get(parameterInput, this);
        
        double lower = lowerInput.get() == null ? 0 : lowerInput.get().getArrayValue();
        double upper = upperInput.get() == null ? 1 : upperInput.get().getArrayValue();

        int i = Randomizer.nextInt(param.getDimension());
        double value = param.getValue(i);
        double scale = kernelDistribution.getScaler(i, value, scaleFactor);
        
        
        if (value < lower || value > upper) {
        	return Double.NEGATIVE_INFINITY;
        }
        
        // transform value
        double y = (upper - value) / (value - lower);
        y *= scale;
        double newValue = (upper + lower * y) / (y + 1.0);
        
        if (newValue < lower || newValue > upper) {
        	throw new RuntimeException("programmer error: new value proposed outside range");
        }
        
        // Ensure that the value is not sitting on the limit (due to numerical issues for example)
        if (!inclusive && (newValue == lower || newValue == upper)) return Double.NEGATIVE_INFINITY;
        
        param.setValue(i, newValue);

        double logHR = Math.log(scale) + 2.0 * Math.log((newValue - lower)/(value - lower));
        return logHR;
    }


    @Override
    public double getCoercableParameterValue() {
        return scaleFactor;
    }

    @Override
    public void setCoercableParameterValue(double value) {
    	scaleFactor = value;
    }

    /**
     * called after every invocation of this operator to see whether
     * a parameter can be optimised for better acceptance hence faster
     * mixing
     *
     * @param logAlpha difference in posterior between previous state & proposed state + hasting ratio
     */

    @Override
    public void optimize(double logAlpha) {
        // must be overridden by operator implementation to have an effect
    	if (optimiseInput.get()) {
	        double delta = calcDelta(logAlpha);
	        double scaleFactor = getCoercableParameterValue();
	        delta += Math.log(scaleFactor);
	        scaleFactor = Math.exp(delta);
	        setCoercableParameterValue(scaleFactor);
    	}
    }
    
    @Override
    public double getTargetAcceptanceProbability() {
    	return 0.3;
    }
    
    
    @Override
	public List<StateNode> listStateNodes() {
		List<StateNode> nodes = new ArrayList<>();
		nodes.add(parameterInput.get());
		return nodes;
	}


    @Override
    public String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        double newWindowSize = getCoercableParameterValue() * ratio;

        DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10 || prob > 0.40) {
            return "Try setting scale factor to about " + formatter.format(newWindowSize);
        } else return "";
    }
    
} // class BactrianIntervalOperator