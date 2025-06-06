package resub.math;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BetaDistributionImpl;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;

@Description("Expansion of the Beta distribution, except it respects the upper/lower limits of the parameter " +
 " Can handle a different prior per parameter if the dimension of 'x' is equal to the dimension of the gamma parameters " +
 " If alpha / beta are set to 0, then the parameter will be skipped (ie. 0 log density) " +
 " If the parameter exceeds lower/upper, the log density is negative infinity")
public class BetaRange extends Distribution  {
	
	final public Input<RealParameter> parameterInput = new Input<>("x", "the parameter at which the density is calculated", Input.Validate.REQUIRED);
	final public Input<Function> alphaInput = new Input<>("alpha", "first shape parameter (default 1)");
	final public Input<Function> betaInput = new Input<>("beta", "the other shape parameter (default 1)");
    final public Input<Function> lowerInput = new Input<>("lower", "lower limit of the parameter (default 0)");
    final public Input<Function> upperInput = new Input<>("upper", "upper limit of the parameter (default 1)");
	
	
	
    org.apache.commons.math.distribution.BetaDistribution m_dist = new BetaDistributionImpl(1, 1);
    
	@Override
    public void initAndValidate() {
		
		//System.out.println("BETARANGE");
  
    }
	
	

	
	@SuppressWarnings("deprecation")
	@Override
	public double calculateLogP() {
		
		//System.out.println("BETARANGE");
		
		logP = 0;
		
		RealParameter param = parameterInput.get();
		for (int i = 0; i < param.getDimension(); i ++) {
			
			//System.out.println("BETARANGE " + i + ": " + logP);
			
			double val = param.getValue(i);
			double alpha = getValOrDefault(alphaInput.get(), i, 1);
			double beta = getValOrDefault(betaInput.get(), i, 1);
			double lower = getValOrDefault(lowerInput.get(), i, 0);
			double upper = getValOrDefault(upperInput.get(), i, 1);
			
			
			//System.out.println("BETARANGE " + lower + ": " + upper);
			
			// Range check
			if (alpha <= 0 || beta <= 0) continue;
			if (val <= lower || val >= upper) {
				logP = Double.NEGATIVE_INFINITY;
				return logP;
			}
			
			// Standardise into (0,1)
			double tval = (val - lower) / (upper - lower);
			
			// Get the density of the standardised vairable
			m_dist.setAlpha(alpha);
		    m_dist.setBeta(beta);
			logP += m_dist.logDensity(tval) - Math.log(upper - lower);
			
		}
		

		
		return logP;
	}
	
	
	/**
	 * Gets the value of this parameter at this index, or returns the default
	 * @param param
	 * @param index
	 * @param defaultVal
	 * @return
	 */
	private double getValOrDefault(Function param, int index, double defaultVal) {
		
		// If the parameter is null, use the default value
		if (param == null) return defaultVal;
		
		// If the parameter has 1 dimension get the 1st value
		if (param.getDimension() == 1) return param.getArrayValue();
		
		// Multiple dimensions -> return the value at 'index'
		return param.getArrayValue(index);
		
	}
	


	@Override
	public List<String> getArguments() {
		List<String> args = new ArrayList<>();
		args.add(parameterInput.get().getID());
		return args;
	}

	@Override
	public List<String> getConditions() {
		List<String> conds = new ArrayList<>();
		if (alphaInput.get() instanceof BEASTObject) conds.add(((BEASTObject)alphaInput.get()).getID());
		if (betaInput.get() instanceof BEASTObject) conds.add(((BEASTObject)betaInput.get()).getID());
		if (lowerInput.get() instanceof BEASTObject) conds.add(((BEASTObject)lowerInput.get()).getID());
		if (upperInput.get() instanceof BEASTObject) conds.add(((BEASTObject)upperInput.get()).getID());
		return conds;
	}

	@Override
	public void sample(State state, Random random) {
		
		if (this.sampledFlag) return;
		this.sampleConditions(state, random);
		this.sampledFlag = true;
		
		RealParameter param = parameterInput.get();
		for (int i = 0; i < param.getDimension(); i ++) {
		
			double alpha = getValOrDefault(alphaInput.get(), i, 1);
			double beta = getValOrDefault(betaInput.get(), i, 1);
			double lower = getValOrDefault(lowerInput.get(), i, 0);
			double upper = getValOrDefault(upperInput.get(), i, 1);
			
			// Samle from beta distribution
			m_dist.setAlpha(alpha);
		    m_dist.setBeta(beta);
		    try {
		    	
				double b = m_dist.inverseCumulativeProbability(random.nextDouble());
				double x = b*(upper-lower) + lower;
				param.setValue(i, x);
				
				
			} catch (MathException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		
		}
		
	}

	
	

}
