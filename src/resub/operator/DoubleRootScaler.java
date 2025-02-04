package resub.operator;

import java.text.DecimalFormat;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.inference.util.InputUtil;

@Description("Scales two tree roots at the same time")
public class DoubleRootScaler extends Operator {
	
	
	public final Input<Tree> tree1Input = new Input<>("tree1", "tree 1", Input.Validate.REQUIRED);
	public final Input<Tree> tree2Input = new Input<>("tree2", "tree 2", Input.Validate.REQUIRED);
	
	public final Input<Double> scaleFactorInput = new Input<>("scaleFactor", "scaling factor: range from 0 to 1. Close to zero is very large jumps, close to 1.0 is very small jumps.", 0.75);
	public final Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);

	
	public final Input<KernelDistribution> kernelDistributionInput = new Input<>("kernelDistribution", "provides sample distribution for proposals", 
    		KernelDistribution.newDefaultKernelDistribution());
    protected KernelDistribution kernelDistribution;
	 
	double scaleFactor;
	double lower = 1e-8;
	double upper = 1 - lower;
	
	
	@Override
	public void initAndValidate() {
		scaleFactor = scaleFactorInput.get();
		kernelDistribution = kernelDistributionInput.get();
	}
	@Override
	public double proposal() {
		
		
		
		// Get the two trees
		final Tree tree1 = (Tree) InputUtil.get(tree1Input, this);
		final Tree tree2 = (Tree) InputUtil.get(tree2Input, this);
		final Node root1 = tree1.getRoot();
		final Node root2 = tree2.getRoot();
		
		
		// Scale them both
		double maxHeight = Math.max(root1.getHeight(), root2.getHeight());
		final double scale = getScaler(0, maxHeight);
		final double newHeight1 = root1.getHeight() * scale;
		final double newHeight2 = root2.getHeight() * scale;

        if (newHeight1 < Math.max(root1.getLeft().getHeight(), root1.getRight().getHeight())) {
            return Double.NEGATIVE_INFINITY;
        }
        if (newHeight2 < Math.max(root2.getLeft().getHeight(), root2.getRight().getHeight())) {
            return Double.NEGATIVE_INFINITY;
        }
        
        
        root1.setHeight(newHeight1);
        root2.setHeight(newHeight2);
        return Math.log(scale) * 2;
		 
	}
	
	
	
	protected double getScaler(int i, double value) {
		return kernelDistribution.getScaler(i, value, getCoercableParameterValue());
	}
	
	
	 /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
        if (optimiseInput.get()) {
            double delta = calcDelta(logAlpha);
            delta += Math.log(1.0 / scaleFactor - 1.0);
            setCoercableParameterValue(1.0 / (Math.exp(delta) + 1.0));
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return scaleFactor;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        scaleFactor = Math.max(Math.min(value, upper), lower);
    }

    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        final double sf = Math.pow(scaleFactor, ratio);

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > 0.40) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }

}
