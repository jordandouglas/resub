package resub.tree;

import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.MRCAPrior;
import beast.base.inference.CalculationNode;


@Description("The height of a node returned as a scalar")
public class MRCATime extends CalculationNode implements Loggable, Function {
	
	public final Input<MRCAPrior> mrcaPriorInput = new Input<>("clade", "the clade", Validate.REQUIRED);

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public int getDimension() {
		return 1;
	}

	@Override
	public double getArrayValue(int dim) {
		return mrcaPriorInput.get().getArrayValue(1);
	}

	@Override
	public void init(PrintStream out) {
		out.print("mrca." + mrcaPriorInput.get().treeInput.get().getDateType() + "(" + mrcaPriorInput.get().taxonsetInput.get().getID() + (mrcaPriorInput.get().useOriginate ? ".originate" : "") +")\t");
	}

	@Override
	public void log(long sample, PrintStream out) {
		out.print(mrcaPriorInput.get().getArrayValue(1) + "\t");
	}

	@Override
	public void close(PrintStream out) {
		
		
	}

}
