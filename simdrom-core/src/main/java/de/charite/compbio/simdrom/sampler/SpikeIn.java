package de.charite.compbio.simdrom.sampler;

import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Set;

import de.charite.compbio.simdrom.sampler.vcf.VCFSampler;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class SpikeIn implements CloseableIterator<VariantContext> {

	private VCFSampler backgroundSampler;
	private VCFSampler mutationSampler;

	private VariantContext backgroundVC = null;
	private VariantContext mutationsVC = null;
	private boolean log;
	private Set<VariantContext> vcLogs;

	public SpikeIn(VCFSampler backgroundSampler, boolean log) {
		this(backgroundSampler, null, log);
	}

	public SpikeIn(VCFSampler backgroundSampler, VCFSampler mutationSampler, boolean log) {
		super();
		this.backgroundSampler = backgroundSampler;
		this.mutationSampler = mutationSampler;
		this.log = log;

		if (backgroundSampler.hasNext())
			backgroundVC = backgroundSampler.next();
		if ((mutationSampler != null && mutationSampler.hasNext())) {
			mutationsVC = mutationSampler.next();
			// FIXME make it nicer!
			if (mutationsVC != null)
				mutationsVC = new VariantContextBuilder(mutationsVC)
						.genotypes(
								GenotypeBuilder.create(backgroundSampler.getSample(),
										(mutationSampler.getSample() == null ? mutationsVC.getAlleles()
												: mutationsVC.getGenotype(mutationSampler.getSample()).getAlleles())))
						.make();
		}
	}

	public VCFHeader getVCFHeader() {
		Set<VCFHeaderLine> metaData = new LinkedHashSet<VCFHeaderLine>();
		metaData.addAll(backgroundSampler.getFileHeader().getMetaDataInInputOrder());
		// FIXME workaround for ExAC and 1000 genome data
		metaData.add(new VCFInfoHeaderLine("OLD_VARIANT", 1, VCFHeaderLineType.String,
				"Flag in 1000 genomes that is not set in the header"));
		if (mutationSampler != null) {
			metaData.addAll(mutationSampler.getFileHeader().getMetaDataInInputOrder());

		}
		return new VCFHeader(metaData, backgroundSampler.getSampleNames());
	}

	@Override
	public boolean hasNext() {
		return (backgroundVC != null || mutationsVC != null);
	}

	@Override
	public VariantContext next() {
		return getNextVariantContext();
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}

	private VariantContext getNextVariantContext() {
		VariantContext output = null;
		boolean backgroundSelection = true;
		if (backgroundVC != null || mutationsVC != null) {

			if (mutationsVC == null)
				output = backgroundVC;
			else if (backgroundVC == null) {
				output = mutationsVC;
				backgroundSelection = false;
			} else if (backgroundVC.getContig().equals(mutationsVC.getContig())) {
				if (backgroundVC.getStart() <= mutationsVC.getStart())
					output = backgroundVC;
				else {
					output = mutationsVC;
					backgroundSelection = false;
				}
			} else {
				output = backgroundVC;
			}

			if (backgroundSelection)
				backgroundVC = backgroundSampler.next();
			else {
				mutationsVC = mutationSampler.next();
				// FIXME make it nicer!
				if (mutationsVC != null)
					mutationsVC = new VariantContextBuilder(mutationsVC)
							.genotypes(GenotypeBuilder.create(backgroundSampler.getSample(),
									mutationsVC.getGenotype(mutationSampler.getSample()).getAlleles()))
							.make();
			}
		}
		if (!backgroundSelection)
			addLog(output);
		return output;
	}

	private void addLog(VariantContext output) {
		if (log)
			getVcLogs().add(output);
	}

	public void close() {
		backgroundSampler.close();
		if (mutationSampler != null)
			mutationSampler.close();
	}

	public Set<VariantContext> getVcLogs() {
		if (vcLogs == null)
			vcLogs = new HashSet<VariantContext>();
		return vcLogs;
	}

}
