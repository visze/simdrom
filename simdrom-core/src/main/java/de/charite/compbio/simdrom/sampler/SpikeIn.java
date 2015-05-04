package de.charite.compbio.simdrom.sampler;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.Set;

import de.charite.compbio.simdrom.sampler.vcf.VCFSampler;

public class SpikeIn implements Iterator<VariantContext> {

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
		if ((mutationSampler != null && mutationSampler.hasNext()))
			mutationsVC = mutationSampler.next();
	}

	public VCFHeader getVCFHeader() {
		if (mutationSampler != null) {
			Set<VCFHeaderLine> metaData = new LinkedHashSet<VCFHeaderLine>();
			metaData.addAll(backgroundSampler.getFileHeader().getMetaDataInInputOrder());
			metaData.addAll(mutationSampler.getFileHeader().getMetaDataInInputOrder());
			return new VCFHeader(metaData, backgroundSampler.getSampleNames());

		} else
			return new VCFHeader(backgroundSampler.getFileHeader().getMetaDataInInputOrder(),
					backgroundSampler.getSampleNames());
	}

	@Override
	public boolean hasNext() {
		return (backgroundVC != null || mutationsVC != null);
	}

	@Override
	public VariantContext next() {
		return getNextVariantContext();
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
			else
				mutationsVC = mutationSampler.next();
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
