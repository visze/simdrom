package de.charite.compbio.simdrom.sampler;

import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Optional;
import java.util.Set;

import de.charite.compbio.simdrom.sampler.vcf.VCFSampler;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class SpikeIn implements CloseableIterator<VariantContext> {

	private VCFSampler backgroundSampler;
	private Optional<VCFSampler> mutationSampler;

	private Optional<VariantContext> backgroundVC;
	private Optional<VariantContext> mutationsVC;

	private boolean log;
	private Set<VariantContext> vcLogs;

	public SpikeIn(VCFSampler backgroundSampler, boolean log) {
		this(backgroundSampler, Optional.empty(), log);
	}

	public SpikeIn(VCFSampler backgroundSampler, Optional<VCFSampler> mutationSampler, boolean log) {
		super();
		this.backgroundSampler = backgroundSampler;
		this.mutationSampler = mutationSampler;
		this.log = log;

		this.backgroundVC = getNextBackground();
		this.mutationsVC = getNextMutation();
	}

	/**
	 * 
	 */
	private Optional<VariantContext> getNextMutation() {
		if ((this.mutationSampler.isPresent() && this.mutationSampler.get().hasNext())) {
			VariantContext mutationsVC = this.mutationSampler.get().next();
			if (mutationsVC != null) {
				String sample = VCFSampler.DEFAULT_SAMPLE_NAME;
				if (this.backgroundSampler.getSample().isPresent())
					sample = this.backgroundSampler.getSample().get();
				Genotype genotype = GenotypeBuilder.create(sample,
						(this.mutationSampler.get().getSample().isPresent()
								? mutationsVC.getGenotype(this.mutationSampler.get().getSample().get()).getAlleles()
								: mutationsVC.getAlleles()));
				mutationsVC = new VariantContextBuilder(mutationsVC).genotypes(genotype).make();
			}
			return Optional.of(mutationsVC);
		}
		return Optional.empty();
	}

	/**
	 * @return
	 */
	private Optional<VariantContext> getNextBackground() {
		if (this.backgroundSampler.hasNext())
			return Optional.of(this.backgroundSampler.next());
		else
			return Optional.empty();
	}

	public VCFHeader getVCFHeader() {
		Set<VCFHeaderLine> metaData = new LinkedHashSet<VCFHeaderLine>();
		metaData.addAll(backgroundSampler.getVCFHeader().getMetaDataInInputOrder());
		// FIXME workaround for ExAC and 1000 genome data
		metaData.add(new VCFInfoHeaderLine("OLD_VARIANT", 1, VCFHeaderLineType.String,
				"Flag in 1000 genomes that is not set in the header"));
		if (mutationSampler.isPresent()) {
			metaData.addAll(mutationSampler.get().getVCFHeader().getMetaDataInInputOrder());

		}
		return new VCFHeader(metaData, backgroundSampler.getSampleNames());
	}

	@Override
	public boolean hasNext() {
		return (backgroundVC.isPresent() || mutationsVC.isPresent());
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
		boolean mutationSelection = false;
		if (backgroundVC.isPresent() && mutationsVC.isPresent()) {
			if (backgroundVC.get().getContig().equals(mutationsVC.get().getContig())) {
				if (backgroundVC.get().getStart() <= mutationsVC.get().getStart()) {
					output = backgroundVC.get();
					mutationSelection = true;
				} else
					output = mutationsVC.get();
			} else
				output = backgroundVC.get();

		} else if (mutationsVC.isPresent()) {
			output = mutationsVC.get();
			mutationSelection = true;
		} else if (backgroundVC.isPresent())
			output = backgroundVC.get();

		if (mutationSelection) {
			mutationsVC = getNextMutation();
			addLog(output);
		} else
			backgroundVC = getNextBackground();

		return output;
	}

	private void addLog(VariantContext output) {
		if (log)
			getVcLogs().add(output);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see htsjdk.samtools.util.CloseableIterator#close()
	 */
	@Override
	public void close() {
		backgroundSampler.close();
		if (mutationSampler.isPresent())
			mutationSampler.get().close();
	}

	public Set<VariantContext> getVcLogs() {
		if (vcLogs == null)
			vcLogs = new HashSet<VariantContext>();
		return vcLogs;
	}

}
