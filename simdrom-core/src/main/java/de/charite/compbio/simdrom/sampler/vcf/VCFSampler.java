/**
 * 
 */
package de.charite.compbio.simdrom.sampler.vcf;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import com.google.common.collect.ImmutableSet;

/**
 * @author Max Schubach <max.schubach@charite.de>
 *
 */
public class VCFSampler implements Iterator<VariantContext> {

	private double probability;
	private int counts;
	private int selected = 0;
	private String afIdentifier;
	private String sample = null;
	private VCFFileReader parser;
	private CloseableIterator<VariantContext> iterator;
	private Random random;

	public VCFSampler(String filePath) {
		this.parser = new VCFFileReader(new File(filePath), false);
	}

	public CloseableIterator<VariantContext> getIterator() {
		if (this.iterator == null)
			this.iterator = this.parser.iterator();
		return iterator;
	}

	@Override
	public boolean hasNext() {
		// FIXME has next can be true, but next can gi9ve back null!
		return getIterator().hasNext();
	}

	@Override
	public VariantContext next() {
		return getNextVariant();
	}

	private VariantContext getNextVariant() {
		VariantContext output = null;
		while (getIterator().hasNext() && output == null) {
			VariantContext candidate = getIterator().next();
			List<Integer> alleles = useAlleles(candidate);
			if (!alleles.isEmpty()) {
				output = selectOneGenotype(candidate, alleles);
				break;
			}
		}
		if (output != null)
			selected += output.getAlternateAlleles().size();
		return output;
	}

	private VariantContext selectOneGenotype(VariantContext candidate, List<Integer> alleles) {
		if (useSample()) {
			if (candidate.hasGenotype(getSample()))
				return new VariantContextBuilder(candidate).alleles(getAlleles(candidate, alleles))
						.genotypes(candidate.getGenotypes(getSample())).make();
			else
				return null;
		} else {
			return new VariantContextBuilder(candidate).genotypes(createGenotype(candidate.getAlleles())).make();
		}
	}

	private Collection<Allele> getAlleles(VariantContext candidate, List<Integer> numAlleles) {
		Collection<Allele> alleles = new ArrayList<Allele>();
		alleles.add(candidate.getReference());
		for (int i : numAlleles) {
			alleles.add(candidate.getAlternateAlleles().get(i));
		}
		return alleles;
	}

	private Genotype createGenotype(List<Allele> alleles) {
		List<Allele> filteredAlleles = new ArrayList<Allele>();
		boolean het = nextDouble() <= 0.5;

		// more the one alternative allele, do not use ref!
		if (alleles.size() > 2) {
			for (int i = 1; i < alleles.size(); i++) {
				filteredAlleles.add(alleles.get(i));
			}
		} else if (het) {
			filteredAlleles = alleles;
		} else {
			filteredAlleles.add(alleles.get(1));
			filteredAlleles.add(alleles.get(1));
		}

		return GenotypeBuilder.create("Sampled", filteredAlleles);
	}

	private List<Integer> useAlleles(VariantContext candidate) {
		List<Integer> candidates = new ArrayList<Integer>();

		if (useAF()) {// AF flag
			Object af = candidate.getCommonInfo().getAttribute(getAFIdentifier());
			if (af instanceof ArrayList<?>) {
				if (((ArrayList<?>) af).get(0) instanceof String) {
					int i = 0;
					for (Object o : (ArrayList<?>) af) {
						if (nextDouble() <= Double.parseDouble((String) o)) {
							candidates.add(i);
							i++;
						}
					}
				}
			} else if (nextDouble() <= candidate.getCommonInfo().getAttributeAsDouble(getAFIdentifier(), 0.0))
				candidates.add(0);
		} else if (useCounts()) { // probability > 1
			for (int i = 0; i < candidate.getAlternateAlleles().size(); i++) {
				if (getProbability() >= counts)
					candidates.add(i);
				if (selected <= (int) Math.ceil(getProbability())
						&& nextDouble() <= Math.ceil(getProbability()) / (double) counts)
					candidates.add(i);
			}

		} else { // probability
			for (int i = 0; i < candidate.getAlternateAlleles().size(); i++) {
				if (nextDouble() <= getProbability())
					candidates.add(i);
			}
		}
		return candidates;

	}

	private boolean useCounts() {
		return probability > 1.0;
	}

	private double nextDouble() {
		if (random == null)
			random = new Random();
		return random.nextDouble();
	}

	public void setProbability(double probability) {
		this.probability = probability;
	}

	public void setSample(String sample) {
		this.sample = sample;
	}

	public String getSample() {
		return sample;
	}

	public double getProbability() {
		return probability;
	}

	public VCFHeader getFileHeader() {
		// List<VCFHeaderLine> headerLines = new ArrayList<VCFHeaderLine>();
		// headerLines.addAll();

		Set<VCFHeaderLine> set = new LinkedHashSet<VCFHeaderLine>();
		set.addAll(parser.getFileHeader().getMetaDataInInputOrder());
		set.add(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"));
		return new VCFHeader(set, getSampleNames());
	}

	public ImmutableSet<String> getSampleNames() {
		if (getSample() == null)
			return ImmutableSet.<String> builder().add("Sampled").build();
		return ImmutableSet.<String> builder().add(getSample()).build();
	}

	public void setCounts(int counts) {
		this.counts = counts;
	}

	public void setAFIdentifier(String afIdentifier) {
		this.afIdentifier = afIdentifier;
	}

	public String getAFIdentifier() {
		return afIdentifier;
	}

	private boolean useAF() {
		return getAFIdentifier() != null;
	}

	private boolean useSample() {
		return getSample() != null;
	}

	public void close() {
		parser.close();
	}

}
