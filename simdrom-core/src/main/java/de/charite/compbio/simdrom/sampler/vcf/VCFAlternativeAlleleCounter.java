/**
 * 
 */
package de.charite.compbio.simdrom.sampler.vcf;

import java.util.HashSet;
import java.util.Set;

import com.google.common.collect.ImmutableSet;

import de.charite.compbio.simdrom.filter.IFilter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * @author <a href="mailto:max.schubach@charite.de">Max Schubach</a>
 *
 */
public class VCFAlternativeAlleleCounter {

	private CloseableIterator<VariantContext> iterator;
	private ImmutableSet<IFilter> filters;
	private String sample;
	private int counts = -1;

	public VCFAlternativeAlleleCounter(CloseableIterator<VariantContext> iterator, ImmutableSet<IFilter> filters,
			String sample) {
		this.iterator = iterator;
		this.filters = filters;
		this.sample = sample;
	}

	/**
	 * @return The number of possible alternative alleles using all positions. If a sample is set only the genotype of
	 *         that sample is used.
	 */
	public int getCounts() {
		if (counts < 0)
			count();
		return counts;
	}

	private void count() {
		counts = 0;
		while (iterator.hasNext()) {
			VariantContext vc = iterator.next();
			for (IFilter iFilter : filters) {
				vc = iFilter.filter(vc);
			}
			if (vc != null) {
				if (sample != null) {
					Set<Allele> alleleSet = new HashSet<>();
					alleleSet.addAll(vc.getGenotype(sample).getAlleles());
					for (Allele allele : alleleSet) {
						if (allele.isNonReference())
							counts += 1;
					}
				} else
					counts += vc.getAlternateAlleles().size();
			}

		}
		// FIXME Really close it? test this!
		iterator.close();
	}
}
